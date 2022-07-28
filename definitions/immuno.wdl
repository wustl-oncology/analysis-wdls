version 1.0


# pipelines
import "germline_exome_hla_typing.wdl" as geht
import "rnaseq_star_fusion.wdl" as rsf
import "somatic_exome.wdl" as se
# others
import "subworkflows/phase_vcf.wdl" as pv
import "subworkflows/pvacseq.wdl" as p
import "tools/extract_hla_alleles.wdl" as eha
import "tools/hla_consensus.wdl" as hc
import "tools/intersect_known_variants.wdl" as ikv
import "tools/pvacfuse.wdl" as pf
import "types.wdl"  # !UnusedImport
import "tools/optitype_dna.wdl" as od 
import "tools/phlat.wdl" as ph

#
# These structs are needed only because MiniWDL, used by some of our
# scripts, has issues with Object types. We avoid this problem by
# enforcing struct boundaries to help the parser understand what we're
# doing.
#
# The reason to use these instead of primitives in outputs is
# to encode intended output directory for each file.
#

struct Rnaseq {
  Array[File] alignments
  Array[File] stringtie_expression
  Array[File] kallisto_expression
  Array[File?] star_fusion
}

struct Qc {
  Array[File?] tumor_rna
  Array[File?] tumor_dna
  Array[File?] normal_dna
  Array[File?] concordance
}

struct Variants {
  Array[File?] mutect
  Array[File?] strelka
  Array[File?] varscan
  Array[File?] docm
}
struct Cnv {
  Array[File?] cnvkit
}
struct Sv {
  Array[File?] manta
}
struct Somatic {
  Variants variants
  Array[File] final
  Cnv cnv
  Sv sv
}

struct Germline {
  Array[File?] variants
}


workflow immuno {
  input {

    # --------- RNAseq Inputs ------------------------------------------

    File reference_annotation
    Array[SequenceData] rna_sequence
    String sample_name

    File trimming_adapters
    String trimming_adapter_trim_end
    Int trimming_adapter_min_overlap
    Int trimming_max_uncalled
    Int trimming_min_readlength

    File kallisto_index
    File gene_transcript_lookup_table
    String? strand  # [first, second, unstranded]
    File refFlat
    File ribosomal_intervals
    File star_aligner_genome_dir_zip
    File star_fusion_genome_dir_zip
    Boolean examine_coding_effect = true
    String? fusioninspector_mode  # enum [inspect validate]
    File cdna_fasta
    File agfusion_database
    Boolean agfusion_annotate_noncanonical = true

    # --------- Somatic Exome Inputs -----------------------------------

    File reference
    File reference_fai
    File reference_dict
    File reference_alt
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_0123

    String tumor_name = "tumor"
    String tumor_sample_name
    Array[SequenceData] tumor_sequence

    String normal_name = "normal"
    String normal_sample_name
    Array[SequenceData] normal_sequence

    Array[File] bqsr_known_sites
    Array[File] bqsr_known_sites_tbi
    File bait_intervals
    File target_intervals
    Int target_interval_padding = 100
    Array[LabelledFile] per_base_intervals
    Array[LabelledFile] per_target_intervals
    Array[LabelledFile] summary_intervals

    File omni_vcf
    File omni_vcf_tbi

    String picard_metric_accumulation_level
    Int qc_minimum_mapping_quality = 0
    Int qc_minimum_base_quality = 0

    Int strelka_cpu_reserved = 8
    Int scatter_count = 50

    Int varscan_strand_filter = 0
    Int varscan_min_coverage = 8
    Float varscan_min_var_freq = 0.05
    Float varscan_p_value = 0.99
    Float? varscan_max_normal_freq

    File docm_vcf
    File docm_vcf_tbi

    Boolean filter_docm_variants = true
    String? gnomad_field_name

    File vep_cache_dir_zip
    String vep_ensembl_assembly
    String vep_ensembl_version
    String vep_ensembl_species
    File? synonyms_file
    Boolean annotate_coding_only = false
    # one of [pick, flag_pick, pick-allele, per_gene, pick_allele_gene, flag_pick_allele, flag_pick_allele_gene]
    String? vep_pick
    Boolean cle_vcf_filter = false
    Array[String] vep_to_table_fields = ["HGVSc", "HGVSp"]
    Array[String] variants_to_table_genotype_fields = ["GT", "AD"]
    Array[String] variants_to_table_fields = ["CHROM", "POS", "ID", "REF", "ALT", "set", "AC", "AF"]
    Array[VepCustomAnnotation] vep_custom_annotations

    File? manta_call_regions
    File? manta_call_regions_tbi
    Boolean manta_non_wgs = true
    Boolean? manta_output_contigs

    File somalier_vcf
    File? validated_variants
    File? validated_variants_tbi

    # --------- Germline Inputs ----------------------------------------

    Array[String] gvcf_gq_bands
    Array[Array[String]] gatk_haplotypecaller_intervals
    Int? ploidy
    String? optitype_name

    # --------- Phase VCF Inputs ---------------------------------------

    Array[String]? clinical_mhc_classI_alleles
    Array[String]? clinical_mhc_classII_alleles

    # --------- PVACseq Inputs -----------------------------------------
    String hla_source_mode
    Int? readcount_minimum_base_quality
    Int? readcount_minimum_mapping_quality
    Array[String] prediction_algorithms
    Array[Int]? epitope_lengths_class_i
    Array[Int]? epitope_lengths_class_ii
    Int? binding_threshold
    Int? percentile_threshold
    Float? minimum_fold_change
    String? top_score_metric  # enum [lowest, median]
    String? additional_report_columns  # enum [sample_name]
    Int? fasta_size
    Int? downstream_sequence_length
    Boolean? exclude_nas
    Int? maximum_transcript_support_level  # enum [1 2 3 4 5]
    Int? normal_cov
    Int? tdna_cov
    Int? trna_cov
    Float? normal_vaf
    Float? tdna_vaf
    Float? trna_vaf
    Float? expn_val
    String? net_chop_method  # enum [cterm 20s]
    Float? net_chop_threshold
    Boolean? netmhc_stab
    Boolean? run_reference_proteome_similarity
    String? blastp_db  # enum [refseq_select_prot refseq_protein]
    Int? pvacseq_threads
    Int? iedb_retries
    Boolean? pvacfuse_keep_tmp_files
  }

  call rsf.rnaseqStarFusion as rna {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    reference_annotation=reference_annotation,
    unaligned=rna_sequence,
    sample_name=sample_name,
    trimming_adapters=trimming_adapters,
    trimming_adapter_trim_end=trimming_adapter_trim_end,
    trimming_adapter_min_overlap=trimming_adapter_min_overlap,
    trimming_max_uncalled=trimming_max_uncalled,
    trimming_min_readlength=trimming_min_readlength,
    kallisto_index=kallisto_index,
    gene_transcript_lookup_table=gene_transcript_lookup_table,
    strand=strand,
    refFlat=refFlat,
    ribosomal_intervals=ribosomal_intervals,
    star_genome_dir_zip=star_aligner_genome_dir_zip,
    star_fusion_genome_dir_zip=star_fusion_genome_dir_zip,
    examine_coding_effect=examine_coding_effect,
    fusioninspector_mode=fusioninspector_mode,
    cdna_fasta=cdna_fasta,
    agfusion_database=agfusion_database,
    agfusion_annotate_noncanonical=agfusion_annotate_noncanonical
  }

  call se.somaticExome {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    reference_alt=reference_alt,
    reference_amb=reference_amb,
    reference_ann=reference_ann,
    reference_bwt=reference_bwt,
    reference_pac=reference_pac,
    reference_0123=reference_0123,
    tumor_sequence=tumor_sequence,
    tumor_name=tumor_name,
    normal_sequence=normal_sequence,
    normal_name=normal_name,
    bqsr_known_sites=bqsr_known_sites,
    bqsr_known_sites_tbi=bqsr_known_sites_tbi,
    bait_intervals=bait_intervals,
    target_intervals=target_intervals,
    target_interval_padding=target_interval_padding,
    per_base_intervals=per_base_intervals,
    per_target_intervals=per_target_intervals,
    summary_intervals=summary_intervals,
    gnomad_field_name=gnomad_field_name,
    omni_vcf=omni_vcf,
    omni_vcf_tbi=omni_vcf_tbi,
    picard_metric_accumulation_level=picard_metric_accumulation_level,
    qc_minimum_mapping_quality=qc_minimum_mapping_quality,
    qc_minimum_base_quality=qc_minimum_base_quality,
    strelka_cpu_reserved=strelka_cpu_reserved,
    scatter_count=scatter_count,
    varscan_strand_filter=varscan_strand_filter,
    varscan_min_coverage=varscan_min_coverage,
    varscan_min_var_freq=varscan_min_var_freq,
    varscan_p_value=varscan_p_value,
    varscan_max_normal_freq=varscan_max_normal_freq,
    docm_vcf=docm_vcf,
    docm_vcf_tbi=docm_vcf_tbi,
    filter_docm_variants=filter_docm_variants,
    vep_cache_dir_zip=vep_cache_dir_zip,
    vep_ensembl_assembly=vep_ensembl_assembly,
    vep_ensembl_version=vep_ensembl_version,
    vep_ensembl_species=vep_ensembl_species,
    synonyms_file=synonyms_file,
    annotate_coding_only=annotate_coding_only,
    vep_pick=vep_pick,
    cle_vcf_filter=cle_vcf_filter,
    variants_to_table_fields=variants_to_table_fields,
    variants_to_table_genotype_fields=variants_to_table_genotype_fields,
    vep_to_table_fields=vep_to_table_fields,
    vep_custom_annotations=vep_custom_annotations,
    manta_call_regions=manta_call_regions,
    manta_call_regions_tbi=manta_call_regions_tbi,
    manta_non_wgs=manta_non_wgs,
    manta_output_contigs=manta_output_contigs,
    somalier_vcf=somalier_vcf,
    tumor_sample_name=tumor_sample_name,
    normal_sample_name=normal_sample_name,
    validated_variants=validated_variants,
    validated_variants_tbi=validated_variants_tbi
  }

  call geht.germlineExomeHlaTyping as germlineExome {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    reference_alt=reference_alt,
    reference_amb=reference_amb,
    reference_ann=reference_ann,
    reference_bwt=reference_bwt,
    reference_pac=reference_pac,
    reference_0123=reference_0123,
    sequence=normal_sequence,
    bqsr_known_sites=bqsr_known_sites,
    bqsr_known_sites_tbi=bqsr_known_sites_tbi,
    bait_intervals=bait_intervals,
    target_intervals=target_intervals,
    target_interval_padding=target_interval_padding,
    per_base_intervals=per_base_intervals,
    per_target_intervals=per_target_intervals,
    summary_intervals=summary_intervals,
    omni_vcf=omni_vcf,
    omni_vcf_tbi=omni_vcf_tbi,
    picard_metric_accumulation_level=picard_metric_accumulation_level,
    gvcf_gq_bands=gvcf_gq_bands,
    intervals=gatk_haplotypecaller_intervals,
    ploidy=ploidy,
    vep_cache_dir_zip=vep_cache_dir_zip,
    vep_ensembl_assembly=vep_ensembl_assembly,
    vep_ensembl_version=vep_ensembl_version,
    vep_ensembl_species=vep_ensembl_species,
    vep_custom_annotations=vep_custom_annotations,
    synonyms_file=synonyms_file,
    annotate_coding_only=annotate_coding_only,
    qc_minimum_mapping_quality=qc_minimum_mapping_quality,
    qc_minimum_base_quality=qc_minimum_base_quality,
    optitype_name=optitype_name
  }

  call od.optitypeDna as optitype {
    input:
    reference=reference,
    reference_fai=reference_fai,
    cram=somaticExome.tumor_cram,
    cram_crai=somaticExome.tumor_cram_crai,
    optitype_name=optitype_name
  }

  call ph.phlat {
    input:
    cram=somaticExome.tumor_cram,
    cram_crai=somaticExome.tumor_cram_crai,
    reference=reference,
    reference_fai=reference_fai
  } 

  call pv.phaseVcf {
    input:
    somatic_vcf=somaticExome.final_filtered_vcf,
    somatic_vcf_tbi=somaticExome.final_filtered_vcf_tbi,
    germline_vcf=germlineExome.final_vcf,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=somaticExome.tumor_cram,
    bam_bai=somaticExome.tumor_cram_crai,
    normal_sample_name=normal_sample_name,
    tumor_sample_name=tumor_sample_name
  }

  call eha.extractHlaAlleles as extractAlleles {
    input:
    file=germlineExome.optitype_tsv, 
    phlat_file=germlineExome.phlat_summary
  }

  call hc.hlaConsensus {
    input:
    hla_source_mode=hla_source_mode,
    optitype_hla_alleles=extractAlleles.allele_string,
    clinical_mhc_classI_alleles=clinical_mhc_classI_alleles,
    clinical_mhc_classII_alleles=clinical_mhc_classII_alleles
  }

  call ikv.intersectKnownVariants as intersectPassingVariants {
    input:
    vcf=somaticExome.final_filtered_vcf,
    vcf_tbi=somaticExome.final_filtered_vcf_tbi
  }

  call p.pvacseq {
    input:
    detect_variants_vcf=intersectPassingVariants.validated_and_pipeline_vcf,
    detect_variants_vcf_tbi=intersectPassingVariants.validated_and_pipeline_vcf_tbi,
    sample_name=tumor_sample_name,
    normal_sample_name=normal_sample_name,
    rnaseq_bam=rna.final_bam,
    rnaseq_bam_bai=rna.final_bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    readcount_minimum_base_quality=readcount_minimum_base_quality,
    readcount_minimum_mapping_quality=readcount_minimum_mapping_quality,
    gene_expression_file=rna.kallisto_gene_abundance,
    transcript_expression_file=rna.kallisto_transcript_abundance_tsv,
    alleles=hlaConsensus.consensus_alleles,
    prediction_algorithms=prediction_algorithms,
    epitope_lengths_class_i=epitope_lengths_class_i,
    epitope_lengths_class_ii=epitope_lengths_class_ii,
    binding_threshold=binding_threshold,
    percentile_threshold=percentile_threshold,
    minimum_fold_change=minimum_fold_change,
    top_score_metric=top_score_metric,
    additional_report_columns=additional_report_columns,
    fasta_size=fasta_size,
    downstream_sequence_length=downstream_sequence_length,
    exclude_nas=exclude_nas,
    phased_proximal_variants_vcf=phaseVcf.phased_vcf,
    phased_proximal_variants_vcf_tbi=phaseVcf.phased_vcf_tbi,
    maximum_transcript_support_level=maximum_transcript_support_level,
    normal_cov=normal_cov,
    tdna_cov=tdna_cov,
    trna_cov=trna_cov,
    normal_vaf=normal_vaf,
    tdna_vaf=tdna_vaf,
    trna_vaf=trna_vaf,
    expn_val=expn_val,
    net_chop_method=net_chop_method,
    net_chop_threshold=net_chop_threshold,
    netmhc_stab=netmhc_stab,
    run_reference_proteome_similarity=run_reference_proteome_similarity,
    n_threads=pvacseq_threads,
    variants_to_table_fields=variants_to_table_fields,
    variants_to_table_genotype_fields=variants_to_table_genotype_fields,
    vep_to_table_fields=vep_to_table_fields
  }

  call pf.pvacfuse {
    input:
    input_fusions_zip=rna.annotated_fusion_predictions_zip,
    sample_name=tumor_sample_name,
    alleles=hlaConsensus.consensus_alleles,
    prediction_algorithms=prediction_algorithms,
    epitope_lengths_class_i=epitope_lengths_class_i,
    epitope_lengths_class_ii=epitope_lengths_class_ii,
    binding_threshold=binding_threshold,
    percentile_threshold=percentile_threshold,
    iedb_retries=iedb_retries,
    keep_tmp_files=pvacfuse_keep_tmp_files,
    net_chop_method=net_chop_method,
    netmhc_stab=netmhc_stab,
    top_score_metric=top_score_metric,
    net_chop_threshold=net_chop_threshold,
    run_reference_proteome_similarity=run_reference_proteome_similarity,
    blastp_db=blastp_db,
    additional_report_columns=additional_report_columns,
    fasta_size=fasta_size,
    downstream_sequence_length=downstream_sequence_length,
    exclude_nas=exclude_nas,
    n_threads=pvacseq_threads
  }

  output {
    # ---------- RNAseq Outputs ----------------------------------------
    Rnaseq rnaseq = object {
      alignments: [
        rna.final_bam,
        rna.final_bam_bai
      ],
      stringtie_expression: [
        rna.stringtie_transcript_gtf,
        rna.stringtie_gene_expression_tsv
      ],
      kallisto_expression: [
        rna.kallisto_transcript_abundance_tsv,
        rna.kallisto_transcript_abundance_h5,
        rna.kallisto_gene_abundance,
        rna.kallisto_fusion_evidence
      ],
      star_fusion: [
        rna.star_fusion_out,
        rna.star_junction_out,
        rna.star_fusion_predict,
        rna.star_fusion_abridge,
        rna.star_fusion_coding_region_effects,
        rna.annotated_fusion_predictions_zip
      ]
    }

    # -------- Somatic Outputs -----------------------------------------

    Qc qc =  object {
      tumor_rna: flatten([
        [ rna.metrics, 
          rna.chart ],
        rna.strand_info  
      ]),
      tumor_dna: flatten([
        [ somaticExome.tumor_mark_duplicates_metrics,
          somaticExome.tumor_insert_size_metrics,
          somaticExome.tumor_alignment_summary_metrics,
          somaticExome.tumor_hs_metrics,
          somaticExome.tumor_flagstats,
          somaticExome.tumor_verify_bam_id_metrics,
          somaticExome.tumor_verify_bam_id_depth ],
        somaticExome.tumor_per_target_coverage_metrics,
        somaticExome.tumor_per_target_hs_metrics,
        somaticExome.tumor_per_base_coverage_metrics,
        somaticExome.tumor_per_base_hs_metrics,
        somaticExome.tumor_summary_hs_metrics
      ]),
      normal_dna: flatten([
        [ somaticExome.normal_mark_duplicates_metrics,
          somaticExome.normal_insert_size_metrics,
          somaticExome.normal_alignment_summary_metrics,
          somaticExome.normal_hs_metrics,
          somaticExome.normal_flagstats,
          somaticExome.normal_verify_bam_id_metrics,
          somaticExome.normal_verify_bam_id_depth ],
        somaticExome.normal_per_target_coverage_metrics,
        somaticExome.normal_per_target_hs_metrics,
        somaticExome.normal_per_base_coverage_metrics,
        somaticExome.normal_per_base_hs_metrics,
        somaticExome.normal_summary_hs_metrics
      ]),
      concordance: [
        somaticExome.somalier_concordance_metrics,
        somaticExome.somalier_concordance_statistics
      ]
    }

    File tumor_cram = somaticExome.tumor_cram
    File tumor_cram_crai = somaticExome.tumor_cram_crai
    File normal_cram = somaticExome.normal_cram
    File normal_cram_crai = somaticExome.normal_cram_crai

    Somatic somatic = object {
      variants: object {
        mutect: [
          somaticExome.mutect_unfiltered_vcf,
          somaticExome.mutect_unfiltered_vcf_tbi,
          somaticExome.mutect_filtered_vcf,
          somaticExome.mutect_filtered_vcf_tbi
        ],
        strelka: [
          somaticExome.strelka_unfiltered_vcf,
          somaticExome.strelka_unfiltered_vcf_tbi,
          somaticExome.strelka_filtered_vcf,
          somaticExome.strelka_filtered_vcf_tbi
        ],
        varscan: [
          somaticExome.varscan_unfiltered_vcf,
          somaticExome.varscan_unfiltered_vcf_tbi,
          somaticExome.varscan_filtered_vcf,
          somaticExome.varscan_filtered_vcf_tbi,
        ],
        docm: [
          somaticExome.docm_filtered_vcf,
          somaticExome.docm_filtered_vcf_tbi
        ],
      },
      final: [
        somaticExome.final_vcf,
        somaticExome.final_vcf_tbi,
        somaticExome.final_filtered_vcf,
        somaticExome.final_filtered_vcf_tbi,
        somaticExome.final_tsv
      ],
      cnv: object {cnvkit: [
        somaticExome.intervals_antitarget,
        somaticExome.intervals_target,
        somaticExome.normal_antitarget_coverage,
        somaticExome.normal_target_coverage,
        somaticExome.reference_coverage,
        somaticExome.cn_diagram,
        somaticExome.cn_scatter_plot,
        somaticExome.tumor_antitarget_coverage,
        somaticExome.tumor_target_coverage,
        somaticExome.tumor_bin_level_ratios,
        somaticExome.tumor_segmented_ratios
      ]},
      sv: object {manta: [
        somaticExome.diploid_variants,
        somaticExome.diploid_variants_tbi,
        somaticExome.somatic_variants,
        somaticExome.somatic_variants_tbi,
        somaticExome.all_candidates,
        somaticExome.all_candidates_tbi,
        somaticExome.small_candidates,
        somaticExome.small_candidates_tbi,
        somaticExome.tumor_only_variants,
        somaticExome.tumor_only_variants_tbi
      ]}
    }

    # ---------- Germline Outputs --------------------------------------

    Germline germline = object {
      variants: [
        germlineExome.final_vcf,
        germlineExome.final_vcf_tbi,
        germlineExome.filtered_vcf,
        germlineExome.filtered_vcf_tbi,
        germlineExome.vep_summary
      ]
    }

    Array[File] hla_typing_normal = flatten([
      [germlineExome.optitype_tsv,
       germlineExome.optitype_plot,
       extractAlleles.allele_file,
       hlaConsensus.consensus_alleles_file],
      hlaConsensus.hla_call_files
    ])

    Array[File] phlat_hla_typing_normal = [germlineExome.phlat_summary]

    # --------- Other Outputs ------------------------------------------

    Array[File] pvactools = flatten([
      [phaseVcf.phased_vcf,
       phaseVcf.phased_vcf_tbi],
      pvacseq.pvacseq_predictions
    ])

    File annotated_vcf = pvacseq.annotated_vcf
    File annotated_tsv = pvacseq.annotated_tsv

    Array[File] hla_typing_tumor = [optitype.optitype_tsv]
    Array[File] pvacfuse_predictions = pvacfuse.pvacfuse_predictions
    Array[File] fusioninspector_evidence = rna.fusioninspector_evidence
    Array[File] phlat_hla_typing_tumor = [phlat.phlat_summary]
  }
}
