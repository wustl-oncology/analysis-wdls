version 1.0

# pipelines
import "germline_exome_hla_typing.wdl" as geht
import "rnaseq.wdl" as r
import "somatic_exome.wdl" as se
# others
import "subworkflows/phase_vcf.wdl" as pv
import "subworkflows/pvacseq.wdl" as p
import "tools/extract_hla_alleles.wdl" as eha
import "tools/hla_consensus.wdl" as hc
import "tools/intersect_known_variants.wdl" as ikv
import "types.wdl"

workflow immuno {
  input {

    # --------- RNAseq Inputs ------------------------------------------

    File reference_index
    File reference_index_1ht2
    File reference_index_2ht2
    File reference_index_3ht2
    File reference_index_4ht2
    File reference_index_5ht2
    File reference_index_6ht2
    File reference_index_7ht2
    File reference_index_8ht2

    File reference_annotation
    Array[SequenceData] rna_sequence
    Array[String] rna_readgroups
    Array[Array[String]] read_group_fields
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
    File? ribosomal_intervals

    # --------- Somatic Exome Inputs -----------------------------------

    File reference
    File reference_fai
    File reference_dict
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_sa

    String tumor_name = "tumor"
    String tumor_sample_name
    Array[SequenceData] tumor_sequence

    String normal_name = "normal"
    String normal_sample_name
    Array[SequenceData] normal_sequence

    Array[File] bqsr_known_sites
    Array[File] bqsr_known_sites_tbi

    Array[String] bqsr_intervals
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
    Float varscan_min_var_freq = 0.1
    Float varscan_p_value = 0.99
    Float? varscan_max_normal_freq

    Int pindel_insert_size = 400

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
    File? phased_proximal_variants_vcf
    File? phased_proximal_variants_vcf_tbi
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
    Int? pvacseq_threads
  }

  call r.rnaseq {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    reference_index=reference_index,
    reference_index_1ht2=reference_index_1ht2,
    reference_index_2ht2=reference_index_2ht2,
    reference_index_3ht2=reference_index_3ht2,
    reference_index_4ht2=reference_index_4ht2,
    reference_index_5ht2=reference_index_5ht2,
    reference_index_6ht2=reference_index_6ht2,
    reference_index_7ht2=reference_index_7ht2,
    reference_index_8ht2=reference_index_8ht2,
    reference_annotation=reference_annotation,
    rna_sequence=rna_sequence,
    read_group_id=rna_readgroups,
    read_group_fields=read_group_fields,
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
    ribosomal_intervals=ribosomal_intervals
  }

  call se.somaticExome as somatic {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    reference_amb=reference_amb,
    reference_ann=reference_ann,
    reference_bwt=reference_bwt,
    reference_pac=reference_pac,
    reference_sa=reference_sa,
    tumor_sequence=tumor_sequence,
    tumor_name=tumor_name,
    normal_sequence=normal_sequence,
    normal_name=normal_name,
    bqsr_known_sites=bqsr_known_sites,
    bqsr_known_sites_tbi=bqsr_known_sites_tbi,
    bqsr_intervals=bqsr_intervals,
    bait_intervals=bait_intervals,
    target_intervals=target_intervals,
    target_interval_padding=target_interval_padding,
    per_base_intervals=per_base_intervals,
    per_target_intervals=per_target_intervals,
    summary_intervals=summary_intervals,
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
    pindel_insert_size=pindel_insert_size,
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

  call geht.germlineExomeHlaTyping as germline {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    reference_amb=reference_amb,
    reference_ann=reference_ann,
    reference_bwt=reference_bwt,
    reference_pac=reference_pac,
    reference_sa=reference_sa,
    sequence=normal_sequence,
    bqsr_known_sites=bqsr_known_sites,
    bqsr_known_sites_tbi=bqsr_known_sites_tbi,
    bqsr_intervals=bqsr_intervals,
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

  call pv.phaseVcf {
    input:
    somatic_vcf=somatic.final_filtered_vcf,
    somatic_vcf_tbi=somatic.final_filtered_vcf_tbi,
    germline_vcf=germline.final_vcf,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=somatic.tumor_cram,
    bam_bai=somatic.tumor_cram_crai,
    normal_sample_name=normal_sample_name,
    tumor_sample_name=tumor_sample_name
  }

  call eha.extractHlaAlleles as extractAlleles {
    input: allele_file=germline.optitype_tsv
  }

  call hc.hlaConsensus {
    input:
    optitype_hla_alleles=extractAlleles.allele_string,
    clinical_mhc_classI_alleles=clinical_mhc_classI_alleles,
    clinical_mhc_classII_alleles=clinical_mhc_classII_alleles
  }

  call ikv.intersectKnownVariants as intersectPassingVariants {
    input:
    vcf=somatic.final_filtered_vcf,
    vcf_tbi=somatic.final_filtered_vcf_tbi
  }

  call p.pvacseq {
    input:
    detect_variants_vcf=intersectPassingVariants.validated_and_pipeline_vcf,
    detect_variants_vcf_tbi=intersectPassingVariants.validated_and_pipeline_vcf_tbi,
    sample_name=tumor_sample_name,
    normal_sample_name=normal_sample_name,
    rnaseq_bam=rnaseq.final_bam,
    rnaseq_bam_bai=rnaseq.final_bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    readcount_minimum_base_quality=readcount_minimum_base_quality,
    readcount_minimum_mapping_quality=readcount_minimum_mapping_quality,
    gene_expression_file=rnaseq.gene_abundance,
    transcript_expression_file=rnaseq.transcript_abundance_tsv,
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

  output {
    # ---------- RNAseq Outputs ----------------------------------------

    File final_bigwig = rnaseq.bamcoverage_bigwig
    File final_bam = rnaseq.final_bam
    File final_bam_bai = rnaseq.final_bam_bai
    File stringtie_transcript_gtf = rnaseq.stringtie_transcript_gtf
    File stringtie_gene_expression_tsv = rnaseq.stringtie_gene_expression_tsv
    File transcript_abundance_tsv = rnaseq.transcript_abundance_tsv
    File transcript_abundance_h5 = rnaseq.transcript_abundance_h5
    File gene_abundance = rnaseq.gene_abundance
    File metrics = rnaseq.metrics
    File? chart = rnaseq.chart

    # -------- Somatic Outputs -----------------------------------------

    File tumor_cram = somatic.tumor_cram
    File tumor_mark_duplicates_metrics = somatic.tumor_mark_duplicates_metrics
    File tumor_insert_size_metrics = somatic.tumor_insert_size_metrics
    File tumor_alignment_summary_metrics = somatic.tumor_alignment_summary_metrics
    File tumor_hs_metrics = somatic.tumor_hs_metrics
    Array[File] tumor_per_target_coverage_metrics = somatic.tumor_per_target_coverage_metrics
    Array[File] tumor_per_target_hs_metrics = somatic.tumor_per_target_hs_metrics
    Array[File] tumor_per_base_coverage_metrics = somatic.tumor_per_base_coverage_metrics
    Array[File] tumor_per_base_hs_metrics = somatic.tumor_per_base_hs_metrics
    Array[File] tumor_summary_hs_metrics = somatic.tumor_summary_hs_metrics
    File tumor_flagstats = somatic.tumor_flagstats
    File tumor_verify_bam_id_metrics = somatic.tumor_verify_bam_id_metrics
    File tumor_verify_bam_id_depth = somatic.tumor_verify_bam_id_depth

    File normal_cram = somatic.normal_cram
    File normal_mark_duplicates_metrics = somatic.normal_mark_duplicates_metrics
    File normal_insert_size_metrics = somatic.normal_insert_size_metrics
    File normal_alignment_summary_metrics = somatic.normal_alignment_summary_metrics
    File normal_hs_metrics = somatic.normal_hs_metrics
    Array[File] normal_per_target_coverage_metrics = somatic.normal_per_target_coverage_metrics
    Array[File] normal_per_target_hs_metrics = somatic.normal_per_target_hs_metrics
    Array[File] normal_per_base_coverage_metrics = somatic.normal_per_base_coverage_metrics
    Array[File] normal_per_base_hs_metrics = somatic.normal_per_base_hs_metrics
    Array[File] normal_summary_hs_metrics = somatic.normal_summary_hs_metrics
    File normal_flagstats = somatic.normal_flagstats
    File normal_verify_bam_id_metrics = somatic.normal_verify_bam_id_metrics
    File normal_verify_bam_id_depth = somatic.normal_verify_bam_id_depth

    File mutect_unfiltered_vcf = somatic.mutect_unfiltered_vcf
    File mutect_unfiltered_vcf_tbi = somatic.mutect_unfiltered_vcf_tbi
    File mutect_filtered_vcf = somatic.mutect_filtered_vcf
    File mutect_filtered_vcf_tbi = somatic.mutect_filtered_vcf_tbi

    File strelka_unfiltered_vcf = somatic.strelka_unfiltered_vcf
    File strelka_unfiltered_vcf_tbi = somatic.strelka_unfiltered_vcf_tbi
    File strelka_filtered_vcf = somatic.strelka_filtered_vcf
    File strelka_filtered_vcf_tbi = somatic.strelka_filtered_vcf_tbi

    File varscan_unfiltered_vcf = somatic.varscan_unfiltered_vcf
    File varscan_unfiltered_vcf_tbi = somatic.varscan_unfiltered_vcf_tbi
    File varscan_filtered_vcf = somatic.varscan_filtered_vcf
    File varscan_filtered_vcf_tbi = somatic.varscan_filtered_vcf_tbi

    File pindel_unfiltered_vcf = somatic.pindel_unfiltered_vcf
    File pindel_unfiltered_vcf_tbi = somatic.pindel_unfiltered_vcf_tbi
    File pindel_filtered_vcf = somatic.pindel_filtered_vcf
    File pindel_filtered_vcf_tbi = somatic.pindel_filtered_vcf_tbi

    File docm_filtered_vcf = somatic.docm_filtered_vcf
    File docm_filtered_vcf_tbi = somatic.docm_filtered_vcf_tbi

    File somatic_final_vcf = somatic.final_vcf
    File somatic_final_vcf_tbi = somatic.final_vcf_tbi
    File final_filtered_vcf = somatic.final_filtered_vcf
    File final_filtered_vcf_tbi = somatic.final_filtered_vcf_tbi

    File final_tsv = somatic.final_tsv
    File somatic_vep_summary = somatic.vep_summary
    File tumor_snv_bam_readcount_tsv = somatic.tumor_snv_bam_readcount_tsv
    File tumor_indel_bam_readcount_tsv = somatic.tumor_indel_bam_readcount_tsv
    File normal_snv_bam_readcount_tsv = somatic.normal_snv_bam_readcount_tsv
    File normal_indel_bam_readcount_tsv = somatic.normal_indel_bam_readcount_tsv

    File? intervals_antitarget = somatic.intervals_antitarget
    File? intervals_target = somatic.intervals_target
    File? normal_antitarget_coverage = somatic.normal_antitarget_coverage
    File? normal_target_coverage = somatic.normal_target_coverage
    File? reference_coverage = somatic.reference_coverage

    File? cn_diagram = somatic.cn_diagram
    File? cn_scatter_plot = somatic.cn_scatter_plot

    File tumor_antitarget_coverage = somatic.tumor_antitarget_coverage
    File tumor_target_coverage = somatic.tumor_target_coverage
    File tumor_bin_level_ratios = somatic.tumor_bin_level_ratios
    File tumor_segmented_ratios = somatic.tumor_segmented_ratios

    File? diploid_variants = somatic.diploid_variants
    File? diploid_variants_tbi = somatic.diploid_variants_tbi
    File? somatic_variants = somatic.somatic_variants
    File? somatic_variants_tbi = somatic.somatic_variants_tbi
    File all_candidates = somatic.all_candidates
    File all_candidates_tbi = somatic.all_candidates_tbi
    File small_candidates = somatic.small_candidates
    File small_candidates_tbi = somatic.small_candidates_tbi
    File? tumor_only_variants = somatic.tumor_only_variants
    File? tumor_only_variants_tbi = somatic.tumor_only_variants_tbi

    File somalier_concordance_metrics = somatic.somalier_concordance_metrics
    File somalier_concordance_statistics = somatic.somalier_concordance_statistics

    # ---------- Germline Outputs --------------------------------------

    File cram = germline.cram
    File mark_duplicates_metrics = germline.mark_duplicates_metrics
    File insert_size_metrics = germline.insert_size_metrics
    File insert_size_histogram = germline.insert_size_histogram
    File alignment_summary_metrics = germline.alignment_summary_metrics
    File hs_metrics = germline.hs_metrics
    Array[File] per_target_coverage_metrics = germline.per_target_coverage_metrics
    Array[File] per_target_hs_metrics = germline.per_target_hs_metrics
    Array[File] per_base_coverage_metrics = germline.per_base_coverage_metrics
    Array[File] per_base_hs_metrics = germline.per_base_hs_metrics
    Array[File] summary_hs_metrics = germline.summary_hs_metrics
    File flagstats = germline.flagstats
    File verify_bam_id_metrics = germline.verify_bam_id_metrics
    File verify_bam_id_depth = germline.verify_bam_id_depth
    File germline_raw_vcf = germline.raw_vcf
    File germline_raw_vcf_tbi = germline.raw_vcf_tbi
    File germline_final_vcf = germline.final_vcf
    File germline_final_vcf_tbi = germline.final_vcf_tbi
    File germline_filtered_vcf = germline.filtered_vcf
    File germline_filtered_vcf_tbi = germline.filtered_vcf_tbi
    File germline_vep_summary = germline.vep_summary
    File optitype_tsv = germline.optitype_tsv
    File optitype_plot = germline.optitype_plot

    # --------- Other Outputs ------------------------------------------

    File phased_vcf = phaseVcf.phased_vcf
    File phased_vcf_tbi = phaseVcf.phased_vcf_tbi
    Array[String] allele_string = extractAlleles.allele_string
    Array[String] consensus_alleles = hlaConsensus.consensus_alleles
    Array[File] hla_call_files = hlaConsensus.hla_call_files
    File annotated_vcf = pvacseq.annotated_vcf
    File annotated_tsv = pvacseq.annotated_tsv
    Array[File] pvacseq_predictions = pvacseq.pvacseq_predictions
  }
}
