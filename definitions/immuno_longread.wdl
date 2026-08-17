version 1.0

import "types.wdl"  # !UnusedImport

import "somatic_germline_wgs_pacbioHifi.wdl" as swgs
import "rnaseq_pacbioHifi.wdl" as rlr
import "pvacseq_longread.wdl" as pslr
import "pvacsplice.wdl" as pspl
import "pvacfuse_longread.wdl" as pflr
import "tools/minimap2_ts_to_xs.wdl" as txs

# Long-read neoantigen pipeline. The counterpart of immuno.wdl.
#
# SUPPORTED INPUT TYPES -- currently PacBio HiFi only:
#   DNA : PacBio HiFi WGS, tumor/normal paired. Either unaligned (or aligned,
#         with skip_align) bams, OR fastq.gz -- see tumor_bams/tumor_fastqs.
#   RNA : PacBio MAS-Iso-Seq FLNC bam(s) from `isoseq refine` -- one or many
#
# Nanopore is NOT supported. tools/minimap2_rnaseq_pacbioHifi.wdl hardcodes
# `-uf`, which is only valid for strand-resolved reads; the DNA arm uses pbmm2
# with --preset HIFI; and DeepSomatic runs its PACBIO model. Adding ONT means
# parameterizing all three, not just swapping a container.
#
# WHAT IS NOT HERE YET
#   * HLA typing. `alleles` is a required input -- supply calls from an
#     orthogonal short-read or clinical source. OptiType/PHLAT/HLA-HD are all
#     short-read tools; a HiFi-native replacement (StarPhase, HLA*LA) is the
#     obvious follow-up.
#   * Nothing else blocking. DNA readcounts for tumor AND normal are handled by
#     add_dna_readcounts (default true) in the somatic workflow; read its
#     normal_sample_column_was_added output on run 1 to learn whether
#     DeepSomatic emitted a normal column on its own.
#
#              tumor/normal HiFi WGS bams OR fastq.gz    FLNC bam(s)
#                              |                          |
#            somatic_germline_wgs_pacbioHifi      rnaseq_pacbioHifi
#                     |            |                 |    |    |
#          VEP vcf ---+            +-- phased -------+    |    +-- aligned bam
#              |                       proximal           |            |
#              |                          |          stringtie         |
#      +-------+--------+                 |          expression        |
#      |                |                 |               |            |
#  pvacseq_longread  pvacsplice <---------+---------------+     pvacfuse_longread
#                    (ts->XS only when
#                     pvacsplice_strand
#                     == "unstranded")
workflow immunoLongread {
  input {
    # =========== DNA: PacBio HiFi WGS ==============================

    # Supply exactly ONE of bams / fastqs per sample. FASTQ mode saves disk:
    # a HiFi uBAM carrying kinetics is several times the size of the equivalent
    # fastq.gz, and nothing downstream here reads the PacBio-only tags that
    # FASTQ drops.
    Array[File] tumor_bams = []
    Array[File] tumor_fastqs = []
    Array[File] normal_bams = []
    Array[File] normal_fastqs = []
    # Optional per-file @RG IDs, positionally matched to the arrays above.
    # e.g. ["NTR004_3591_h2", "NTR004_3591_h3"] so each replicate stays
    # identifiable by @RG ID inside the merged bam.
    Array[String] tumor_read_group_ids = []
    Array[String] normal_read_group_ids = []
    String tumor_sample_name
    String normal_sample_name
    Boolean skip_align = false
    Boolean strip_kinetics = false

    File reference
    File reference_fai
    File reference_dict
    File? reference_mmi

    Int pbmm2_cores = 24
    Int samtools_cores = 8
    Int deepsomatic_cores = 16
    Int clair3_cores = 16
    Int hiphase_cores = 24
    Int chunk_size = 75000000
    String clair3_model = "hifi_revio"
    Boolean run_mosdepth = true
    Boolean rename_vaf_to_af = true
    # Runs bam-readcount on the phased tumor AND normal bams to produce DNA
    # AD/AF/DP for both. On by default: DeepSomatic does not reliably emit a
    # normal genotype column, which would leave --normal-cov / --normal-vaf
    # with nothing to filter on.
    Boolean add_dna_readcounts = true

    # =========== RNA: PacBio MAS-Iso-Seq ===========================

    # One or more unaligned FLNC bams from `isoseq refine`. Multiple entries are
    # replicates of one biological sample -- each is aligned with its own @RG and
    # merged into a single bam, exactly like the DNA arm handles its replicates.
    Array[File] flnc_bams
    # Optional per-replicate @RG IDs, positionally matched to flnc_bams.
    Array[String] flnc_read_group_ids = []
    # Gene model GTF. MUST be the same Ensembl release as the VEP cache below,
    # or StringTie's transcript ids will not join against VEP's CSQ and TX will
    # come out empty with no error.
    File reference_annotation

    # Optional prebuilt minimap2 .mmi for the RNA arm. Saves the genome-indexing
    # step on every run. Must be built with matching -k/-w -- see the caveat in
    # tools/minimap2_rnaseq_pacbioHifi.wdl.
    File? rna_reference_mmi
    String minimap2_preset = "splice:hq"
    String minimap2_additional_args = ""
    Int minimap2_cores = 16
    String? minimap2_docker_image
    String stringtie_strand = "unstranded"
    Boolean stringtie_reference_only = true
    Int stringtie_cores = 12
    String? stringtie_docker_image

    # =========== VEP ===============================================

    # Must be an ENSEMBL cache, not --refseq. pVACtools expects ENST/ENSG ids.
    File vep_cache_dir_zip
    String vep_ensembl_assembly = "GRCh38"
    String vep_ensembl_version
    String vep_ensembl_species = "homo_sapiens"
    File? vep_synonyms_file
    Array[VepCustomAnnotation] vep_custom_annotations = []  # !UnverifiedStruct
    Boolean vep_annotate_coding_only = false
    Boolean vep_everything = true
    Array[String] vep_plugins = ["Frameshift", "Wildtype"]
    String vep_pick = "flag_pick"

    # =========== HLA ===============================================

    # 2-field allele strings, e.g. ["HLA-A*01:01", "HLA-B*07:02", ...]
    Array[String] alleles

    # =========== pVACtools, shared across all three arms ============

    Array[String] prediction_algorithms
    Array[Int]? epitope_lengths_class_i
    Array[Int]? epitope_lengths_class_ii
    Int? binding_threshold
    Int? percentile_threshold
    String? percentile_threshold_strategy
    String? top_score_metric   # enum [lowest, median]
    String? top_score_metric2  # enum [ic50, percentile]
    String? additional_report_columns  # enum [sample_name]
    Int? fasta_size
    Int? downstream_sequence_length
    Boolean? exclude_nas
    String? net_chop_method  # enum [cterm 20s]
    Float? net_chop_threshold
    Boolean? netmhc_stab
    Boolean? run_reference_proteome_similarity
    File? peptide_fasta
    File? genes_of_interest_file
    Int? iedb_retries
    String? netmhciipan_version  # enum [4.3, 4.2, 4.1, 4.0]
    Float? tumor_purity
    Boolean? allele_specific_binding_thresholds
    Int? aggregate_inclusion_binding_threshold
    Int? aggregate_inclusion_count_limit
    Array[String]? problematic_amino_acids
    Boolean? allele_specific_anchors
    Float? anchor_contribution_threshold
    Array[String]? biotypes
    Boolean? allow_incomplete_transcripts
    Array[String]? transcript_prioritization_strategy
    Int? maximum_transcript_support_level

    # Coverage / VAF gates. Left unset by default so pVACtools defaults apply.
    # HiFi WGS runs ~20-30x versus 100x+ for exome, and Iso-Seq read counts are
    # one to two orders of magnitude below short-read counts, so the pVACtools
    # defaults are likely to over-filter on both the DNA and RNA side. Read the
    # mosdepth summary and the RDP column from run 1 before setting these.
    Int? normal_cov
    Int? tdna_cov
    Int? trna_cov
    Float? normal_vaf
    Float? tdna_vaf
    Float? trna_vaf
    Float? expn_val

    Int? readcount_minimum_base_quality
    Int? readcount_minimum_mapping_quality
    Float? minimum_fold_change

    # =========== pVACseq ============================================

    Int? pvacseq_threads

    # =========== pVACsplice / regtools ==============================

    Boolean run_pvacsplice = true

    # regtools strand mode. tools/regtools.wdl maps these onto regtools' -s flag:
    #   "first"      -> RF
    #   "second"     -> FR
    #   "unstranded" -> XS   (reads the XS:A tag)
    #
    # Defaults to "first" (RF) because that is what has worked in practice on
    # long-read data here. Note the knock-on effect: with RF or FR, regtools
    # never looks at XS, so the ts:A -> XS:A conversion below is skipped
    # entirely. Set this to "unstranded" and the conversion turns itself back on
    # (minimap2 emits ts:A, never XS, so regtools would otherwise see no strand
    # information at all).
    String pvacsplice_strand = "first"

    String? regtools_output_filename_tsv
    String? regtools_output_filename_vcf
    String? regtools_output_filename_bed
    Int? regtools_window_size
    Int? max_distance_exon
    Int? max_distance_intron
    Boolean annotate_intronic_variant = false
    Boolean annotate_exonic_variant = false
    Boolean not_skipping_single_exon_transcripts = false
    Boolean singecell_barcode = false
    Boolean intron_motif_priority = false
    Int? pvacsplice_threads
    # Minimum junction read support. tools/pvacsplice.wdl defaults this to 10,
    # which is a short-read number -- long-read junction support is one to two
    # orders of magnitude lower, so 10 discards nearly everything. Overridable
    # from the inputs yml like any other input.
    Int junction_score = 1
    Int? variant_distance
    Boolean? save_gtf
    Array[String]? junction_anchor_types
    Boolean? pvacsplice_keep_tmp_files

    # =========== pVACfuse / CTAT-LR-fusion ==========================

    Boolean run_pvacfuse = true
    # The genome inside this lib must match `reference` above -- ctat-LR-fusion's
    # --LR_bam kickstart mode assumes the bam was aligned to the same genome, so
    # contig names and lengths have to agree.
    File star_fusion_genome_lib_zip
    File agfusion_database
    Boolean agfusion_annotate_noncanonical = true
    # Resume entry point. Set this to an existing AGFusion output directory,
    # zipped, and CTAT-LR-Fusion + AGFusion are both skipped -- pVACfuse runs
    # straight off it. Use this to re-run pVACfuse with different alleles or
    # thresholds without paying for CTAT-LR-Fusion again.
    File? agfusion_dir_zip
    Int ctat_cpu = 30
    Boolean ctat_examine_coding_effect = true
    Boolean ctat_vis = true
    Int? pvacfuse_threads
    Boolean pvacfuse_keep_tmp_files = false
    Int? pvacfuse_read_support
    Float? pvacfuse_expn_val
  }

  # ---------------- DNA arm ----------------

  call swgs.somaticPacbioHifi as somatic {
    input:
    tumor_bams=tumor_bams,
    tumor_fastqs=tumor_fastqs,
    tumor_read_group_ids=tumor_read_group_ids,
    normal_bams=normal_bams,
    normal_fastqs=normal_fastqs,
    normal_read_group_ids=normal_read_group_ids,
    tumor_sample_name=tumor_sample_name,
    normal_sample_name=normal_sample_name,
    skip_align=skip_align,
    strip_kinetics=strip_kinetics,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    reference_mmi=reference_mmi,
    pbmm2_cores=pbmm2_cores,
    samtools_cores=samtools_cores,
    deepsomatic_cores=deepsomatic_cores,
    clair3_cores=clair3_cores,
    hiphase_cores=hiphase_cores,
    chunk_size=chunk_size,
    clair3_model=clair3_model,
    run_mosdepth=run_mosdepth,
    vep_cache_dir_zip=vep_cache_dir_zip,
    vep_ensembl_assembly=vep_ensembl_assembly,
    vep_ensembl_version=vep_ensembl_version,
    vep_ensembl_species=vep_ensembl_species,
    vep_synonyms_file=vep_synonyms_file,
    vep_custom_annotations=vep_custom_annotations,
    vep_annotate_coding_only=vep_annotate_coding_only,
    vep_everything=vep_everything,
    vep_plugins=vep_plugins,
    vep_pick=vep_pick,
    rename_vaf_to_af=rename_vaf_to_af,
    add_dna_readcounts=add_dna_readcounts,
    readcount_minimum_base_quality=readcount_minimum_base_quality,
    readcount_minimum_mapping_quality=readcount_minimum_mapping_quality
  }

  # ---------------- RNA arm ----------------

  call rlr.rnaseqPacbioHifi as rna {
    input:
    flnc_bams=flnc_bams,
    flnc_read_group_ids=flnc_read_group_ids,
    sample_name=tumor_sample_name,
    reference=select_first([rna_reference_mmi, reference]),
    reference_fai=reference_fai,
    reference_annotation=reference_annotation,
    minimap2_preset=minimap2_preset,
    minimap2_additional_args=minimap2_additional_args,
    minimap2_cores=minimap2_cores,
    samtools_cores=samtools_cores,
    minimap2_docker_image=minimap2_docker_image,
    stringtie_strand=stringtie_strand,
    stringtie_reference_only=stringtie_reference_only,
    stringtie_cores=stringtie_cores,
    stringtie_docker_image=stringtie_docker_image
  }

  # ---------------- pVACseq ----------------

  call pslr.pvacseqLongread as pvacseq {
    input:
    detect_variants_vcf=somatic.final_annotated_vcf,
    detect_variants_vcf_tbi=somatic.final_annotated_vcf_tbi,
    sample_name=tumor_sample_name,
    normal_sample_name=normal_sample_name,
    rnaseq_bam=rna.final_bam,
    rnaseq_bam_bai=rna.final_bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    gene_expression_file=rna.stringtie_gene_expression_tsv,
    transcript_expression_file=rna.stringtie_transcript_gtf,
    expression_tool="stringtie",
    phased_proximal_variants_vcf=somatic.phased_proximal_variants_vcf,
    phased_proximal_variants_vcf_tbi=somatic.phased_proximal_variants_vcf_tbi,
    alleles=alleles,
    prediction_algorithms=prediction_algorithms,
    epitope_lengths_class_i=epitope_lengths_class_i,
    epitope_lengths_class_ii=epitope_lengths_class_ii,
    binding_threshold=binding_threshold,
    percentile_threshold=percentile_threshold,
    percentile_threshold_strategy=percentile_threshold_strategy,
    minimum_fold_change=minimum_fold_change,
    top_score_metric=top_score_metric,
    top_score_metric2=top_score_metric2,
    additional_report_columns=additional_report_columns,
    fasta_size=fasta_size,
    downstream_sequence_length=downstream_sequence_length,
    exclude_nas=exclude_nas,
    transcript_prioritization_strategy=transcript_prioritization_strategy,
    maximum_transcript_support_level=maximum_transcript_support_level,
    normal_cov=normal_cov,
    tdna_cov=tdna_cov,
    trna_cov=trna_cov,
    normal_vaf=normal_vaf,
    tdna_vaf=tdna_vaf,
    trna_vaf=trna_vaf,
    expn_val=expn_val,
    readcount_minimum_base_quality=readcount_minimum_base_quality,
    readcount_minimum_mapping_quality=readcount_minimum_mapping_quality,
    net_chop_method=net_chop_method,
    net_chop_threshold=net_chop_threshold,
    netmhc_stab=netmhc_stab,
    run_reference_proteome_similarity=run_reference_proteome_similarity,
    peptide_fasta=peptide_fasta,
    genes_of_interest_file=genes_of_interest_file,
    n_threads=pvacseq_threads,
    iedb_retries=iedb_retries,
    netmhciipan_version=netmhciipan_version,
    tumor_purity=tumor_purity,
    allele_specific_binding_thresholds=allele_specific_binding_thresholds,
    aggregate_inclusion_binding_threshold=aggregate_inclusion_binding_threshold,
    aggregate_inclusion_count_limit=aggregate_inclusion_count_limit,
    problematic_amino_acids=problematic_amino_acids,
    allele_specific_anchors=allele_specific_anchors,
    anchor_contribution_threshold=anchor_contribution_threshold,
    biotypes=biotypes,
    allow_incomplete_transcripts=allow_incomplete_transcripts,
    prefix="variants.final"
  }

  # ---------------- pVACsplice ----------------
  #
  # Only needed when pvacsplice_strand is "unstranded", which sends regtools
  # `-s XS`. minimap2 emits ts:A and never XS, so without this conversion
  # regtools would see no strand information and find no junctions. With the
  # default "first" (RF) regtools ignores XS and this whole step is skipped.

  if (run_pvacsplice && pvacsplice_strand == "unstranded") {
    call txs.minimap2TsToXs as xsTag {
      input:
      bam=rna.final_bam,
      bam_bai=rna.final_bam_bai,
      output_basename=tumor_sample_name + ".rna.xs_tagged",
      cores=samtools_cores
    }
  }

  if (run_pvacsplice) {
    File pvacsplice_bam = select_first([xsTag.xs_tagged_bam, rna.final_bam])
    File pvacsplice_bam_bai = select_first([xsTag.xs_tagged_bam_bai, rna.final_bam_bai])

    call pspl.pvacsplice as pvacsplice {
      input:
      detect_variants_vcf=somatic.final_annotated_vcf,
      detect_variants_vcf_tbi=somatic.final_annotated_vcf_tbi,
      sample_name=tumor_sample_name,
      normal_sample_name=normal_sample_name,
      rnaseq_bam=pvacsplice_bam,
      rnaseq_bam_bai=pvacsplice_bam_bai,
      reference=reference,
      reference_fai=reference_fai,
      reference_dict=reference_dict,
      reference_annotation=reference_annotation,
      gene_expression_file=rna.stringtie_gene_expression_tsv,
      transcript_expression_file=rna.stringtie_transcript_gtf,
      expression_tool="stringtie",
      readcount_minimum_base_quality=readcount_minimum_base_quality,
      readcount_minimum_mapping_quality=readcount_minimum_mapping_quality,
      peptide_fasta=peptide_fasta,
      genes_of_interest_file=genes_of_interest_file,
      output_filename_tsv=regtools_output_filename_tsv,
      output_filename_vcf=regtools_output_filename_vcf,
      output_filename_bed=regtools_output_filename_bed,
      strand=pvacsplice_strand,
      window_size=regtools_window_size,
      max_distance_exon=max_distance_exon,
      max_distance_intron=max_distance_intron,
      annotate_intronic_variant=annotate_intronic_variant,
      annotate_exonic_variant=annotate_exonic_variant,
      not_skipping_single_exon_transcripts=not_skipping_single_exon_transcripts,
      singecell_barcode=singecell_barcode,
      intron_motif_priority=intron_motif_priority,
      alleles=alleles,
      prediction_algorithms=prediction_algorithms,
      epitope_lengths_class_i=epitope_lengths_class_i,
      epitope_lengths_class_ii=epitope_lengths_class_ii,
      binding_threshold=binding_threshold,
      percentile_threshold=percentile_threshold,
      percentile_threshold_strategy=percentile_threshold_strategy,
      iedb_retries=iedb_retries,
      top_score_metric=top_score_metric,
      top_score_metric2=top_score_metric2,
      additional_report_columns=additional_report_columns,
      fasta_size=fasta_size,
      exclude_nas=exclude_nas,
      transcript_prioritization_strategy=transcript_prioritization_strategy,
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
      n_threads=pvacsplice_threads,
      netmhciipan_version=netmhciipan_version,
      tumor_purity=tumor_purity,
      allele_specific_binding_thresholds=allele_specific_binding_thresholds,
      aggregate_inclusion_binding_threshold=aggregate_inclusion_binding_threshold,
      aggregate_inclusion_count_limit=aggregate_inclusion_count_limit,
      problematic_amino_acids=problematic_amino_acids,
      biotypes=biotypes,
      allow_incomplete_transcripts=allow_incomplete_transcripts,
      junction_score=junction_score,
      variant_distance=variant_distance,
      save_gtf=save_gtf,
      junction_anchor_types=junction_anchor_types,
      keep_tmp_files=pvacsplice_keep_tmp_files
    }
  }

  # ---------------- pVACfuse ----------------

  if (run_pvacfuse) {
    call pflr.pvacfuse_longread as pvacfuse {
      input:
      lr_bam=rna.final_bam,
      lr_bam_bai=rna.final_bam_bai,
      agfusion_dir_zip=agfusion_dir_zip,
      star_fusion_genome_lib_zip=star_fusion_genome_lib_zip,
      ctat_cpu=ctat_cpu,
      examine_coding_effect=ctat_examine_coding_effect,
      vis=ctat_vis,
      agfusion_database=agfusion_database,
      agfusion_annotate_noncanonical=agfusion_annotate_noncanonical,
      sample_name=tumor_sample_name,
      alleles=alleles,
      prediction_algorithms=prediction_algorithms,
      epitope_lengths_class_i=epitope_lengths_class_i,
      epitope_lengths_class_ii=epitope_lengths_class_ii,
      binding_threshold=binding_threshold,
      percentile_threshold=percentile_threshold,
      percentile_threshold_strategy=percentile_threshold_strategy,
      iedb_retries=iedb_retries,
      pvacfuse_keep_tmp_files=pvacfuse_keep_tmp_files,
      net_chop_method=net_chop_method,
      net_chop_threshold=net_chop_threshold,
      top_score_metric=top_score_metric,
      top_score_metric2=top_score_metric2,
      run_reference_proteome_similarity=run_reference_proteome_similarity,
      additional_report_columns=additional_report_columns,
      fasta_size=fasta_size,
      downstream_sequence_length=downstream_sequence_length,
      pvacfuse_threads=select_first([pvacfuse_threads, 8]),
      allele_specific_binding_thresholds=select_first([allele_specific_binding_thresholds, false]),
      aggregate_inclusion_binding_threshold=aggregate_inclusion_binding_threshold,
      aggregate_inclusion_count_limit=aggregate_inclusion_count_limit,
      problematic_amino_acids=problematic_amino_acids,
      peptide_fasta=peptide_fasta,
      genes_of_interest_file=genes_of_interest_file,
      read_support=pvacfuse_read_support,
      expn_val=pvacfuse_expn_val,
      netmhciipan_version=netmhciipan_version
    }
  }

  output {
    # ---------- DNA ----------
    String tumor_input_mode = somatic.tumor_input_mode
    String normal_input_mode = somatic.normal_input_mode
    File tumor_bam = somatic.tumor_bam
    File tumor_bam_bai = somatic.tumor_bam_bai
    File normal_bam = somatic.normal_bam
    File normal_bam_bai = somatic.normal_bam_bai
    File tumor_bam_phased = somatic.tumor_bam_phased
    File normal_bam_phased = somatic.normal_bam_phased

    File somatic_vcf = somatic.somatic_vcf
    File somatic_vcf_phased = somatic.somatic_vcf_phased
    File tumor_germline_vcf = somatic.tumor_germline_vcf
    File normal_germline_vcf = somatic.normal_germline_vcf
    File normal_germline_vcf_phased = somatic.normal_germline_vcf_phased

    File final_annotated_vcf = somatic.final_annotated_vcf
    File final_annotated_vcf_tbi = somatic.final_annotated_vcf_tbi
    File vep_summary = somatic.vep_summary
    File phased_proximal_variants_vcf = somatic.phased_proximal_variants_vcf

    # ---------- RNA ----------
    File rna_bam = rna.final_bam
    File rna_bam_bai = rna.final_bam_bai
    File rna_flagstats = rna.flagstats
    File stringtie_transcript_gtf = rna.stringtie_transcript_gtf
    File stringtie_gene_expression_tsv = rna.stringtie_gene_expression_tsv
    File? rna_bam_xs_tagged = xsTag.xs_tagged_bam

    # ---------- QC ----------
    File? mosdepth_tumor_summary = somatic.mosdepth_tumor_summary
    File? mosdepth_normal_summary = somatic.mosdepth_normal_summary
    File hiphase_tumor_stats = somatic.hiphase_tumor_stats
    File hiphase_tumor_summary = somatic.hiphase_tumor_summary
    File hiphase_tumor_blocks = somatic.hiphase_tumor_blocks

    # ---------- pVACseq ----------
    Array[File] pvacseq_mhc_i = pvacseq.mhc_i
    Array[File] pvacseq_mhc_ii = pvacseq.mhc_ii
    Array[File] pvacseq_combined = pvacseq.combined
    File? pvacseq_mhc_i_log = pvacseq.mhc_i_log
    File? pvacseq_mhc_ii_log = pvacseq.mhc_ii_log
    File pvacseq_annotated_expression_vcf = pvacseq.annotated_vcf
    File pvacseq_annotated_expression_vcf_tbi = pvacseq.annotated_vcf_tbi
    File variants_final_annotated_tsv = pvacseq.annotated_tsv

    # ---------- pVACsplice ----------
    Array[File]? pvacsplice_mhc_i = pvacsplice.mhc_i
    Array[File]? pvacsplice_mhc_ii = pvacsplice.mhc_ii
    Array[File]? pvacsplice_combined = pvacsplice.combined
    File? pvacsplice_splice_transcript_combined_report = pvacsplice.splice_transcript_combined_report
    File? pvacsplice_splice_fasta = pvacsplice.splice_fasta
    File? pvacsplice_splice_gtf = pvacsplice.splice_gtf
    File? pvacsplice_annotated_vcf = pvacsplice.annotated_vcf

    # ---------- pVACfuse ----------
    Array[File]? pvacfuse_mhc_i = pvacfuse.mhc_i_all
    Array[File]? pvacfuse_mhc_ii = pvacfuse.mhc_ii_all
    Array[File]? pvacfuse_combined = pvacfuse.combined_all
    # Absent when agfusion_dir_zip was supplied -- CTAT-LR-Fusion never ran.
    File? pvacfuse_ctat_fusion_predictions = pvacfuse.ctat_fusion_predictions
    File? pvacfuse_ctat_fusion_predictions_abridged = pvacfuse.ctat_fusion_predictions_abridged
    File? pvacfuse_agfusion_zip = pvacfuse.agfusion_annotated_predictions_zip
  }
}
