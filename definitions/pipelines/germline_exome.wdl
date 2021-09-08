version 1.0

import "../types.wdl"

import "../pipelines/alignment_exome.wdl" as ae
import "../tools/interval_list_expand.wdl" as ile
import "../subworkflows/germline_detect_variants.wdl" as gdv
import "../tools/bam_to_cram.wdl" as btc
import "../tools/index_cram.wdl" as ic

workflow germlineExome {
  input {
    File reference
    File reference_fai
    File reference_dict
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_sa
    Array[SequenceData] sequence
    TrimmingOptions? trimming
    Array[File] bqsr_known_sites
    Array[File] bqsr_known_sites_tbi
    Array[String]? bqsr_intervals
    File bait_intervals
    File target_intervals
    Int target_interval_padding = 100
    Array[LabelledFile] per_base_intervals
    Array[LabelledFile] per_target_intervals
    Array[LabelledFile] summary_intervals
    File omni_vcf
    File omni_vcf_tbi
    String picard_metric_accumulation_level
    Array[String] gvcf_gq_bands
    Array[Array[String]] intervals
    Int? ploidy
    File vep_cache_dir_zip
    String vep_ensembl_assembly
    String vep_ensembl_version
    String vep_ensembl_species
    Array[String]? vep_plugins
    File? synonyms_file
    Boolean? annotate_coding_only
    Int? qc_minimum_mapping_quality
    Int? qc_minimum_base_quality
    Array[VepCustomAnnotation] vep_custom_annotations
    Array[String]? variants_to_table_fields
    Array[String]? variants_to_table_genotype_fields
    Array[String]? vep_to_table_fields
  }

  call ae.alignmentExome as alignmentAndQc {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    reference_amb=reference_amb,
    reference_ann=reference_ann,
    reference_bwt=reference_bwt,
    reference_pac=reference_pac,
    reference_sa=reference_sa,
    sequence=sequence,
    trimming=trimming,
    bqsr_known_sites=bqsr_known_sites,
    bqsr_known_sites_tbi=bqsr_known_sites_tbi,
    bqsr_intervals=bqsr_intervals,
    bait_intervals=bait_intervals,
    target_intervals=target_intervals,
    per_base_intervals=per_base_intervals,
    per_target_intervals=per_target_intervals,
    summary_intervals=summary_intervals,
    omni_vcf=omni_vcf,
    omni_vcf_tbi=omni_vcf_tbi,
    picard_metric_accumulation_level=picard_metric_accumulation_level,
    qc_minimum_mapping_quality=qc_minimum_mapping_quality,
    qc_minimum_base_quality=qc_minimum_base_quality
  }

  call ile.intervalListExpand as padTargetIntervals {
    input:
    interval_list=target_intervals,
    roi_padding=target_interval_padding
  }

  call gdv.germlineDetectVariants as detectVariants {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=alignmentAndQc.bam,
    bai=alignmentAndQc.bam_bai,
    gvcf_gq_bands=gvcf_gq_bands,
    intervals=intervals,
    verify_bam_id_metrics=alignmentAndQc.verify_bam_id_metrics,
    ploidy=ploidy,
    vep_cache_dir_zip=vep_cache_dir_zip,
    synonyms_file=synonyms_file,
    annotate_coding_only=annotate_coding_only,
    limit_variant_intervals=padTargetIntervals.expanded_interval_list,
    vep_ensembl_assembly=vep_ensembl_assembly,
    vep_ensembl_version=vep_ensembl_version,
    vep_ensembl_species=vep_ensembl_species,
    vep_plugins=vep_plugins,
    vep_to_table_fields=vep_to_table_fields,
    vep_custom_annotations=vep_custom_annotations,
    variants_to_table_fields=variants_to_table_fields,
    variants_to_table_genotype_fields=variants_to_table_genotype_fields
  }

  call btc.bamToCram {
    input:
    bam=alignmentAndQc.bam,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict
  }

  call ic.indexCram {
    input: cram=bamToCram.cram
  }

  output {
    File cram = indexCram.indexed_cram
    File cram_crai = indexCram.indexed_cram_crai
    File mark_duplicates_metrics = alignmentAndQc.mark_duplicates_metrics
    File insert_size_metrics = alignmentAndQc.insert_size_metrics
    File insert_size_histogram = alignmentAndQc.insert_size_histogram
    File alignment_summary_metrics = alignmentAndQc.alignment_summary_metrics
    File hs_metrics = alignmentAndQc.hs_metrics
    Array[File] per_target_coverage_metrics = alignmentAndQc.per_target_coverage_metrics
    Array[File] per_target_hs_metrics = alignmentAndQc.per_target_hs_metrics
    Array[File] per_base_coverage_metrics = alignmentAndQc.per_base_coverage_metrics
    Array[File] per_base_hs_metrics = alignmentAndQc.per_base_hs_metrics
    Array[File] summary_hs_metrics = alignmentAndQc.summary_hs_metrics
    File flagstats = alignmentAndQc.flagstats
    File verify_bam_id_metrics = alignmentAndQc.verify_bam_id_metrics
    File verify_bam_id_depth = alignmentAndQc.verify_bam_id_depth
    File raw_vcf = detectVariants.raw_vcf
    File raw_vcf_tbi = detectVariants.raw_vcf_tbi
    File final_vcf = detectVariants.final_vcf
    File final_vcf_tbi = detectVariants.final_vcf_tbi
    File filtered_vcf = detectVariants.filtered_vcf
    File filtered_vcf_tbi = detectVariants.filtered_vcf_tbi
    File vep_summary = detectVariants.vep_summary
    File final_tsv = detectVariants.final_tsv
    File filtered_tsv = detectVariants.filtered_tsv
  }
}
