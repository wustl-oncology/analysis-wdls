version 1.0

import "types.wdl"
import "subworkflows/sequence_to_bqsr.wdl" as s2b
import "subworkflows/qc_exome.wdl" as qe

workflow alignmentExome {
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
    String? final_name
    File bait_intervals
    File target_intervals
    Array[LabelledFile] per_base_intervals
    Array[LabelledFile] per_target_intervals
    Array[LabelledFile] summary_intervals
    File omni_vcf
    File omni_vcf_tbi
    String picard_metric_accumulation_level
    Int? qc_minimum_mapping_quality
    Int? qc_minimum_base_quality
  }

  call s2b.sequenceToBqsr as alignment {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    reference_amb=reference_amb,
    reference_ann=reference_ann,
    reference_bwt=reference_bwt,
    reference_pac=reference_pac,
    reference_sa=reference_sa,
    unaligned=sequence,
    trimming=trimming,
    bqsr_known_sites=bqsr_known_sites,
    bqsr_known_sites_tbi=bqsr_known_sites_tbi,
    bqsr_intervals=bqsr_intervals,
    final_name=final_name
  }

  call qe.qcExome as qc {
    input:
    bam=alignment.final_bam,
    bam_bai=alignment.final_bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bait_intervals=bait_intervals,
    target_intervals=target_intervals,
    per_base_intervals=per_base_intervals,
    per_target_intervals=per_target_intervals,
    summary_intervals=summary_intervals,
    omni_vcf=omni_vcf,
    omni_vcf_tbi=omni_vcf_tbi,
    picard_metric_accumulation_level=picard_metric_accumulation_level,
    minimum_mapping_quality=qc_minimum_mapping_quality,
    minimum_base_quality=qc_minimum_base_quality
  }

  output {
    File bam = alignment.final_bam
    File bam_bai = alignment.final_bam_bai
    File mark_duplicates_metrics = alignment.mark_duplicates_metrics_file
    File insert_size_metrics = qc.insert_size_metrics
    File insert_size_histogram = qc.insert_size_histogram
    File alignment_summary_metrics = qc.alignment_summary_metrics
    File hs_metrics = qc.hs_metrics
    Array[File] per_target_coverage_metrics = qc.per_target_coverage_metrics
    Array[File] per_target_hs_metrics = qc.per_target_hs_metrics
    Array[File] per_base_coverage_metrics = qc.per_base_coverage_metrics
    Array[File] per_base_hs_metrics = qc.per_base_hs_metrics
    Array[File] summary_hs_metrics = qc.summary_hs_metrics
    File flagstats = qc.flagstats
    File verify_bam_id_metrics = qc.verify_bam_id_metrics
    File verify_bam_id_depth = qc.verify_bam_id_depth
  }
}
