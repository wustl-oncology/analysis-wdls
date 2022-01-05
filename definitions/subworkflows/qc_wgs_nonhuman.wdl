version 1.0

import "../types.wdl"

import "../subworkflows/hs_metrics.wdl" as hm
import "../tools/bam_to_bigwig.wdl" as btb
import "../tools/collect_alignment_summary_metrics.wdl" as casm
import "../tools/collect_gc_bias_metrics.wdl" as cgbm
import "../tools/collect_insert_size_metrics.wdl" as cism
import "../tools/collect_wgs_metrics.wdl" as cwm
import "../tools/samtools_flagstat.wdl" as sf

workflow qcWgsNonhuman {
  input {
    File bam
    File bam_bai
    File reference
    File reference_fai
    File reference_dict
    String picard_metric_accumulation_level
    Int? minimum_mapping_quality
    Int? minimum_base_quality
    Array[LabelledFile] per_base_intervals
    Array[LabelledFile] per_target_intervals
    Array[LabelledFile] summary_intervals
  }

  call cism.collectInsertSizeMetrics {
    input:
    bam=bam,
    bam_bai=bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    metric_accumulation_level=picard_metric_accumulation_level
  }

  call casm.collectAlignmentSummaryMetrics {
    input:
    bam=bam,
    bam_bai=bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    metric_accumulation_level=picard_metric_accumulation_level
  }

  call cgbm.collectGcBiasMetrics {
    input:
    bam=bam,
    bam_bai=bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    metric_accumulation_level=picard_metric_accumulation_level
  }

  call cwm.collectWgsMetrics {
    input:
    bam=bam,
    bam_bai=bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict
  }

  call sf.samtoolsFlagstat {
    input:
    bam=bam,
    bam_bai=bam_bai
  }

  call hm.hsMetrics as collectHsMetrics {
    input:
    bam=bam,
    bam_bai=bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    minimum_mapping_quality=minimum_mapping_quality,
    minimum_base_quality=minimum_base_quality,
    per_base_intervals=per_base_intervals,
    per_target_intervals=per_target_intervals,
    summary_intervals=summary_intervals
  }

  call btb.bamToBigwig as cgpbigwigBamcoverage {
    input:
    bam=bam,
    bam_bai=bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict
  }

  output {
    File insert_size_metrics = collectInsertSizeMetrics.insert_size_metrics
    File insert_size_histogram = collectInsertSizeMetrics.insert_size_histogram
    File alignment_summary_metrics = collectAlignmentSummaryMetrics.alignment_summary_metrics
    File gc_bias_metrics = collectGcBiasMetrics.gc_bias_metrics
    File wgs_metrics = collectWgsMetrics.wgs_metrics
    File gc_bias_metrics_chart = collectGcBiasMetrics.gc_bias_metrics_chart
    File gc_bias_metrics_summary = collectGcBiasMetrics.gc_bias_metrics_summary
    File flagstats = samtoolsFlagstat.flagstats
    Array[File] per_base_coverage_metrics = collectHsMetrics.per_base_coverage_metrics
    Array[File] per_base_hs_metrics = collectHsMetrics.per_base_hs_metrics
    Array[File] per_target_coverage_metrics = collectHsMetrics.per_target_coverage_metrics
    Array[File] per_target_hs_metrics = collectHsMetrics.per_target_hs_metrics
    Array[File] summary_hs_metrics = collectHsMetrics.summary_hs_metrics
    File bamcoverage_bigwig = cgpbigwigBamcoverage.outfile
  }
}
