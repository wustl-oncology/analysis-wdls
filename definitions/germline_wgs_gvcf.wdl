version 1.0

import "types.wdl"

import "alignment_wgs.wdl" as aw

import "subworkflows/gatk_haplotypecaller_iterator.wdl" as ghi
import "tools/index_cram.wdl" as ic
import "tools/bam_to_cram.wdl" as btc
import "tools/freemix.wdl" as f

workflow germlineWgsGvcf {
  input {
    File reference
    File reference_fai
    File reference_dict
    File reference_alt
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_0123
    Array[SequenceData] sequence
    TrimmingOptions? trimming
    File omni_vcf
    File omni_vcf_tbi
    String picard_metric_accumulation_level
    String emit_reference_confidence  # enum ["NONE", "BP_RESOLUTION", "GVCF"]
    Array[String] gvcf_gq_bands
    Array[Array[String]] intervals
    Int? ploidy
    File qc_intervals
    File? synonyms_file
    Boolean? annotate_coding_only
    Array[File] bqsr_known_sites
    Array[File] bqsr_known_sites_tbi
    Array[String]? bqsr_intervals
    Int? minimum_mapping_quality
    Int? minimum_base_quality
    Array[LabelledFile] per_base_intervals
    Array[LabelledFile] per_target_intervals
    Array[LabelledFile] summary_intervals
  }

  call aw.alignmentWgs as alignmentAndQc {
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
    sequence=sequence,
    trimming=trimming,
    bqsr_known_sites=bqsr_known_sites,
    bqsr_known_sites_tbi=bqsr_known_sites_tbi,
    omni_vcf=omni_vcf,
    omni_vcf_tbi=omni_vcf_tbi,
    intervals=qc_intervals,
    picard_metric_accumulation_level=picard_metric_accumulation_level,
    bqsr_intervals=bqsr_intervals,
    minimum_mapping_quality=minimum_mapping_quality,
    minimum_base_quality=minimum_base_quality,
    per_base_intervals=per_base_intervals,
    per_target_intervals=per_target_intervals,
    summary_intervals=summary_intervals
  }

  call f.freemix {
    input:
    verify_bam_id_metrics=alignmentAndQc.verify_bam_id_metrics
  }

  call ghi.gatkHaplotypecallerIterator as generateGvcfs {
    input:
    bam=alignmentAndQc.bam,
    bai=alignmentAndQc.bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    emit_reference_confidence=emit_reference_confidence,
    gvcf_gq_bands=gvcf_gq_bands,
    intervals=intervals,
    ploidy=ploidy,
    contamination_fraction=freemix.out
  }

  call btc.bamToCram {
    input:
    bam=alignmentAndQc.bam,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict
  }

  call ic.indexCram {
    input:
    cram=bamToCram.cram
  }

  output {
    File cram = indexCram.indexed_cram
    File mark_duplicates_metrics = alignmentAndQc.mark_duplicates_metrics
    File insert_size_metrics = alignmentAndQc.insert_size_metrics
    File insert_size_histogram = alignmentAndQc.insert_size_histogram
    File alignment_summary_metrics = alignmentAndQc.alignment_summary_metrics
    File gc_bias_metrics = alignmentAndQc.gc_bias_metrics
    File gc_bias_metrics_chart = alignmentAndQc.gc_bias_metrics_chart
    File gc_bias_metrics_summary = alignmentAndQc.gc_bias_metrics_summary
    File wgs_metrics = alignmentAndQc.wgs_metrics
    File flagstats = alignmentAndQc.flagstats
    File verify_bam_id_metrics = alignmentAndQc.verify_bam_id_metrics
    File verify_bam_id_depth = alignmentAndQc.verify_bam_id_depth
    Array[File] per_base_coverage_metrics = alignmentAndQc.per_base_coverage_metrics
    Array[File] per_base_hs_metrics = alignmentAndQc.per_base_hs_metrics
    Array[File] per_target_coverage_metrics = alignmentAndQc.per_target_coverage_metrics
    Array[File] per_target_hs_metrics = alignmentAndQc.per_target_hs_metrics
    Array[File] summary_hs_metrics = alignmentAndQc.summary_hs_metrics
    Array[File] gvcf = generateGvcfs.gvcf
  }
}
