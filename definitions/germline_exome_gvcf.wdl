version 1.0

import "types.wdl"
import "alignment_exome.wdl" as ae
import "tools/index_cram.wdl" as ic
import "tools/bam_to_cram.wdl" as btc
import "subworkflows/gatk_haplotypecaller_iterator.wdl" as ghi

workflow germlineExomeGvcf {
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
    Array[LabelledFile] per_base_intervals
    Array[LabelledFile] per_target_intervals
    Array[LabelledFile] summary_intervals
    File omni_vcf
    File omni_vcf_tbi
    String picard_metric_accumulation_level
    String emit_reference_confidence  # enum ['NONE', 'BP_RESOLUTION', 'GVCF']
    Array[String] gvcf_gq_bands
    Array[Array[String]] intervals
    Int? ploidy
    Int? qc_minimum_mapping_quality
    Int? qc_minimum_base_quality
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
    verify_bam_id_metrics=alignmentAndQc.verify_bam_id_metrics,
    ploidy=ploidy
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
    Array[File] per_targetCoverage_metrics = alignmentAndQc.per_target_coverage_metrics
    Array[File] per_target_hs_metrics = alignmentAndQc.per_target_hs_metrics
    Array[File] per_baseCoverage_metrics = alignmentAndQc.per_base_coverage_metrics
    Array[File] per_base_hs_metrics = alignmentAndQc.per_base_hs_metrics
    Array[File] summary_hs_metrics = alignmentAndQc.summary_hs_metrics
    File flagstats = alignmentAndQc.flagstats
    File verify_bam_id_metrics = alignmentAndQc.verify_bam_id_metrics
    File verify_bam_id_depth = alignmentAndQc.verify_bam_id_depth
    Array[File] gvcf = generateGvcfs.gvcf
  }
}
