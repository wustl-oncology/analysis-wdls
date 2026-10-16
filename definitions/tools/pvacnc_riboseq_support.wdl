version 1.0

task pvacncRiboseqSupport {
  input {
    String sample_name
    String caller
    File mapping_tsv
    File noncanonical_orf_blocks
    File? reference_annotation
    File? filtered_report
    Array[File] riboseq_bams
    Array[File] riboseq_bais
    Array[String] riboseq_sample_labels
    Int minimum_mapping_quality = 0
    Int minimum_base_quality = 0
    Int minimum_alt_reads = 1
    Boolean exclude_secondary = false
    Boolean exclude_supplementary = false
    Boolean exclude_qcfail = false
    Boolean exclude_duplicates = false
    String docker_image = "jinglunli/pvacnc:latest"
  }

  String output_prefix = sample_name + "." + caller
  Int space_needed_gb = 10 + round(size(riboseq_bams, "GB"))

  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "16GB"
    cpu: 2
    docker: docker_image
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    set -e
    pvacnc-riboseq-support \
      --sample-name ~{sample_name} \
      --caller ~{caller} \
      --mapping ~{mapping_tsv} \
      ~{if defined(filtered_report) then "--filtered-report " + select_first([filtered_report]) else ""} \
      --bam ~{sep=' --bam ' riboseq_bams} \
      --bam-index ~{sep=' --bam-index ' riboseq_bais} \
      --label ~{sep=' --label ' riboseq_sample_labels} \
      --minimum-mapping-quality ~{minimum_mapping_quality} \
      --minimum-base-quality ~{minimum_base_quality} \
      --minimum-alt-reads ~{minimum_alt_reads} \
      ~{if exclude_secondary then "--exclude-secondary" else ""} \
      ~{if exclude_supplementary then "--exclude-supplementary" else ""} \
      ~{if exclude_qcfail then "--exclude-qcfail" else ""} \
      ~{if exclude_duplicates then "--exclude-duplicates" else ""} \
      --output-prefix ~{output_prefix}

    pvacnc-riboseq-orf-support \
      --candidate-summary ~{output_prefix}.riboseq_candidate_summary.tsv \
      --noncanonical-blocks ~{noncanonical_orf_blocks} \
      ~{if defined(reference_annotation) then "--reference-gtf " + select_first([reference_annotation]) else ""} \
      --bam ~{sep=' --bam ' riboseq_bams} \
      --bam-index ~{sep=' --bam-index ' riboseq_bais} \
      --label ~{sep=' --label ' riboseq_sample_labels} \
      --minimum-mapping-quality ~{minimum_mapping_quality} \
      ~{if exclude_secondary then "--exclude-secondary" else ""} \
      ~{if exclude_supplementary then "--exclude-supplementary" else ""} \
      ~{if exclude_qcfail then "--exclude-qcfail" else ""} \
      ~{if exclude_duplicates then "--exclude-duplicates" else ""} \
      --output-prefix ~{output_prefix}
  >>>

  output {
    File orf_support_by_bam = "~{output_prefix}.riboseq_orf_support_by_bam.tsv"
    File orf_candidate_summary = "~{output_prefix}.riboseq_orf_candidate_summary.tsv"
    File support_by_bam = "~{output_prefix}.riboseq_support_by_bam.tsv"
    File allele_summary = "~{output_prefix}.riboseq_allele_summary.tsv"
    File candidate_summary = "~{output_prefix}.riboseq_candidate_summary.tsv"
    File run_summary = "~{output_prefix}.riboseq_run_summary.tsv"
  }
}

workflow wf {
  input {
    String sample_name
    String caller
    File mapping_tsv
    File noncanonical_orf_blocks
    File? reference_annotation
    File? filtered_report
    Array[File] riboseq_bams
    Array[File] riboseq_bais
    Array[String] riboseq_sample_labels
    Int minimum_mapping_quality = 0
    Int minimum_base_quality = 0
    Int minimum_alt_reads = 1
    Boolean exclude_secondary = false
    Boolean exclude_supplementary = false
    Boolean exclude_qcfail = false
    Boolean exclude_duplicates = false
    String docker_image = "jinglunli/pvacnc:latest"
  }

  call pvacncRiboseqSupport {
    input:
      sample_name = sample_name,
      caller = caller,
      mapping_tsv = mapping_tsv,
      noncanonical_orf_blocks = noncanonical_orf_blocks,
      reference_annotation = reference_annotation,
      filtered_report = filtered_report,
      riboseq_bams = riboseq_bams,
      riboseq_bais = riboseq_bais,
      riboseq_sample_labels = riboseq_sample_labels,
      minimum_mapping_quality = minimum_mapping_quality,
      minimum_base_quality = minimum_base_quality,
      minimum_alt_reads = minimum_alt_reads,
      exclude_secondary = exclude_secondary,
      exclude_supplementary = exclude_supplementary,
      exclude_qcfail = exclude_qcfail,
      exclude_duplicates = exclude_duplicates,
      docker_image = docker_image
  }

  output {
    File orf_support_by_bam = pvacncRiboseqSupport.orf_support_by_bam
    File orf_candidate_summary = pvacncRiboseqSupport.orf_candidate_summary
    File support_by_bam = pvacncRiboseqSupport.support_by_bam
    File allele_summary = pvacncRiboseqSupport.allele_summary
    File candidate_summary = pvacncRiboseqSupport.candidate_summary
    File run_summary = pvacncRiboseqSupport.run_summary
  }
}
