version 1.0

task applyBqsr {
  input {
    File reference
    File reference_fai
    File reference_dict
    File bam
    File bam_bai
    File bqsr_table
    String output_name = "final"
  }

  Int space_needed_gb = 10 + round(size([bqsr_table, reference, reference_fai, reference_dict], "GB") + size([bam, bam_bai], "GB") * 2)
  runtime {
    docker: "broadinstitute/gatk:4.1.8.1"
    memory: "18GB"
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  command <<<
    /gatk/gatk --java-options -Xmx16g ApplyBQSR -O ~{output_name}.bam ~{sep=" " prefix("--static-quantized-quals ", [10, 20, 30])} -R ~{reference} -I ~{bam} -bqsr ~{bqsr_table}
  >>>

  output {
    File bqsr_bam = "~{output_name}.bam"
    File bqsr_bam_bai = "~{output_name}.bai"
  }
}

workflow wf {
  input {
    File reference
    File reference_fai
    File reference_dict
    File bam
    File bam_bai
    File bqsr_table
    String output_name = "final"
  }
  call applyBqsr {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=bam,
    bam_bai=bam_bai,
    bqsr_table=bqsr_table,
    output_name=output_name
  }
}
