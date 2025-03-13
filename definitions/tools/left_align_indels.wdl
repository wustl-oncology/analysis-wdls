version 1.0

task leftAlignIndels {
  input {
    File reference
    File reference_fai
    File reference_dict

    File bam
    File bam_bai

    String output_bam_basename = "left_align_indels"
  }

  Int space_needed_gb = 10 + round(size([bam, bam_bai], "GB")*3 + size([reference, reference_fai, reference_dict], "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "broadinstitute/gatk:4.6.1.0"
    memory: "8GB"
    bootDiskSizeGb: 25
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String outfile_bam = "~{output_bam_basename}.bam"
  String outfile_bam_bai = "~{output_bam_basename}.bam.bai"

  command <<<
    /gatk/gatk --java-options -Xmx6g LeftAlignIndels -O ~{outfile_bam} \
    -R ~{reference} \
    -I ~{bam} \
    --create-output-bam-index
  >>>

  output {
    File left_align_indels_bam = outfile_bam
    File left_align_indels_bam_bai = outfile_bam_bai
  }
}

workflow wf {
  input {
    File reference
    File reference_fai
    File reference_dict
    File bam
    File bam_bai
    String output_bam_basename = "left_align_indels"
  }

  call leftAlignIndels {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=bam,
    bam_bai=bam_bai,
    output_bam_basename=output_bam_basename,
  }
}
