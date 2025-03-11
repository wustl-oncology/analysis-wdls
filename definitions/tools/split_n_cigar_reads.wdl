version 1.0

task splitNCigarReads {
  input {
    File reference
    File reference_fai
    File reference_dict

    File bam
    File bam_bai

    String output_bam_basename = "split_n_cigar"

  Int space_needed_gb = 10 + round(size([bam, bam_bai], "GB")*3 + size([reference, reference_fai, reference_dict], "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "broadinstitute/gatk:4.6.1.0"
    memory: "32GB"
    bootDiskSizeGb: 25
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String outfile = "~{output_bam_basename}.bam"

  command <<<
    /gatk/gatk --java-options -Xmx28g SplitNCigarReads -O ~{outfile} \
    -R ~{reference} \
    -I ~{bam}
  >>>

  output {
    File split_n_cigar_bam = outfile
  }
}

workflow wf {
  input {
    File reference
    File reference_fai
    File reference_dict
    File bam
    File bam_bai
    String output_bam_basename = "split_n_cigar"
  }

  call splitNCigarReads {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=bam,
    bam_bai=bam_bai,
    output_bam_basename=output_bam_basename,
  }
}
