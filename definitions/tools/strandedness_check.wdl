version 1.0

task strandednessCheck {
  input {
    File gtf_file
    File kallisto_index
    File cdna_fasta
    File reads1
    File reads2
  }

  Int space_needed_gb = 10 + round(2*size([gtf_file, kallisto_index, cdna_fasta, reads1, reads2], "GB"))
  runtime {
    memory: "16GB"
    bootDiskSizeGb: space_needed_gb  # default
    cpu: 1
    docker: "mgibio/checkstrandedness:v1"
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  String outfile = basename(reads1, ".fastq") + "strandness_check.txt"
  command <<<
    check_strandedness --print_commands \
        --gtf ~{gtf_file} --kallisto_index ~{kallisto_index} --transcripts ~{cdna_fasta} \
        --reads_1 ~{reads1} --reads_2 ~{reads2} -n 100000 > ~{outfile}
  >>>

  output {
    File strandedness_check = outfile
  }
}

workflow wf {
  call strandednessCheck { input: }
}
