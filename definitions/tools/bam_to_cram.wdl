version 1.0

task bamToCram {
  input {
    File reference
    File reference_fai
    File reference_dict
    File bam
  }

  Int size_needed_gb = 10 + round(size([reference, reference_fai, reference_dict, bam], "GB"))
  runtime {
    docker: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
    memory: "4GB"
    disks: "local-disk ~{size_needed_gb} HDD"
  }

  String outfile = basename(bam, ".bam") + ".cram"
  command <<<
    /usr/local/bin/samtools view -C -T ~{reference} ~{bam} > ~{outfile}
  >>>

  output {
    File cram = outfile
  }
}
