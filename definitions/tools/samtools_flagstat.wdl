version 1.0

task samtoolsFlagstat {
  input {
    File bam
    File bam_bai
  }

  Int space_needed_gb = 10 + round(size([bam, bam_bai], "GB")*2)
  runtime {
    docker: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
    memory: "4GB"
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  String outfile = basename(bam) + ".flagstat"
  command <<<
    /usr/local/bin/samtools flagstat ~{bam} > ~{outfile}
  >>>

  output {
    File flagstats = outfile
  }
}

workflow wf {
  input {
    File bam
    File bam_bai
  }

  call samtoolsFlagstat {
    input:
    bam=bam,
    bam_bai=bam_bai
  }
}
