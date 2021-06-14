version 1.0

task indexBam {
  input { File bam }

  Int space_needed_gb = 10 + round(size(bam, "GB")*2)
  runtime {
    docker: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
    memory: "4GB"
    disks: "local-disk ~{space_needed_gb} HDD"
  }
  command <<<
    /usr/local/bin/samtools index ~{bam} "~{basename(bam)}.bai"
  >>>
  output {
    File indexed_bam = bam
    File indexed_bam_bai = "~{basename(bam)}.bai"
  }
}

workflow wf {
  input { File bam }
  call indexBam { input: bam=bam }
}
