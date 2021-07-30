version 1.0

task indexBam {
  input { File bam }

  Int space_needed_gb = 10 + round(size(bam, "GB")*3)
  runtime {
    docker: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
    memory: "4GB"
    disks: "local-disk ~{space_needed_gb} SSD"
  }
  command <<<
    mv ~{bam} ~{basename(bam)}
    /usr/local/bin/samtools index ~{basename(bam)} ~{basename(bam)}.bai
  >>>
  output {
    File indexed_bam = basename(bam)
    File indexed_bam_bai = "~{basename(bam)}.bai"
  }
}

workflow wf {
  input { File bam }
  call indexBam { input: bam=bam }
}
