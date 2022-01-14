version 1.0

task bisulfiteQcCpgRetentionDistribution {
  input {
    File vcf
    File bam
    File reference
    File reference_fai
    File QCannotation
  }

  Int space_needed_gb = 10 + round(size([vcf, bam, reference, QCannotation], "GB"))
  runtime {
    cpu: 1
    memory: "16GB"
    bootDiskSizeGb: 10
    docker: "mgibio/biscuit:0.3.8.2"
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  command <<<
    /bin/bash /opt/biscuit/scripts/Bisulfite_QC_CpGretentiondistribution.sh \
    ~{vcf} ~{bam} ~{reference} ~{QCannotation}
  >>>

  output {
    File cpg_retention_dist = "CpGRetentionDist.txt"
  }
}
