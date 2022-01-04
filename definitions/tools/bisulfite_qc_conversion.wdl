version 1.0

task bisulfiteQcConversion {
  input {
    File vcf
    File bam
    File reference
    File QCannotation
  }

  Int space_needed_gb = 10 + round(size([vcf, bam, reference, QCannotation], "GB"))
  runtime {
    cpu: 1
    memory: "16GB"
    bootDiskSizeGb: 20
    docker: "mgibio/biscuit:0.3.8.2"
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  command <<<
    /bin/bash /opt/biscuit/scripts/Bisulfite_QC_bisulfiteconversion.sh \
    ~{vcf} ~{bam} ~{reference} ~{QCannotation}
  >>>

  output {
    File base_conversion = "totalBaseConversionRate.txt"
    File read_conversion = "totalReadConversionRate.txt"
    File cph_retention = "CpHRetentionByReadPos.txt"
    File cpg_retention = "CpGRetentionByReadPos.txt"
  }
}
