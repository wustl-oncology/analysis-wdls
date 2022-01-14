version 1.0

task bisulfiteQcCoverageStats {
  input {
    File vcf
    File bam
    File reference
    File reference_fai
    File QCannotation
  }

  Int space_needed_gb = 10 + round(size([vcf, bam, reference, reference_fai, QCannotation], "GB"))
  runtime {
    cpu: 1
    memory: "16GB"
    bootDiskSizeGb: 20
    docker: "mgibio/biscuit:0.3.8.2"
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  command <<<
    /bin/bash /opt/biscuit/scripts/Bisulfite_QC_Coveragestats.sh \
    ~{vcf} ~{bam} ~{reference} ~{QCannotation}
  >>>

  output {
    File bga_bed = "bga.bed"
    File cov_dist = "covdist_table.txt"
    File bga_bed_dup = "bga_dup.bed"
    File dup_report = "dup_report.txt"
    File cpg_bed = "cpg.bed"
    File cov_dist_cpg = "covdist_cpg_table.txt"
    File cpg_dist = "cpg_dist_table.txt"
  }
}
