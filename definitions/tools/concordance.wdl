version 1.0

task concordance {
  input {
    File vcf

    File reference
    File reference_fai
    File reference_dict

    File bam_1
    File bam_1_bai
    File bam_2
    File bam_2_bai
    File? bam_3
    File? bam_3_bai
  }

  Int space_needed_gb = 10 + round(size([vcf, reference, reference_fai, reference_dict, bam_1, bam_1_bai, bam2, bam_2_bai, bam_3, bam_3_bai], "GB"))
  runtime {
    cpu: 1
    memory: "8GB"
    bootDiskSizeGb: 10
    docker: "brentp/somalier:v0.1.5"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    /usr/bin/somalier -s ~{vcf} -f ~{reference} ~{bam1} ~{bam2} ~{bam3}
  >>>

  output {
    File somalier_pairs = "concordance.somalier.pairs.tsv"
    File somalier_samples = "concordance.somalier.samples.tsv"
  }
}
