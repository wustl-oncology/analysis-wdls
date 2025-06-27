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

  Int space_needed_gb = 10 + round(size([vcf, reference, reference_fai, reference_dict, bam_1, bam_1_bai, bam_2, bam_2_bai, bam_3, bam_3_bai], "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    cpu: 1
    memory: "8GB"
    bootDiskSizeGb: 10
    docker: "brentp/somalier:v0.2.19"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    /usr/bin/somalier extract -d extracted/ -s ~{vcf} -f ~{reference} ~{bam_1} 
    /usr/bin/somalier extract -d extracted/ -s ~{vcf} -f ~{reference} ~{bam_2} 
    /usr/bin/somalier extract -d extracted/ -s ~{vcf} -f ~{reference} ~{bam_3} 
    /usr/bin/somalier relate -o concordance extracted/*.somalier

    mv concordance.pairs.tsv concordance.somalier.pairs.tsv
    mv concordance.samples.tsv concordance.somalier.samples.tsv
  >>>

  output {
    File somalier_pairs = "concordance.somalier.pairs.tsv"
    File somalier_samples = "concordance.somalier.samples.tsv"
  }
}

workflow wf {
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

  call concordance {
    input:
    vcf=vcf,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam_1=bam_1,
    bam_1_bai=bam_1_bai,
    bam_2=bam_2,
    bam_2_bai=bam_2_bai,
    bam_3=bam_3,
    bam_3_bai=bam_3_bai
  }
}
