version 1.0



task concordance {
  input {
    File vcf
    File reference
    File reference_fai
    File reference_dict

    Array[File] bams
    Array[File] bais
  }

  Int space_needed_gb = 10 + round(size(flatten([[vcf, reference, reference_fai, reference_dict], bams, bais]), "GB"))
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
    set -eou pipefail

    for bam in ~{sep=' ' bams}; do
      /usr/bin/somalier extract -d extracted/ -s ~{vcf} -f ~{reference} "$bam"
    done

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
  
  Array[File] bams = select_all([bam_1, bam_2, bam_3])
  Array[File] bais = select_all([bam_1_bai, bam_2_bai, bam_3_bai])

  call concordance {
    input:
    vcf=vcf,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bams = bams,
    bais = bais
  }
}