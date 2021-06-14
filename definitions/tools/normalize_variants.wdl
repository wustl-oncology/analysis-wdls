version 1.0

task normalizeVariants {
  input {
    File reference
    File reference_fai
    File reference_dict

    File vcf
    File vcf_tbi
  }

  runtime {
    memory: "9GB"
    docker: "broadinstitute/gatk:4.1.8.1"
  }

  command <<<
    /gatk/gatk --java-options -Xmx8g LeftAlignAndTrimVariants -O normalized.vcf.gz -R ~{reference} -V ~{vcf}
  >>>

  output {
    File normalized_vcf = "normalized.vcf.gz"
    File normalized_vcf_tbi = "normalized.vcf.gz.tbi"
  }
}
