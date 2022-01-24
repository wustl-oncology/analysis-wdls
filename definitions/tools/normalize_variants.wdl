version 1.0

task normalizeVariants {
  input {
    File reference
    File reference_fai
    File reference_dict

    File vcf
    File vcf_tbi
  }

  Int space_needed_gb = 10 + round(size([vcf, vcf_tbi], "GB") + size([reference, reference_fai, reference_dict], "GB"))
  runtime {
    memory: "9GB"
    docker: "broadinstitute/gatk:4.1.8.1"
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  command <<<
    /gatk/gatk --java-options -Xmx8g LeftAlignAndTrimVariants -O normalized.vcf.gz -R ~{reference} -V ~{vcf}
  >>>

  output {
    File normalized_vcf = "normalized.vcf.gz"
  }
}
