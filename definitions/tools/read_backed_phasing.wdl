version 1.0

task {
  input {
    File bam
    File bam_index
    File reference
    File reference_fai
    File reference_dict
    File vcf
    File vcf_tbi
  }

  Int space_needed_gb = 10 + round(size([bam, bam_index, reference, reference_fai, reference_dict, vcf, vcf_tbi], "GB"))
  runtime {
    docker: "mgibio/gatk-cwl:3.6.0"
    memory: "9GB"
    bootDiskSizeGb: 25
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    /usr/bin/java -Xmx8g -jar /opt/GenomeAnalysisTK.jar -T ReadBackedPhasing -L ~{vcf} -o phased.vcf
  >>>

  output {
    File phased_vcf = glob("phased.vcf")
  }
}
