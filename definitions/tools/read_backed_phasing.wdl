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
  runtime {
      docker: "mgibio/gatk-cwl:3.6.0"
      memory: "9GB"
      bootDiskSizeGb: 25
  }
  command <<<
  /usr/bin/java -Xmx8g -jar /opt/GenomeAnalysisTK.jar -T ReadBackedPhasing -L ~{vcf} -o phased.vcf
  >>>
  output {
    File phased_vcf = glob("phased.vcf")
  }
}
