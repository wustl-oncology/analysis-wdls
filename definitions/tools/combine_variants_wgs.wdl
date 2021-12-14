version 1.0

task combineVariantsWgs {
  input {
    File reference
    File reference_fai
    File reference_dict
    File mutect_vcf
    File mutect_vcf_tbi
    File varscan_vcf
    File varscan_vcf_tbi
    File strelka_vcf
    File strelka_vcf_tbi
  }

  String outfile = "combined.vcf.gz"
  command <<<
    /usr/bin/java -Xmx8g -jar /opt/GenomeAnalysisTK.jar -T CombineVariants \
    -genotypeMergeOptions PRIORITIZE \
    --rod_priority_list mutect,varscan,strelka \
    -o ~{outfile} \
    -R ~{reference} \
    --variant:mutect ~{mutect_vcf} \
    --variant:varscan ~{varscan_vcf} \
    --variant:strelka_vcf ~{strelka_vcf}
  >>>

  output {
    File combined_vcf = outfile
  }
}
