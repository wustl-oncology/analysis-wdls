version 1.0

task mergeVcf {
  input {
    Array[File] vcfs
    Array[File] vcf_tbis
    String merged_vcf_basename = "merged"
  }

  runtime {
    docker: "mgibio/bcftools-cwl:1.3.1"
    memory: "4GB"
  }

  String output_file = merged_vcf_basename + ".vcf.gz"
  command <<<
    /opt/bcftools/bin/bcftools concat --allow-overlaps --remove-duplicates --output-type z -o ~{output_file} ~{sep=" " vcfs}
  >>>

  output {
    File merged_vcf = output_file
  }
}
