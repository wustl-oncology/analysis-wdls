version 1.0

task filterVcfMapq0 {
  input {
    File vcf
    File tumor_bam
    File tumor_bam_bai
    File reference
    File reference_fai
    File reference_dict
    Float threshold
  }

  Float reference_size = size([reference, reference_fai, reference_dict], "GB")
  Float bam_size = size([tumor_bam, tumor_bam_bai], "GB")
  Int space_needed_gb = 10 + round(reference_size + bam_size + 2*size(vcf, "GB"))
  runtime {
    docker: "mgibio/mapq0-filter:v0.3.1"
    memory: "8GB"
    bootDiskSizeGb: 10
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  String outfile = "mapq_filtered.vcf.gz"
  command <<<
    /bin/bash /usr/bin/mapq0_vcf_filter.sh ~{outfile} ~{vcf} ~{tumor_bam} ~{reference} ~{threshold}
  >>>

  output {
    File mapq0_filtered_vcf = outfile
    File mapq0_filtered_vcf_tbi = outfile + ".tbi"
  }
}

workflow wf { call filterVcfMapq0 { input: } }
