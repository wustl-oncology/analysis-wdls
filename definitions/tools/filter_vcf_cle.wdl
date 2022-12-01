version 1.0

task filterVcfCle {
  input {
    File vcf
    Boolean filter
  }

  Int space_needed_gb = 10 + round(size(vcf, "GB")*2)
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "mgibio/cle:v1.3.1"
    memory: "4GB"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    /usr/bin/perl /usr/bin/docm_and_coding_indel_selection.pl ~{vcf} "$PWD" ~{filter}
  >>>

  output {
    File cle_filtered_vcf = "annotated_filtered.vcf"
  }
}
