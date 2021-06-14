version 1.0

task indexVcf {
  input {
    File vcf
  }

  Int space_needed_gb = 10 + round(2*size(vcf, "GB"))
  runtime {
    docker: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
    memory: "4GB"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  # TODO: check how to run locally
  command <<<
    cp ~{vcf} "test.vcf.gz"
    /usr/local/bin/tabix -p vcf "test.vcf.gz"
  >>>
  output {
    File indexed_vcf = "test.vcf.gz"
    File indexed_vcf_tbi = "test.vcf.gz.tbi"
  }
}

workflow wf {
  input { File vcf }
  call indexVcf {
    input: vcf=vcf
  }
}
