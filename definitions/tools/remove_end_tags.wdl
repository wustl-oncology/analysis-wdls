version 1.0

task removeEndTags {
  input {
    File vcf
    File vcf_tbi
  }

  runtime {
    memory: "4GB"
    docker: "mgibio/bcftools-cwl:1.3.1"
  }

  String outfile = "pindel.noend.vcf.gz"
  command <<<
    /opt/bcftools/bin/bcftools annotate -x INFO/END -Oz -o ~{outfile} ~{vcf}
  >>>

  output {
    File processed_vcf = outfile
  }
}

workflow wf {
  input {
    File vcf
    File vcf_tbi
  }

  call removeEndTags {
    input:
    vcf=vcf,
    vcf_tbi=vcf_tbi
  }
}
