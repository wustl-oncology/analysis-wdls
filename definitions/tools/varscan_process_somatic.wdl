version 1.0

task varscanProcessSomatic {
  input {
    File variants
    Float? max_normal_freq
  }

  Int space_needed_gb = 10 + round(size(variants, "GB")*2)
  runtime {
    memory: "4GB"
    docker: "mgibio/cle:v1.3.1"
  }

  String variantsBase = basename(variants, ".vcf")
  command <<<
    mv ~{variants} ~{basename(variants)}
    java -jar /opt/varscan/VarScan.jar processSomatic \
    ~{basename(variants)} ~{if defined(max_normal_freq) then "--max-normal-freq {max_normal_freq}" else ""}
  >>>

  output {
    File somatic_hc = variantsBase + ".Somatic.hc.vcf"
    File somatic = variantsBase + ".Somatic.vcf"
    File germline_hc = variantsBase + ".Germline.hc.vcf"
    File germline = variantsBase + ".Germline.vcf"
    File loh_hc = variantsBase + ".LOH.hc.vcf"
    File loh = variantsBase + ".LOH.vcf"
  }
}

workflow wf {
  input {
    File variants
    Float? max_normal_freq
  }

  call varscanProcessSomatic {
    input: variants=variants, max_normal_freq=max_normal_freq
  }
}
