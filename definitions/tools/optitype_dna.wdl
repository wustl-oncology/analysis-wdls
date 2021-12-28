version 1.0

task optitypeDna {
  input {
    String optitype_name = "optitype"
    File cram
    File cram_crai
    File reference
    File reference_fai
  }

  Int space_needed_gb = 10 + round(5*size([cram, cram_crai, reference, reference_fai], "GB"))
  runtime {
    memory: "64GB"
    docker: "mgibio/immuno_tools-cwl:1.0.1"
    disks: "local-disk ~{space_needed_gb} SSD"
    bootDiskSizeGb: 3*space_needed_gb
  }

  command <<<
    /bin/bash /usr/bin/optitype_script.sh /tmp . \
    ~{optitype_name} ~{cram} ~{reference}
  >>>

  output {
    File optitype_tsv = optitype_name + "_result.tsv"
    File optitype_plot = optitype_name + "_coverage_plot.pdf"
  }
}

workflow wf {
  input {
    String? optitype_name
    File cram
    File cram_crai
    File reference
    File reference_fai
  }
  call optitypeDna {
    input:
    optitype_name=optitype_name,
    cram=cram,
    cram_crai=cram_crai,
    reference=reference,
    reference_fai=reference_fai,
  }
}
