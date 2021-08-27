version 1.0

task optitypeDna {
  input {
    String optitype_name = "optitype"
    File cram
    File reference
    File reference_fai
  }

  Int space_needed_gb = 10 + round(size([cram, reference, reference_fai], "GB"))
  runtime {
    memory: "64GB"
    bootDiskSizeGb: 20
    docker: "mgibio/immuno_tools-cwl:1.0.1"
    disks: "local-disk ~{space_needed_gb} SSD"
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
