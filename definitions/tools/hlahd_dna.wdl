version 1.0

task hlahdDna {
  input {
    String hlahd_name = "hlahd"
    File cram
    File cram_crai
    File reference
    File reference_fai
    Int threads = 8
    Int mem = 50
  }

  Int space_needed_gb = 10 + round(5*size([cram, cram_crai, reference, reference_fai], "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "~{mem}GB"
    cpu: threads 
    docker: "jinglunli/hlahd:1.2"
    disks: "local-disk ~{space_needed_gb} HDD"
    bootDiskSizeGb: 3*space_needed_gb
  }

  command <<<
    /bin/bash /usr/bin/hlahd_script_wdl.sh /tmp . \
    ~{hlahd_name} ~{cram} ~{reference} ~{threads} ~{mem}
  >>>

  output {
    File hlahd_result_txt = hlahd_name + "_DNA/result/" + hlahd_name + "_DNA_final.result.txt"
  }
}

workflow wf {
  input {
    String? hlahd_name
    File cram
    File cram_crai
    File reference
    File reference_fai
    Int? threads
    Int? mem
  }
  call hlahdDna {
    input:
    hlahd_name=hlahd_name,
    cram=cram,
    cram_crai=cram_crai,
    reference=reference,
    reference_fai=reference_fai,
    threads=threads,
    mem=mem
  }
}
