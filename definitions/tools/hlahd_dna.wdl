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
    docker: "griffithlab/hlahd:1.0"
    disks: "local-disk ~{space_needed_gb} HDD"
    bootDiskSizeGb: 3*space_needed_gb
  }

  command <<<
    # hlahd_script_wdl.sh internally runs Picard SamToFastq to convert the 
    # input CRAM before HLA typing, but it never applies the `mem` argument to 
    # that Java invocation - it defaults to Java 8's unconfigured heap (~4GB)
    # regardless of how much memory this task actually requested/received.
    #
    # Workaround: _JAVA_OPTIONS is honored by every JVM launched in this shell,
    # so appending -Xmx here forces SamToFastq to actually use the memory reserved for this task.
    # Leave ~10GB of headroom for hlahd/optitype's own (non-JVM) memory use after the conversion step.
    export _JAVA_OPTIONS="${_JAVA_OPTIONS} -Xmx~{mem - 10}g"

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
