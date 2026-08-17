version 1.0

task optitypeDna {
  input {
    String optitype_name = "optitype"
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
    docker: "mgibio/immuno_tools-cwl:1.0.2"
    disks: "local-disk ~{space_needed_gb} HDD"
    bootDiskSizeGb: 3*space_needed_gb
  }

  command <<<
    # optitype_script_wdl.sh (baked into the mgibio/immuno_tools-cwl image) internally
    # runs Picard SamToFastq to convert the input CRAM before HLA typing, but it never
    # applies the `mem` argument to that Java invocation - it defaults to Java 8's
    # unconfigured heap (~4GB) regardless of how much memory this task actually
    # requested/received. On larger CRAMs this causes SamToFastq to die with
    # "OutOfMemoryError: GC overhead limit exceeded" even though the LSF/backend job
    # itself had ~{mem}GB available.
    # See: https://github.com/wustl-oncology/analysis-wdls/issues/225
    #
    # Workaround: _JAVA_OPTIONS is honored by every JVM launched in this shell (Cromwell
    # already relies on this for -Djava.io.tmpdir, exported above by Cromwell's generated
    # script), so appending -Xmx here forces SamToFastq - and any other java call inside
    # the wrapper script - to actually use the memory reserved for this task. Append
    # rather than overwrite: _JAVA_OPTIONS already holds Cromwell's -Djava.io.tmpdir
    # setting, which must be preserved so temp files still land on the disk provisioned
    # for this task rather than falling back to the container's default /tmp. Leave
    # ~10GB of headroom for optitype's own (non-JVM) memory use after the conversion
    # step.
    export _JAVA_OPTIONS="${_JAVA_OPTIONS} -Xmx~{mem - 10}g"

    /bin/bash /usr/bin/optitype_script_wdl.sh /tmp . \
    ~{optitype_name} ~{cram} ~{reference} ~{threads} ~{mem}
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
    Int? threads
    Int? mem
  }
  call optitypeDna {
    input:
    optitype_name=optitype_name,
    cram=cram,
    cram_crai=cram_crai,
    reference=reference,
    reference_fai=reference_fai,
    threads=threads,
    mem=mem
  }
}
