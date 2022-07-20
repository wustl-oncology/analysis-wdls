version 1.0

task phlat {
  input {
    String phlat_name = "phlat" 
    File cram
    File reference
    Int nthreads = 8
    Int mem = 20
    String index_dir = "" # optional if indexes are within run.b38.sh default dir
   }

  Int space_needed_gb = 10 + round(5*size([cram, cram_crai, reference, reference_fai], "GB"))
  runtime {
    memory: "~{mem}GB"
    cpu: nthreads
    docker: "mgibio/phlat:1.1_withindex"
    disks: "local-disk ~{space_needed_gb} HDD"
    bootDiskSizeGb: 3*space_needed_gb
  }

  command <<<
    /bin/bash /usr/bin/run.b38.sh       \
    --tag ~{phlat_name} --bam ~{cram}   \ 
    --ref-fasta ~{reference}            \
    --index-dir ~{index_dir}
  >>>

  output {
    File phlat_summary = "/usr/bin/phlat-release/example/results/${phlat_name}_result.sum"
  }
}

workflow wf {
  input {
    String? phlat_name
    File cram
    File reference
    Int? nthreads
    Int? mem
    String? index_dir 
  }
  call phlat {
    input:
    phlat_name=phlat_name,
    cram=cram,
    reference=reference,
    nthreads=nthreads,
    mem=mem,
    index_dir=index_dir,
  }
} 
