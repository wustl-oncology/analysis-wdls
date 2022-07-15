version 1.0

task phlat {
  input {
    String phlat_name = "phlat"
    File cram
    File cram_crai
    File reference
    File reference_fai
    String name = basename(cram, ".*")
    String index_dir = "" # optional if indexes are within run.b38.sh default dir
   }

  Int space_needed_gb = 10 + round(5*size([cram, cram_crai, reference, reference_fai], "GB"))
  runtime {
    memory: "20GB"
    docker: "laljorani/phlat:latest"
    disks: "local-disk ~{space_needed_gb} HDD"
    bootDiskSizeGb: 3*space_needed_gb
  }

  command <<<
    /bin/bash /usr/bin/run.b38.sh \
    --tag ~(name)                 \
    --bam ~(cram)                 \ 
    --ref-fasta ~(reference)      \
    --index-dir ~(index_dir)
  >>>

  output {
    ...
  }
}

workflow wf {
  input {
    String? phlat
    File cram
    File cram_crai
    File reference
    File reference_fai
  }
  call optitypeDna {
    input:
    phlat_name=phlat_name,
    cram=cram,
    cram_crai=cram_crai,
    reference=reference,
    reference_fai=reference_fai,
  }
} 
