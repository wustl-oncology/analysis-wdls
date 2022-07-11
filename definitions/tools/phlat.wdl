version 1.0

task phlat {
  input {
    String name = "phlat"
    File cram
    File cram_crai
    File reference
    File reference_fai
   }

  runtime {
    memory: "20GB"
    docker: "laljorani/phlat:latest"
    disks: ""
    bootDiskSizeGb:
  }

  command <<<
    /bin/bash /usr/bin/run.b38.sh \
    --phlat-dir ~()               \
    --data-dir ~()                \
    --tag ~()                     \
    --bam ~()                     \
    --index-dir ~()               \
    --rs-dir ~()                  \
    --b2url ~()                   \
    --fastq1 ~()                  \
    --fastq2 ~()                  \
    --ref-fasta ~() 
  >>>

  output {
  }
}

workflow wf {
  input {
   }

  call phlat{
  }
} 
