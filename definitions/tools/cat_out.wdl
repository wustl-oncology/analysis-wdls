version 1.0

task catOut {
  input {
    Array[File] pindel_outs
  }

  Int space_needed_gb = 10 + round(size(pindel_outs, "GB")*2)
  runtime {
    memory: "4GB"
    docker: "ubuntu:xenial"
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  command <<<
    /bin/cat ~{sep=" " pindel_outs} > "per_chromosome_pindel.out"
  >>>

  output {
    File pindel_out = "per_chromosome_pindel.out"
  }
}
