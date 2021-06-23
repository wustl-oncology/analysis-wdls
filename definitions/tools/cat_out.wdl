version 1.0

task catOut {
  input {
    Array[File] pindel_outs
  }

  runtime {
    memory: "4GB"
    docker: "ubuntu:xenial"
  }

  command <<<
    /bin/cat ~{sep=" " pindel_outs} > "per_chromosome_pindel.out"
  >>>

  output {
    File pindel_out = "per_chromosome_pindel.out"
  }
}
