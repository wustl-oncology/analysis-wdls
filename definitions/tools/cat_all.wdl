version 1.0

task catAll {
  input {
    Array[File] region_pindel_outs
  }

  runtime {
    memory: "4GB"
    docker: "ubuntu:xenial"
  }

  command <<<
    /bin/cat ~{sep=" " region_pindel_outs} | /bin/grep "ChrID" /dev/stdin > all_region_pindel.head
  >>>

  output {
    File all_region_pindel_head = "all_region_pindel.head"
  }
}
