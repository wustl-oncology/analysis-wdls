version 1.0

task extractHlaAlleles {
  input {
    File file
  }

  Int space_needed_gb = 10 + round(size(file, "GB"))
  runtime {
    memory: "2GB"
    docker: "ubuntu:xenial"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String outname = "helper.txt"
  command <<<
    /usr/bin/awk '{FS="\t";if($2==""){}else{printf "HLA-"$2"\n"};if($3==""){}else{printf "HLA-"$3"\n"}}' ~{file} | sed '1, 2d' > ~{outname}
  >>>

  output {
    Array[String] allele_string = read_lines(outname)
    File allele_file = outname
  }
}

workflow wf {
  input { File file }
  call extractHlaAlleles { input: file=file }
}
