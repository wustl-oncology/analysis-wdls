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
    /usr/bin/awk '{FS="\t";getline;for(n=2;n<=NF-2;n++){if($n==""){}else{printf "HLA-"$n"\n"}}}' ~{file} > ~{outname}
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
