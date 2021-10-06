version 1.0

task extractHlaAlleles {
  input {
    File allele_file
  }

  Int space_needed_gb = 10 + round(size(allele_file, "GB"))
  runtime {
    memory: "2GB"
    docker: "ubuntu:xenial"
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  String outname = "helper.txt"
  command <<<
    /usr/bin/awk '{getline; printf "HLA-"$2 "\nHLA-"$3 "\nHLA-"$4 "\nHLA-"$5 "\nHLA-"$6 "\nHLA-"$7}' ~{allele_file} > ~{outname}
  >>>

  output {
    Array[String] allele_string = read_lines(outname)
    File outfile = outname
  }
}

workflow wf { call extractHlaAlleles { input: } }
