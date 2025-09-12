version 1.0

task extractHlaAlleles {
  input {
    File optitype_file
    File phlat_file
    File hlahd_file 
  }

  Int space_needed_gb = 10 + round(size([optitype_file, phlat_file, hlahd_file], "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "2GB"
    docker: "ubuntu:xenial"
    disks: "local-disk ~{space_needed_gb} HDD"
  }
  
 # Steps:
  # 1) Extract HLA class I (OptiType)
  # 2) Extract HLA class II (PHLAT)
  # 3) From HLA-HD results, take only first 3 columns for DPA1/DPB1/DRB3/DRB4/DRB5
 
  String outname = "hla_calls_newline.txt"
  String temp = "temp.txt"

command <<<
  /usr/bin/awk '{FS="\t";getline;for(n=2;n<=NF-2;n++){if($n==""){}else{printf "HLA-"$n"\n"}}}' ~{optitype_file} > ~{temp}
  grep "HLA_D" ~{phlat_file} | /usr/bin/awk '{FS="\t";if($2==""){}else{printf $2"\n"};if($3==""){}else{printf $3"\n"}}' >> ~{temp}
  /usr/bin/awk -F'\t' '$1 ~ /^(DPA1|DPB1|DRB3|DRB4|DRB5)$/ {
    a=$2; b=$3;
    for(i=1;i<=2;i++){
      val = (i==1 ? a : b);
      if(val=="" || val=="-" || val ~ /^Not typed/) { next }
      sub(/^HLA-/, "", val);  # strip HLA- if present
      print val;
    }
  }' ~{hlahd_file} >> ~{temp}
  /usr/bin/awk -F":" '{print $1 ":" $2}' ~{temp} > ~{outname}
>>>

  output {
    Array[String] allele_string = read_lines(outname)
    File allele_file = outname
  }
}

workflow wf {
  input {
    File optitype_file 
    File phlat_file 
    File hlahd_file   
  }
  call extractHlaAlleles { 
    input: 
    optitype_file=optitype_file,
    phlat_file=phlat_file,
    hlahd_file=hlahd_file,
  }
}
