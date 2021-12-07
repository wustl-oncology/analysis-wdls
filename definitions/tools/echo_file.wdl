version 1.0

task echoFile {
  input {}
  runtime { docker: "ubuntu:bionic" }
  command <<<
    echo TEST > outfile.txt
  >>>
  output { File out = "outfile.txt" }
}

workflow wf {
  input {}
  call echoFile {}
  output { File out = echoFile.out }
}
