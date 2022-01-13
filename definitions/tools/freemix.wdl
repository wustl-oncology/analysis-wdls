version 1.0

task freemix {
  input {
    File verify_bam_id_metrics
  }

  runtime {
    docker: "ubuntu:xenial"
  }

  command <<<
    python <<CODE
    with open("~{verify_bam_id_metrics}", "r") as f:
        header = f.readline().split("\t")
        if len(header) >= 7 and header[6] == "FREEMIX":
            print(f.readline().split("\t")[6])
    CODE
  >>>

  output {
    # Contains EITHER the FREEMIX contamination fraction value, or is empty
    # WDL doesn't have a good way to extract `Float?` from a command...
    String out = stdout()
  }
}
