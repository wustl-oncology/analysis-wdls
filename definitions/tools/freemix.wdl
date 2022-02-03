version 1.0

task freemix {
  input {
    File verify_bam_id_metrics
  }

  runtime {
    docker: "python:3.10"
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

workflow wf {
  input {
    File verify_bam_id_metrics
  }

  call freemix {
    input:
    verify_bam_id_metrics=verify_bam_id_metrics
  }
}
