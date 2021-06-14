version 1.0

task vtDecompose {
  input {
    File vcf
    File vcf_tbi
  }

  runtime {
    memory: "4GB"
    docker: "quay.io/biocontainers/vt:0.57721--hf74b74d_1"
  }

  command <<<
    vt decompose -s -o decomposed.vcf.gz ~{vcf}
  >>>

  output {
    File decomposed_vcf = "decomposed.vcf.gz"
  }
}
