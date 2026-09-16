version 1.0

task orfanageBest {
  input {
    File stringtie_gtf
    File reference_annotation
    File reference
    String output_prefix
    String docker_image = "quay.io/biocontainers/orfanage:1.2.0--heaafb18_2"
  }

  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "16GB"
    cpu: 4
    docker: docker_image
    disks: "local-disk 50 HDD"
  }

  command <<<
    orfanage --reference ~{reference} --query ~{stringtie_gtf} \
      --output ~{output_prefix}.gtf --mode BEST ~{reference_annotation}
  >>>

  output {
    File best_gtf = "~{output_prefix}.gtf"
  }
}
