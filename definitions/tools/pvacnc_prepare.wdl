version 1.0

task pvacncPrepare {
  input {
    String sample_name
    String tumor_sample_name
    File source_vcf
    File reference
    File reference_fai
    File reference_annotation
    File with_e_orfanage_gtf
    File no_e_orfanage_gtf
    Array[File] canonical_pvacseq_outputs
    String docker_image = "jinglunli/pvacnc:0.2.2"
  }

  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "32GB"
    cpu: 4
    docker: docker_image
    disks: "local-disk 100 HDD"
  }

  command <<<
    set -e
    mkdir -p canonical_outputs
    for output in ~{sep=' ' canonical_pvacseq_outputs}; do
      case "$(basename "$output")" in
        *Combined.all_epitopes.tsv)
          cp -L "$output" "canonical_outputs/$(basename "$output")"
          ;;
      esac
    done
    test -f canonical_outputs/*Combined.all_epitopes.tsv
    pvacnc-prepare \
      ~{sample_name} ~{tumor_sample_name} ~{source_vcf} ~{reference} \
      ~{reference_fai} ~{reference_annotation} ~{with_e_orfanage_gtf} \
      ~{no_e_orfanage_gtf} canonical_outputs
  >>>

  output {
    File pvacnc_input_vcf_gz = "pvacseq_input/~{sample_name}.pvacnc_input.vcf.gz"
    File pvacnc_input_vcf_tbi = "pvacseq_input/~{sample_name}.pvacnc_input.vcf.gz.tbi"
    File mapping_tsv = "pvacseq_input/~{sample_name}.pvacnc_mapping.tsv"
    File routing_tsv = "results/~{sample_name}.e_first_no_e_rescue_routing.tsv"
    File effects_tsv = "results/~{sample_name}.e_first_no_e_rescue_orfanage_effects.tsv"
    Array[File] preparation_outputs = flatten([
      glob("with_e/results/*"),
      glob("no_e/results/*"),
      glob("results/*"),
      glob("pvacseq_input/*")
    ])
  }
}
