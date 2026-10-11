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
    String docker_image = "jinglunli/pvacnc:0.3.1"
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
    pvacnc-prepare \
      ~{sample_name} ~{tumor_sample_name} ~{source_vcf} ~{reference} \
      ~{reference_fai} ~{reference_annotation} ~{with_e_orfanage_gtf} \
      ~{no_e_orfanage_gtf}

    # Apply the caller-specific publication prefix after preparation finishes.
    for directory in with_e/results no_e/results results pvacseq_input; do
      for file in "$directory/~{sample_name}."*; do
        [ -f "$file" ] || continue
        name=$(basename "$file")
        suffix=${name#~{sample_name}.}
        suffix=${suffix#orfanage_}
        mv "$file" "$directory/~{sample_name}.orfanage.$suffix"
      done
    done
  >>>

  output {
    File noncanonical_orf_blocks = "results/~{sample_name}.orfanage.noncanonical_orf_blocks.tsv"
    File pvacnc_input_vcf_gz = "pvacseq_input/~{sample_name}.orfanage.pvacnc_input.vcf.gz"
    File pvacnc_input_vcf_tbi = "pvacseq_input/~{sample_name}.orfanage.pvacnc_input.vcf.gz.tbi"
    File mapping_tsv = "pvacseq_input/~{sample_name}.orfanage.pvacnc_mapping.tsv"
    File skipped_candidates_tsv = "pvacseq_input/~{sample_name}.orfanage.pvacnc_skipped_candidates.tsv"
    File routing_tsv = "results/~{sample_name}.orfanage.e_first_no_e_rescue_routing.tsv"
    File effects_tsv = "results/~{sample_name}.orfanage.e_first_no_e_rescue_orfanage_effects.tsv"
    Array[File] preparation_outputs = flatten([
      glob("with_e/results/*"),
      glob("no_e/results/*"),
      glob("results/*"),
      glob("pvacseq_input/*")
    ])
  }
}
