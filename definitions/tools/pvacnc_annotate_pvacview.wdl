version 1.0

task pvacncAnnotatePvacview {
  input {
    File mapping_tsv
    Array[File] source_mhc_i
    Array[File] source_mhc_ii
    Array[File] source_combined
    String output_prefix = "pvacnc"
    String docker_image = "jinglunli/pvacnc:latest"
  }

  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "4GB"
    cpu: 1
    docker: docker_image
    disks: "local-disk 20 HDD"
  }

  command <<<
    set -e
    mkdir -p pvacseq_predictions/MHC_Class_I pvacseq_predictions/MHC_Class_II pvacseq_predictions/combined

    for input in ~{sep=' ' source_mhc_i}; do
      cp -L "$input" "pvacseq_predictions/MHC_Class_I/$(basename "$input")"
    done
    for input in ~{sep=' ' source_mhc_ii}; do
      cp -L "$input" "pvacseq_predictions/MHC_Class_II/$(basename "$input")"
    done
    for input in ~{sep=' ' source_combined}; do
      cp -L "$input" "pvacseq_predictions/combined/$(basename "$input")"
    done

    pvacnc-annotate-pvacview \
      --mapping ~{mapping_tsv} \
      --directory pvacseq_predictions/MHC_Class_I \
      --directory pvacseq_predictions/MHC_Class_II \
      --directory pvacseq_predictions/combined \
      > ~{output_prefix}.pvacview_annotation.log
  >>>

  output {
    Array[File] mhc_i = glob("pvacseq_predictions/MHC_Class_I/*")
    Array[File] mhc_ii = glob("pvacseq_predictions/MHC_Class_II/*")
    Array[File] combined = glob("pvacseq_predictions/combined/*")
    File annotation_log = "~{output_prefix}.pvacview_annotation.log"
  }
}
