version 1.0

task ctat_lr_fusion {
  input {
    File lr_bam
    File star_fusion_genome_lib_zip
    Int cpu = 30
    Boolean examine_coding_effect = true
    Boolean vis = true
  }

  Int space_needed_gb = 50 + round(size([lr_bam, star_fusion_genome_lib_zip], "GB") * 3)
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "trinityctat/ctat_lr_fusion:1.2.1"
    memory: "300GB"
    cpu: cpu
    disks: "local-disk ~{space_needed_gb} HDD"
    bootDiskSizeGb: 50
  }

  command <<<
    set -eou pipefail

    # Unzip the genome library
    mkdir -p genome_lib
    unzip -qq ~{star_fusion_genome_lib_zip} -d genome_lib

    # Run CTAT LR Fusion
    ctat-LR-fusion \
      --LR_bam ~{lr_bam} \
      --genome_lib_dir genome_lib \
      --CPU ~{cpu} \
      ~{true="--vis" false="" vis} \
      ~{true="--examine_coding_effect" false="" examine_coding_effect}
  >>>

  output {
    File fusion_predictions = "ctat_LR_fusion_outdir/ctat-LR-fusion.fusion_predictions.tsv"
    File fusion_predictions_abridged = "ctat_LR_fusion_outdir/ctat-LR-fusion.fusion_predictions.abridged.tsv"
    Directory ctat_lr_fusion_outdir = "ctat_LR_fusion_outdir"
  }
}
