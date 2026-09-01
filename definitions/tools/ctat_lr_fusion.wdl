version 1.0

task ctat_lr_fusion {
  input {
    File lr_bam
    # Required. Cromwell localizes each File independently, so the index is not
    # guaranteed to land beside the bam -- the command block below moves both to
    # the working directory under their basenames so ctat-LR-fusion can do
    # random access. Without this the run fails on a missing index.
    File lr_bam_bai
    File star_fusion_genome_lib_zip
    Int cpu = 30
    Boolean examine_coding_effect = true
    Boolean vis = true
  }

  Int space_needed_gb = 50 + round(size([lr_bam, lr_bam_bai, star_fusion_genome_lib_zip], "GB") * 3)
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

    # Put the bam and its index next to each other under their plain basenames.
    # mv rather than cp: same filesystem, so it is free, and it avoids doubling
    # disk for a large alignment. This mirrors tools/bam_readcount.wdl and
    # tools/strelka.wdl, which do the same thing for the same reason.
    mv ~{lr_bam} ~{basename(lr_bam)}
    mv ~{lr_bam_bai} ~{basename(lr_bam_bai)}

    # Unzip the genome library
    mkdir -p genome_lib
    unzip -qq ~{star_fusion_genome_lib_zip} -d genome_lib

    # Run CTAT LR Fusion
    ctat-LR-fusion \
      --LR_bam ~{basename(lr_bam)} \
      --genome_lib_dir genome_lib \
      --CPU ~{cpu} \
      ~{true="--vis" false="" vis} \
      ~{true="--examine_coding_effect" false="" examine_coding_effect}

    # WDL 1.0 has no portable Directory type (see docs/common_errors.md), so
    # the full output tree crosses the task boundary as an archive rather
    # than a Directory output. tar rather than zip: unlike tools/agfusion.wdl
    # and tools/pvacfuse.wdl, this container isn't known to ship `zip`, and
    # tar is present on essentially every Linux image without extra deps.
    tar czf ctat_LR_fusion_outdir.tar.gz -C ctat_LR_fusion_outdir .
  >>>

  output {
    File fusion_predictions = "ctat_LR_fusion_outdir/ctat-LR-fusion.fusion_predictions.tsv"
    File fusion_predictions_abridged = "ctat_LR_fusion_outdir/ctat-LR-fusion.fusion_predictions.abridged.tsv"
    File ctat_lr_fusion_outdir_tar_gz = "ctat_LR_fusion_outdir.tar.gz"
    # Only written when vis=true (the default). Optional so a vis=false run
    # doesn't fail on a missing glob.
    File? fusion_inspector_web_html = "ctat_LR_fusion_outdir/ctat-LR-fusion.fusion_inspector_web.html"
  }
}
