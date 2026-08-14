version 1.0

# Coverage summary for HiFi alignments.
#
# There is no bait/target interval list on the WGS arm, so Picard CollectHsMetrics
# does not apply. mosdepth gives the per-base depth bed and summary that you need
# to (a) sanity check depth before picking pVACseq's --tdna-cov / --normal-cov
# thresholds, and (b) feed a TMB calculation later if you want one.
task mosdepth {
  input {
    File bam
    File bam_bai
    String output_basename
    Int cores = 4
    # 0 disables per-base output and is much faster; leave at 0 unless you
    # actually need the per-base bed.
    Int window_size = 0
  }

  Int space_needed_gb = 10 + round(1.5*size(bam, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "quay.io/biocontainers/mosdepth:0.3.4--hd299d5a_0"
    memory: "16GB"
    cpu: cores
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    set -euxo pipefail
    mosdepth --version

    mosdepth \
      -t ~{cores} \
      ~{if window_size > 0 then "--by " + window_size else ""} \
      ~{output_basename} \
      ~{bam}
  >>>

  output {
    File summary = output_basename + ".mosdepth.summary.txt"
    File global_dist = output_basename + ".mosdepth.global.dist.txt"
    File? regions_bed = output_basename + ".regions.bed.gz"
  }
}

workflow wf {
  input {
    File bam
    File bam_bai
    String output_basename
  }
  call mosdepth {
    input: bam=bam, bam_bai=bam_bai, output_basename=output_basename
  }
}
