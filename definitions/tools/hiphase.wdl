version 1.0

# HiPhase: read-backed phasing for PacBio HiFi.
#
# This is the long-read replacement for GATK3 ReadBackedPhasing in
# subworkflows/phase_vcf.wdl, and it is strictly better for the pVACseq use
# case: phase blocks span the full read length (15-25 kb) rather than a
# short-read insert (~300-500 bp), so far more germline/somatic variant pairs
# get a usable phase relationship.
#
# hiphaseSomatic phases TWO vcfs against the same bam in one pass, so the
# germline and somatic calls come out sharing the same PS phase-set ids. That
# shared PS is what makes them safe to concatenate into a single proximal
# variants vcf downstream (see tools/combine_phased_proximal_vcf.wdl).

task hiphaseSomatic {
  input {
    File bam
    File bam_bai
    File germline_vcf
    File germline_vcf_tbi
    File somatic_vcf
    File somatic_vcf_tbi
    File reference
    File reference_fai
    String output_basename
    Int cores = 24
  }

  Int space_needed_gb = 20 + round(3*size(bam, "GB") + size(reference, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "quay.io/pacbio/hiphase@sha256:353b4ffdae4281bdd5daf5a73ea3bb26ea742ef2c36e9980cb1f1ed524a07482"
    memory: "~{cores * 4}GB"
    cpu: cores
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    set -euxo pipefail
    hiphase --version

    hiphase \
      --bam ~{bam} \
      -t ~{cores} \
      --output-bam ~{output_basename}.hiphase.bam \
      --vcf ~{germline_vcf} \
      --output-vcf ~{output_basename}.germline.hiphase.vcf.gz \
      --vcf ~{somatic_vcf} \
      --output-vcf ~{output_basename}.somatic.hiphase.vcf.gz \
      -r ~{reference} \
      --stats-file ~{output_basename}.hiphase.stats.csv \
      --summary-file ~{output_basename}.hiphase.summary.tsv \
      --blocks-file ~{output_basename}.hiphase.blocks.tsv \
      --ignore-read-groups
  >>>

  output {
    File phased_bam = output_basename + ".hiphase.bam"
    File phased_bam_bai = output_basename + ".hiphase.bam.bai"
    File phased_germline_vcf = output_basename + ".germline.hiphase.vcf.gz"
    File phased_germline_vcf_tbi = output_basename + ".germline.hiphase.vcf.gz.tbi"
    File phased_somatic_vcf = output_basename + ".somatic.hiphase.vcf.gz"
    File phased_somatic_vcf_tbi = output_basename + ".somatic.hiphase.vcf.gz.tbi"
    File stats = output_basename + ".hiphase.stats.csv"
    File summary = output_basename + ".hiphase.summary.tsv"
    File blocks = output_basename + ".hiphase.blocks.tsv"
  }
}

# Germline-only phasing, used on the normal bam.
task hiphaseGermline {
  input {
    File bam
    File bam_bai
    File vcf
    File vcf_tbi
    File reference
    File reference_fai
    String output_basename
    Int cores = 24
  }

  Int space_needed_gb = 20 + round(3*size(bam, "GB") + size(reference, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "quay.io/pacbio/hiphase@sha256:353b4ffdae4281bdd5daf5a73ea3bb26ea742ef2c36e9980cb1f1ed524a07482"
    memory: "~{cores * 4}GB"
    cpu: cores
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    set -euxo pipefail
    hiphase --version

    hiphase \
      --bam ~{bam} \
      -t ~{cores} \
      --output-bam ~{output_basename}.hiphase.bam \
      --vcf ~{vcf} \
      --output-vcf ~{output_basename}.hiphase.vcf.gz \
      -r ~{reference} \
      --stats-file ~{output_basename}.hiphase.stats.csv \
      --summary-file ~{output_basename}.hiphase.summary.tsv \
      --blocks-file ~{output_basename}.hiphase.blocks.tsv \
      --ignore-read-groups
  >>>

  output {
    File phased_bam = output_basename + ".hiphase.bam"
    File phased_bam_bai = output_basename + ".hiphase.bam.bai"
    File phased_vcf = output_basename + ".hiphase.vcf.gz"
    File phased_vcf_tbi = output_basename + ".hiphase.vcf.gz.tbi"
    File stats = output_basename + ".hiphase.stats.csv"
    File summary = output_basename + ".hiphase.summary.tsv"
    File blocks = output_basename + ".hiphase.blocks.tsv"
  }
}
