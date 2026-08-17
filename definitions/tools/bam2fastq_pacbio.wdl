version 1.0

# Converts an unaligned PacBio bam (e.g. FLNC reads from `isoseq refine`) to
# fastq using PacBio's own pbtk toolkit, as a replacement for `samtools fastq`.
#
# Why this exists: `samtools fastq` routes any read carrying the
# paired-in-sequencing SAM flag (0x1) to its singleton/`-s` output rather than
# the main output, even when the data is biologically single-end -- the
# routing is flag-driven, not data-type-driven. tools/minimap2_rnaseq_pacbioHifi.wdl's
# bamToFastqLr sends `-s` to /dev/null on the (usually safe) assumption that
# FLNC reads are never paired. When that assumption is wrong for a given bam,
# every read is silently discarded: samtools reports "processed N reads" on
# stderr while the fastq it writes is 0 bytes, and nothing downstream notices
# because the task still exits 0 -- the failure only surfaces much later, as
# an empty alignment feeding ctat-LR-fusion.
#
# bam2fastq is PacBio's own tool for exactly this conversion. It converts
# every record in the bam with no flag-based routing, so it isn't exposed to
# the same failure mode. It requires a .pbi index (built with pbindex)
# sitting next to the input bam, which this task builds locally since
# Cromwell doesn't guarantee the localized input path is writable or that the
# index travels with it.
#
# Docker and both commands proven manually against a real PacBio Kinnex FLNC
# bam in human/long_read_analysis/bam2fastq/B627Ctopup/ before being wired in
# here -- see build_index.sh / bam2fastq.sh / bam2fastq_job.sh there. No
# thread flags are passed to either tool because none were used in that
# proven invocation.
task bam2fastqPacbio {
  input {
    File flnc_bam
    String output_basename = "flnc"
    Int cores = 4
  }

  # bam2fastq output is gzip-compressed, unlike samtools fastq's plain-text
  # output, so a smaller multiplier than bamToFastqLr's is enough headroom.
  Int space_needed_gb = 10 + round(3 * size(flnc_bam, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "quay.io/biocontainers/pbtk:3.5.0--h9ee0642_0"
    memory: "8GB"
    cpu: cores
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  # pbindex writes <bam>.pbi next to its input, so the bam is moved into the
  # task's own cwd first -- Cromwell's localized input path is not
  # guaranteed to be writable. Same reasoning as tools/ctat_lr_fusion.wdl's
  # bam+bai relocation.
  String local_bam = basename(flnc_bam)
  command <<<
    set -eou pipefail

    mv ~{flnc_bam} ~{local_bam}
    pbindex ~{local_bam}
    bam2fastq -o ~{output_basename} ~{local_bam}
  >>>

  output {
    File fastq = output_basename + ".fastq.gz"
  }
}
