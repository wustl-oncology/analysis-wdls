version 1.0

# Add an XS:A strand tag to a minimap2 spliced bam, derived from minimap2's ts:A.
#
# WHY THIS EXISTS
# tools/regtools.wdl maps strand "unstranded" -> "XS" and passes `-s XS`, which
# tells regtools to read the XS:A tag off each alignment to decide junction
# strand. STAR provides that tag because tools/star_align_fusion.wdl runs with
# `--outSAMstrandField intronMotif`.
#
# minimap2 does not emit XS. In splice mode it emits `ts:A:+/-`, documented as
# "Transcript strand (splice mode only)", and there is no flag to emit XS
# instead. Without this conversion regtools finds no strand information on
# minimap2 junctions.
#
# THE SIGN FLIP
# ts is the transcript strand RELATIVE TO THE READ; XS is relative to the
# REFERENCE. For a read aligned to the reverse strand (FLAG 0x10 set) the two
# are opposite, so ts must be inverted. Forward-strand reads pass through.
#
# Reads that already carry XS are left untouched, so running this twice is safe.
task minimap2TsToXs {
  input {
    File bam
    File bam_bai  # !UnusedDeclaration -- forces localization beside the bam
    String output_basename = "rna.xs_tagged"
    Int cores = 8
  }

  Int space_needed_gb = 20 + round(4*size(bam, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "quay.io/biocontainers/samtools:1.17--hd87286a_1"
    memory: "16GB"
    cpu: cores
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    set -eou pipefail
    samtools --version

    samtools view -h -@ ~{cores} ~{bam} \
      | awk 'BEGIN{OFS="\t"}
          /^@/ { print; next }
          {
            ts = ""; has_xs = 0;
            for (i = 12; i <= NF; i++) {
              if ($i ~ /^ts:A:/) { ts = substr($i, 6, 1) }
              if ($i ~ /^XS:A:/) { has_xs = 1 }
            }
            if (has_xs == 1 || (ts != "+" && ts != "-")) { print; next }
            # int($2/16)%2 tests FLAG 0x10 without needing gawk-only and()
            rev = int($2 / 16) % 2;
            if (rev == 1) { xs = (ts == "+" ? "-" : "+") } else { xs = ts }
            print $0, "XS:A:" xs
          }' \
      | samtools view -b -@ ~{cores} -o ~{output_basename}.bam -

    samtools index -@ ~{cores} ~{output_basename}.bam

    # Sanity check: regtools only sees junctions on spliced reads, so a count of
    # zero here means pvacsplice will find nothing.
    echo "alignments carrying XS:A after conversion:"
    samtools view ~{output_basename}.bam | grep -c "XS:A:" || true
  >>>

  output {
    File xs_tagged_bam = output_basename + ".bam"
    File xs_tagged_bam_bai = output_basename + ".bam.bai"
  }
}

workflow wf {
  input {
    File bam
    File bam_bai
    String? output_basename
  }
  call minimap2TsToXs {
    input:
    bam=bam,
    bam_bai=bam_bai,
    output_basename=select_first([output_basename, "rna.xs_tagged"])
  }
}
