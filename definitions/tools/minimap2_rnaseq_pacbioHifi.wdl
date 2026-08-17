version 1.0

import "bam2fastq_pacbio.wdl" as b2f

# Spliced alignment of PacBio Iso-Seq FLNC reads with minimap2.
#
# FLNC ("full-length non-concatemer") reads come out of `isoseq refine` as an
# UNALIGNED bam. They are already full length, primer trimmed, polyA trimmed and
# strand resolved -- that last property is why `-uf` (force the transcript strand
# to match the read orientation) is correct here. Do NOT copy this preset for raw
# cDNA or for genomic HiFi reads.
#
# bam -> fastq is done by bam2fastq_pacbio.wdl's bam2fastqPacbio (pbtk), not
# samtools fastq: samtools fastq routes any read carrying the
# paired-in-sequencing flag to its singleton output, which this pipeline sent
# to /dev/null on the assumption FLNC reads are never paired. When that
# assumption doesn't hold, every read is silently discarded -- samtools
# reports "processed N reads" while writing a 0-byte fastq, exits 0, and the
# failure only surfaces downstream as an empty alignment. bam2fastq has no
# such flag-based routing. See bam2fastq_pacbio.wdl for the full writeup.

task minimap2Lr {
  input {
    File fastq
    # Either the genome fasta or a prebuilt .mmi. A prebuilt index skips the
    # several minutes and ~10GB of RAM minimap2 otherwise spends indexing the
    # genome on every run.
    #
    # CAVEAT: an .mmi stores the -k/-w it was built with, and those OVERRIDE the
    # preset's values. If the index was built with plain `minimap2 -d` rather
    # than with `-x splice:hq`, minimap2 prints a compatibility warning and uses
    # the index's parameters, which is not the same alignment the preset would
    # otherwise give. Check stderr for that warning the first time.
    File reference
    File reference_fai
    # FLNC -> fastq loses the @RG that `isoseq refine` wrote, so it is rebuilt
    # here. Distinct IDs per replicate are what let `samtools merge` combine the
    # shards cleanly and what keeps replicates separable afterwards:
    #   samtools view -r <id> merged.bam
    # SM is identical across replicates -- they are one biological sample.
    String sample_name
    String read_group_id
    # `splice:hq` is the minimap2 preset for PacBio Iso-Seq / HiFi cDNA.
    String preset = "splice:hq"
    # Extra minimap2 flags. Empty by default so the command matches a
    # configuration proven against ctat-LR-fusion in practice. See the note on
    # --secondary=no below before adding anything here.
    String additional_args = ""
    Int cores = 16
    # NOTE: verify this build hash against
    # https://quay.io/repository/biocontainers/minimap2?tab=tags
    String docker_image = "quay.io/biocontainers/minimap2:2.28--he4a0461_0"
  }

  Int space_needed_gb = 20 + round(6*size(fastq, "GB") + size([reference, reference_fai], "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: docker_image
    memory: "~{cores * 4}GB"
    cpu: cores
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String outfile = "aligned.sam"
  command <<<
    set -eou pipefail
    minimap2 --version

    # This single alignment feeds THREE consumers -- StringTie3, regtools
    # (pvacsplice) and ctat-LR-fusion (pvacfuse) -- so the flags are chosen to
    # satisfy all three rather than optimised for any one.
    #
    # -uf   : transcript strand from read orientation. Valid for FLNC, which is
    #         strand resolved; would be wrong for unoriented cDNA.
    # --MD  : MD tags, so bam-readcount can work off this bam later.
    # -R    : read group. Pure metadata -- it adds an @RG header line and an RG:Z
    #         tag per record, and does not change alignment in any way, so this
    #         is not a meaningful deviation from the proven ctat config.
    #
    # NOT set, deliberately:
    #
    #   --secondary=no  Tempting, to keep multimappers away from StringTie. But
    #                   StringTie reads minimap2's ts tag natively (the manual:
    #                   "StringTie can recognize the ts tag as well, if the XS
    #                   tag is missing"), and the flag is NOT part of the
    #                   configuration that has been shown to work with
    #                   ctat-LR-fusion here. Fusion evidence lives in
    #                   supplementary alignments (0x800), which --secondary=no
    #                   does not touch -- but there is no reason to deviate from
    #                   a proven fusion config on a speculative quantification
    #                   benefit. Pass it via additional_args if StringTie
    #                   abundances ever look inflated.
    #
    #   -Y              Soft-clips supplementary alignments. Also absent from the
    #                   proven ctat config; hard clipping is minimap2's default
    #                   and ctat has been run against it successfully.
    minimap2 \
      -t ~{cores} \
      -ax ~{preset} \
      -uf \
      --MD \
      -R '@RG\tID:~{read_group_id}\tSM:~{sample_name}\tPL:PACBIO\tPU:~{read_group_id}' \
      ~{additional_args} \
      ~{reference} \
      ~{fastq} \
      > ~{outfile}
  >>>

  output {
    File aligned_sam = outfile
  }
}

task samtoolsSortIndexLr {
  input {
    File aligned_sam
    String output_basename = "flnc.aligned"
    Int cores = 8
  }

  Int space_needed_gb = 20 + round(3*size(aligned_sam, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "quay.io/biocontainers/samtools:1.17--hd87286a_1"
    memory: "32GB"
    cpu: cores
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    set -eou pipefail
    samtools sort -@ ~{cores} -m 2G -o ~{output_basename}.bam ~{aligned_sam}
    samtools index -@ ~{cores} ~{output_basename}.bam
  >>>

  output {
    File sorted_bam = output_basename + ".bam"
    File sorted_bam_bai = output_basename + ".bam.bai"
  }
}

task samtoolsMergeIndexLr {
  input {
    Array[File] bams
    String output_basename
    Int cores = 8
  }

  Int space_needed_gb = 20 + round(2.5*size(bams, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "quay.io/biocontainers/samtools:1.17--hd87286a_1"
    memory: "~{cores * 4}GB"
    cpu: cores
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    set -eou pipefail
    # -c/-p merge identical @RG/@PG header lines rather than erroring. The shards
    # have distinct RG IDs by construction, so nothing is actually collapsed.
    samtools merge -@ ~{cores} -c -p -o ~{output_basename}.bam ~{sep=' ' bams}
    samtools index -@ ~{cores} ~{output_basename}.bam

    echo "read groups in merged bam:"
    samtools view -H ~{output_basename}.bam | grep "^@RG" || echo "(none)"
  >>>

  output {
    File merged_bam = output_basename + ".bam"
    File merged_bam_bai = output_basename + ".bam.bai"
  }
}

workflow minimap2Align {
  input {
    # One or more unaligned FLNC bams from `isoseq refine`. Multiple entries are
    # treated as replicates of ONE biological sample: each is aligned
    # independently (so it can carry its own @RG) and the results are merged
    # into a single sorted, indexed bam.
    Array[File] flnc_bams
    File reference
    File reference_fai
    String sample_name
    # Optional per-replicate @RG IDs, positionally matched to flnc_bams.
    # Leave empty and they are auto-generated as <sample>_rna_rg<index>.
    Array[String] read_group_ids = []
    String output_basename = "flnc.aligned"
    String preset = "splice:hq"
    String additional_args = ""
    Int minimap2_cores = 16
    Int samtools_cores = 8
    String? minimap2_docker_image
  }

  # Align each replicate on its own so it can carry a distinct @RG. This also
  # parallelises the alignment across replicates rather than serialising one
  # large concatenated fastq.
  scatter (idx in range(length(flnc_bams))) {
    call b2f.bam2fastqPacbio as bamToFastqLr {
      input: flnc_bam=flnc_bams[idx]
    }

    call minimap2Lr {
      input:
      fastq=bamToFastqLr.fastq,
      reference=reference,
      reference_fai=reference_fai,
      sample_name=sample_name,
      read_group_id=if length(read_group_ids) == length(flnc_bams) then read_group_ids[idx] else sample_name + "_rna_rg" + idx,
      preset=preset,
      additional_args=additional_args,
      cores=minimap2_cores,
      docker_image=select_first([minimap2_docker_image, "quay.io/biocontainers/minimap2:2.28--he4a0461_0"])
    }

    call samtoolsSortIndexLr {
      input:
      aligned_sam=minimap2Lr.aligned_sam,
      output_basename=output_basename + ".shard" + idx,
      cores=samtools_cores
    }
  }

  if (length(flnc_bams) > 1) {
    call samtoolsMergeIndexLr {
      input:
      bams=samtoolsSortIndexLr.sorted_bam,
      output_basename=output_basename,
      cores=samtools_cores
    }
  }

  output {
    # Single replicate -> the one sorted shard. Multiple -> the merged bam.
    File aligned_bam = select_first([samtoolsMergeIndexLr.merged_bam, samtoolsSortIndexLr.sorted_bam[0]])
    File aligned_bam_bai = select_first([samtoolsMergeIndexLr.merged_bam_bai, samtoolsSortIndexLr.sorted_bam_bai[0]])
  }
}
