version 1.0

# pbmm2 alignment of PacBio HiFi genomic reads (the DNA arm).
#
# pbmm2 is minimap2 with PacBio-tuned presets and, importantly, it preserves the
# PacBio-specific bam tags (kinetics, base modifications) that plain minimap2
# would drop. Use this for WGS/genomic HiFi; use tools/minimap2_rnaseq_pacbioHifi.wdl for Iso-Seq.

task pbmm2Align {
  input {
    File bam_file
    # May be either the plain fasta or a pre-built .mmi index.
    File reference
    File reference_fai
    String sample_name
    Boolean strip_kinetics = false
    # -A 2 raises the match score; this is what HiFi-somatic-WDL uses and it
    # improves indel representation in low-complexity regions.
    String additional_args = "-A 2"
    Int cores = 24
  }

  String ofile_name = sub(basename(bam_file), "\\.bam$", ".aligned.bam")
  Int space_needed_gb = 20 + round(2*size(bam_file, "GB") + size(reference, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "quay.io/biocontainers/pbmm2@sha256:c1ec77296850cbdb02621bca1addcc25e510aacdabbad753ab3b0b8ba43ccd52"
    memory: "~{cores * 4}GB"
    cpu: cores
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    set -euxo pipefail
    pbmm2 --version

    pbmm2 align \
      ~{reference} \
      ~{bam_file} \
      ~{ofile_name} \
      --sample ~{sample_name} \
      --sort -j ~{cores} \
      --unmapped \
      --preset HIFI \
      --log-level INFO --log-file pbmm2.log \
      ~{if strip_kinetics then "--strip" else ""} \
      ~{additional_args}
  >>>

  output {
    File aligned_bam = ofile_name
    File aligned_bam_bai = ofile_name + ".bai"
    File align_log = "pbmm2.log"
  }
}

task samtoolsMergeIndex {
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
    set -euxo pipefail
    samtools merge -@ ~{cores} -o ~{output_basename}.bam ~{sep=' ' bams}
    samtools index -@ ~{cores} ~{output_basename}.bam
  >>>

  output {
    File merged_bam = output_basename + ".bam"
    File merged_bam_bai = output_basename + ".bam.bai"
  }
}

task samtoolsIndexHifi {
  input {
    File bam
    Int cores = 8
  }

  Int space_needed_gb = 10 + round(1.2*size(bam, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "quay.io/biocontainers/samtools:1.17--hd87286a_1"
    memory: "16GB"
    cpu: cores
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    set -euxo pipefail
    cp ~{bam} ~{basename(bam)}
    samtools index -@ ~{cores} ~{basename(bam)}
  >>>

  output {
    File indexed_bam = basename(bam)
    File indexed_bam_bai = basename(bam) + ".bai"
  }
}

# Align (or just index) one or more HiFi bams and hand back a single bam.
workflow alignHifiBams {
  input {
    Array[File] bams
    File reference
    File reference_fai
    String sample_name
    Boolean skip_align = false
    Boolean strip_kinetics = false
    String additional_pbmm2_args = "-A 2"
    Int pbmm2_cores = 24
    Int samtools_cores = 8
  }

  if (!skip_align) {
    scatter (bam in bams) {
      call pbmm2Align {
        input:
        bam_file=bam,
        reference=reference,
        reference_fai=reference_fai,
        sample_name=sample_name,
        strip_kinetics=strip_kinetics,
        additional_args=additional_pbmm2_args,
        cores=pbmm2_cores
      }
    }
  }

  Array[File] to_combine = if skip_align then bams else select_first([pbmm2Align.aligned_bam])

  if (length(to_combine) > 1) {
    call samtoolsMergeIndex {
      input:
      bams=to_combine,
      output_basename=sample_name + ".aligned",
      cores=samtools_cores
    }
  }

  if (length(to_combine) == 1) {
    call samtoolsIndexHifi {
      input: bam=to_combine[0], cores=samtools_cores
    }
  }

  output {
    File final_bam = select_first([samtoolsMergeIndex.merged_bam, samtoolsIndexHifi.indexed_bam])
    File final_bam_bai = select_first([samtoolsMergeIndex.merged_bam_bai, samtoolsIndexHifi.indexed_bam_bai])
  }
}
