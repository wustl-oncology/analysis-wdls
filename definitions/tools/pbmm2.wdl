version 1.0

# pbmm2 alignment of PacBio HiFi genomic reads (the DNA arm).
#
# pbmm2 is minimap2 with PacBio-tuned presets and, importantly, it preserves the
# PacBio-specific bam tags (kinetics, base modifications) that plain minimap2
# would drop. Use this for WGS/genomic HiFi; use
# tools/minimap2_rnaseq_pacbioHifi.wdl for Iso-Seq.
#
# TWO INPUT FORMATS
#
#   bams   -- unaligned HiFi bam (the native format)
#   fastqs -- .fastq.gz / .fq.gz
#
# Supply exactly one. FASTQ exists to save disk: a HiFi uBAM carrying kinetics
# is several times the size of the equivalent fastq.gz.
#
# WHAT YOU GIVE UP WITH FASTQ: the PacBio tags. Kinetics (fi/fp/ri/rp) and base
# modification tags (MM/ML, i.e. 5mC calls) do not exist in FASTQ and cannot be
# recovered. For THIS pipeline that costs nothing -- pbmm2 -> DeepSomatic ->
# Clair3 -> HiPhase -> VEP reads none of them. It would matter if you later
# wanted the methylation/DMR arm that hifisomatic.wdl has, which is built
# entirely on MM/ML.
#
# READ GROUPS: with bam input the @RG rides along. With FASTQ there is none, and
# pbmm2's docs are explicit that "--rg sets the read group" for FASTA/Q input --
# --sample alone is not sufficient. This file therefore constructs an @RG line
# for FASTQ input so the output bam is not read-group-less.

task pbmm2Align {
  input {
    File reads
    # true when `reads` is fastq/fastq.gz rather than bam. Controls read-group
    # construction and --strip, both of which behave differently.
    Boolean is_fastq = false
    # May be either the plain fasta or a pre-built .mmi index.
    File reference
    File reference_fai
    String sample_name
    String output_basename
    # Read group ID used only for FASTQ input. Must be unique per input file,
    # which is why the workflow below suffixes it with the scatter index.
    String read_group_id = "rg0"
    # BAM input only -- FASTQ has no kinetic tags to strip, and passing --strip
    # with FASTQ input is meaningless.
    Boolean strip_kinetics = false
    # -A 2 raises the match score; this is what HiFi-somatic-WDL uses and it
    # improves indel representation in low-complexity regions.
    String additional_args = "-A 2"
    Int cores = 24
  }

  String ofile_name = output_basename + ".aligned.bam"
  # fastq.gz expands considerably on the way to a sorted bam, so it needs a
  # larger multiplier than an already-binary uBAM.
  Float size_multiplier = if is_fastq then 6.0 else 2.0
  Int space_needed_gb = 20 + round(size_multiplier*size(reads, "GB") + size(reference, "GB"))
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
      ~{reads} \
      ~{ofile_name} \
      --sample ~{sample_name} \
      ~{if is_fastq then "--rg '@RG\\tID:" + read_group_id + "\\tSM:" + sample_name + "\\tPL:PACBIO\\tPU:" + read_group_id + "'" else ""} \
      --sort -j ~{cores} \
      --unmapped \
      --preset HIFI \
      --log-level INFO --log-file pbmm2.log \
      ~{if (strip_kinetics && !is_fastq) then "--strip" else ""} \
      ~{additional_args}

    echo "read groups in output:"
    samtools view -H ~{ofile_name} | grep "^@RG" || echo "(none)"
  >>>

  output {
    File aligned_bam = ofile_name
    File aligned_bam_bai = ofile_name + ".bai"
    File align_log = "pbmm2.log"
  }
}

# Fails the workflow within seconds on a bad input combination, instead of
# letting it surface as an opaque select_first error much later.
#
# Nothing depends on this task's output, so Cromwell may start pbmm2 alongside
# it -- but the workflow still fails immediately and no further calls are made.
task checkAlignInputs {
  input {
    Int n_bams
    Int n_fastqs
    Boolean skip_align
    String sample_name
  }

  runtime {
    preemptible: 1
    maxRetries: 1
    docker: "quay.io/biocontainers/samtools:1.17--hd87286a_1"
    memory: "2GB"
    cpu: 1
    disks: "local-disk 10 HDD"
  }

  command <<<
    set -eu

    if [ ~{n_bams} -eq 0 ] && [ ~{n_fastqs} -eq 0 ]; then
      echo "ERROR (~{sample_name}): no input reads. Set either bams or fastqs." >&2
      exit 1
    fi

    if [ ~{n_bams} -gt 0 ] && [ ~{n_fastqs} -gt 0 ]; then
      echo "ERROR (~{sample_name}): both bams (~{n_bams}) and fastqs (~{n_fastqs}) were supplied. Set exactly one." >&2
      exit 1
    fi

    if [ ~{n_fastqs} -gt 0 ] && [ "~{skip_align}" = "true" ]; then
      echo "ERROR (~{sample_name}): skip_align=true is only meaningful for already-aligned bams, but fastqs were supplied." >&2
      exit 1
    fi

    if [ ~{n_fastqs} -gt 0 ]; then
      echo "~{sample_name}: FASTQ mode, ~{n_fastqs} file(s). PacBio kinetic and base-modification tags are not available in this mode."
    else
      echo "~{sample_name}: BAM mode, ~{n_bams} file(s), skip_align=~{skip_align}."
    fi
  >>>

  output {
    String mode_message = read_string(stdout())
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

# Align (or just index) one or more HiFi read files and hand back a single bam.
workflow alignHifiBams {
  input {
    # Supply exactly ONE of these.
    Array[File] bams = []
    Array[File] fastqs = []
    File reference
    File reference_fai
    String sample_name
    # Only meaningful with `bams`: treat them as already aligned and just
    # merge/index. Rejected when fastqs are supplied.
    Boolean skip_align = false
    Boolean strip_kinetics = false
    # Optional per-input-file read group IDs, positionally matched to bams/fastqs.
    # Use this to keep replicates distinguishable inside the merged bam: the @RG
    # SM stays `sample_name` for every shard (they are one biological sample, and
    # every downstream caller keys off SM), while ID records which run each read
    # came from. Leave empty and IDs are auto-generated as <sample>_rg<index>.
    #
    # Only takes effect in FASTQ mode -- bam input carries its own @RG already.
    Array[String] read_group_ids = []
    String additional_pbmm2_args = "-A 2"
    Int pbmm2_cores = 24
    Int samtools_cores = 8
  }

  Boolean use_fastq = length(fastqs) > 0

  call checkAlignInputs {
    input:
    n_bams=length(bams),
    n_fastqs=length(fastqs),
    skip_align=skip_align,
    sample_name=sample_name
  }

  Array[File] input_reads = if use_fastq then fastqs else bams

  if (!skip_align) {
    # Scatter over indices rather than files so each shard gets a deterministic,
    # collision-free output name and read group id regardless of what the input
    # files happen to be called.
    scatter (idx in range(length(input_reads))) {
      call pbmm2Align {
        input:
        reads=input_reads[idx],
        is_fastq=use_fastq,
        reference=reference,
        reference_fai=reference_fai,
        sample_name=sample_name,
        output_basename=sample_name + ".shard" + idx,
        read_group_id=if length(read_group_ids) == length(input_reads) then read_group_ids[idx] else sample_name + "_rg" + idx,
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
    String input_mode = checkAlignInputs.mode_message
  }
}
