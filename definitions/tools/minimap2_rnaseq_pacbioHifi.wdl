version 1.0

# Spliced alignment of PacBio Iso-Seq FLNC reads with minimap2.
#
# FLNC ("full-length non-concatemer") reads come out of `isoseq refine` as an
# UNALIGNED bam. They are already full length, primer trimmed, polyA trimmed and
# strand resolved -- that last property is why `-uf` (force the transcript strand
# to match the read orientation) is correct here. Do NOT copy this preset for raw
# cDNA or for genomic HiFi reads.
#
# Split into three tasks so each one can use a stock, single-tool container.
# If you have an image with both minimap2 and samtools you can collapse
# bamToFastqLr + minimap2Lr + samtoolsSortIndexLr into one piped command and save
# a lot of intermediate disk.

task bamToFastqLr {
  input {
    File flnc_bam
    Int cores = 4
  }

  Int space_needed_gb = 10 + round(4*size(flnc_bam, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "quay.io/biocontainers/samtools:1.17--hd87286a_1"
    memory: "8GB"
    cpu: cores
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String outfile = "flnc.fastq"
  command <<<
    set -eou pipefail
    samtools --version
    # -0 /dev/null -s /dev/null keeps samtools from silently dropping reads it
    # thinks are paired; FLNC is single-end so everything should land in stdout.
    samtools fastq -@ ~{cores} -0 /dev/null -s /dev/null ~{flnc_bam} > ~{outfile}
  >>>

  output {
    File fastq = outfile
  }
}

task minimap2Lr {
  input {
    File fastq
    File reference
    File reference_fai
    # `splice:hq` is the minimap2 preset for PacBio Iso-Seq / HiFi cDNA.
    String preset = "splice:hq"
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

    # -uf          : transcript strand comes from read orientation (valid for FLNC)
    # --secondary=no: StringTie should not see multimapping records
    # --MD          : emit MD tags so bam-readcount can work off this bam later
    # -Y            : soft clip supplementary alignments, matching the short-read arm
    minimap2 \
      -t ~{cores} \
      -ax ~{preset} \
      -uf \
      --secondary=no \
      --MD \
      -Y \
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

workflow minimap2Align {
  input {
    File flnc_bam
    File reference
    File reference_fai
    String output_basename = "flnc.aligned"
    String preset = "splice:hq"
    Int minimap2_cores = 16
    Int samtools_cores = 8
    String? minimap2_docker_image
  }

  call bamToFastqLr {
    input: flnc_bam=flnc_bam
  }

  call minimap2Lr {
    input:
    fastq=bamToFastqLr.fastq,
    reference=reference,
    reference_fai=reference_fai,
    preset=preset,
    cores=minimap2_cores,
    docker_image=select_first([minimap2_docker_image, "quay.io/biocontainers/minimap2:2.28--he4a0461_0"])
  }

  call samtoolsSortIndexLr {
    input:
    aligned_sam=minimap2Lr.aligned_sam,
    output_basename=output_basename,
    cores=samtools_cores
  }

  output {
    File aligned_bam = samtoolsSortIndexLr.sorted_bam
    File aligned_bam_bai = samtoolsSortIndexLr.sorted_bam_bai
  }
}
