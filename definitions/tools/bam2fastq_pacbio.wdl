version 1.0

# Converts an unaligned, demultiplexed PacBio bam (e.g. FLNC reads from `isoseq refine`) to fastq


task bam2fastqPacbio {
  input {
    File flnc_bam
    String output_basename = "flnc"
    Int cores = 4
  }


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
  # guaranteed to be writable.

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
