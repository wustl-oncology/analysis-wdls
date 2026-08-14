version 1.0

# Chop the reference into equal-size chunks for scattering variant calling.
#
# Unlike the short-read exome arm, there is no bait/target interval list to
# scatter over here, so we tile the genome instead. 75 Mbp gives ~42 chunks for
# hg38, which is a reasonable balance between parallelism and per-task overhead.
# alt/random/unplaced contigs, chrM and chrEBV are dropped -- they are noise for
# somatic calling and chrM in particular will blow up a somatic caller's runtime.
task splitContigs {
  input {
    File reference_fai
    Int chunk_size = 75000000
    Boolean exclude_nonprimary = true
  }

  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "quay.io/biocontainers/bedtools:2.31.0--hf5e1c6e_2"
    memory: "4GB"
    cpu: 2
    disks: "local-disk 10 HDD"
  }

  command <<<
    set -euxo pipefail
    bedtools --version

    bedtools makewindows -g ~{reference_fai} -w ~{chunk_size} > contigs.bed

    if [[ "~{exclude_nonprimary}" == "true" ]]; then
      grep -v -E "random|chrUn|chrM|chrEBV|_alt|HLA-" contigs.bed > primary.bed
    else
      cp contigs.bed primary.bed
    fi

    # one bed file per chunk
    split -l 1 primary.bed contigs_split.
    for file in contigs_split.*; do mv -- "$file" "$file.bed"; done

    echo "chunks created:"
    find . -maxdepth 1 -name 'contigs_split.*.bed' | wc -l
  >>>

  output {
    Array[File] contigs = glob("contigs_split.*.bed")
    File all_contigs_bed = "primary.bed"
  }
}

workflow wf {
  input {
    File reference_fai
    Int? chunk_size
  }
  call splitContigs {
    input:
    reference_fai=reference_fai,
    chunk_size=select_first([chunk_size, 75000000])
  }
}
