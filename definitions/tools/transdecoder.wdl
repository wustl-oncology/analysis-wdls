version 1.0

task transdecoder {
  input {
    File stringtie_gtf
    File reference
    File reference_fai
    String output_prefix
    Int minimum_orf_length = 30
    String docker_image = "quay.io/biocontainers/transdecoder:6.0.0--pl5321hdfd78af_0"
  }

  String gtf_base = basename(stringtie_gtf, ".gtf")

  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "32GB"
    cpu: 4
    docker: docker_image
    disks: "local-disk 100 HDD"
  }

  command <<<
    set -e
    export PATH="/usr/local/bin:/usr/local/opt/transdecoder:/opt/conda/bin:/usr/bin:/bin"
    cp -L ~{reference} genome.fa
    cp -L ~{reference_fai} genome.fa.fai
    TransDecoder \
      --genome genome.fa \
      --gtf ~{stringtie_gtf} \
      -m ~{minimum_orf_length} \
      -S \
      -O transdecoder_output

    cp "transdecoder_output/~{gtf_base}.cDNA.fasta" "~{output_prefix}.cDNA.fasta"
    cp "transdecoder_output/~{gtf_base}.cDNA.fasta.transdecoder.pep" "~{output_prefix}.pep"
    cp "transdecoder_output/~{gtf_base}.cDNA.fasta.transdecoder.cds" "~{output_prefix}.cds"
    cp "transdecoder_output/~{gtf_base}.cDNA.fasta.transdecoder.gff3" "~{output_prefix}.transcript.gff3"
    cp "transdecoder_output/~{gtf_base}.cDNA.fasta.transdecoder.genome.gff3" "~{output_prefix}.genome.gff3"

    printf 'metric\tvalue\n' > "~{output_prefix}.summary.tsv"
    printf 'minimum_orf_length_aa\t~{minimum_orf_length}\n' >> "~{output_prefix}.summary.tsv"
    printf 'input_transcripts\t%s\n' "$(grep -c '^>' "~{output_prefix}.cDNA.fasta")" >> "~{output_prefix}.summary.tsv"
    printf 'predicted_orf_proteins\t%s\n' "$(grep -c '^>' "~{output_prefix}.pep")" >> "~{output_prefix}.summary.tsv"
    printf 'predicted_orf_cds_sequences\t%s\n' "$(grep -c '^>' "~{output_prefix}.cds")" >> "~{output_prefix}.summary.tsv"
  >>>

  output {
    File cdna_fasta = "~{output_prefix}.cDNA.fasta"
    File peptide_fasta = "~{output_prefix}.pep"
    File cds_fasta = "~{output_prefix}.cds"
    File transcript_gff3 = "~{output_prefix}.transcript.gff3"
    File genome_gff3 = "~{output_prefix}.genome.gff3"
    File summary = "~{output_prefix}.summary.tsv"
  }
}
