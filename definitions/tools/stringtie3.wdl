version 1.0

# StringTie 3 quantification of long-read (Iso-Seq / ONT) alignments.
# compare to StringTie2, StringTie3 is better at handling pacbio hifi isoseq data
# StringTie3 features an optimized module designed to better handle and discard incomplete transcripts caused by internal poly(A)-priming artifacts (common in pacbio hifi).

# this wld uses reference guided and transcript estimation mode. (to be used for vcf annotation for pvacseq input)
# might want to add more options (short read , denovo modes) in this tool wdl in the future to support other use cases.

task stringtie3 {
  input {
    File bam
    File bam_bai  # !UnusedDeclaration -- present to force localization beside the bam
    File reference_annotation
    String sample_name
    # Iso-Seq FLNC is strand resolved, so the alignments are already oriented.
    # "unstranded" is the right default when minimap2 was run with -uf.
    String strand = "unstranded"  # [first, second, unstranded]
    Boolean reference_only = true
    Int cores = 12
    String docker_image = "quay.io/biocontainers/stringtie:3.0.3--h29c0135_0"
  }

  Int space_needed_gb = 10 + round(3*size(bam, "GB") + size(reference_annotation, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: docker_image
    memory: "32GB"
    cpu: cores
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String transcripts = "stringtie_transcripts.gtf"
  String expression = "stringtie_gene_expression.tsv"
  Map[String, String] strandness = {
    "first": "--rf", "second": "--fr", "unstranded": ""
  }

  command <<<
    set -eou pipefail
    stringtie --version

    stringtie \
      -L \
      ~{if reference_only then "-e" else ""} \
      -p ~{cores} \
      ~{strandness[strand]} \
      -G ~{reference_annotation} \
      -o ~{transcripts} \
      -A ~{expression} \
      -l ~{sample_name} \
      ~{bam}
  >>>

  output {
    File transcript_gtf = transcripts
    File gene_expression_tsv = expression
  }
}

workflow wf {
  input {
    File bam
    File bam_bai
    File reference_annotation
    String sample_name
    String? strand
    Boolean? reference_only
  }

  call stringtie3 {
    input:
    bam=bam,
    bam_bai=bam_bai,
    reference_annotation=reference_annotation,
    sample_name=sample_name,
    strand=select_first([strand, "unstranded"]),
    reference_only=select_first([reference_only, true])
  }
}
