version 1.0

workflow regtools_workflow {
  input {
    File input_vcf
    File input_bam
    File input_bam_bai
    File input_reference_dna_fasta
    File input_reference_gtf

    String output_filename_tsv = "splice_junction.tsv"
    String? output_filename_vcf = "splice_variant.vcf"
    String? output_filename_bed = "splice_variant.bed"

    String strand = "unstranded" # [first, second, unstranded]
    Int? window_size
    Int? max_distance_exon # max distance from exon/intron boundary to annotate a variant in exonic region as splicing variant
    Int? max_distance_intron

    Boolean annotate_intronic_variant = false
    Boolean annotate_exonic_variant = false
    Boolean not_skipping_single_exon_transcripts = false
    Boolean singecell_barcode = false
    Boolean intron_motif_priority = false
  }

  Float input_size = size([input_vcf, input_bam, input_reference_dna_fasta, input_reference_gtf, input_bam_bai], "GB")
  Int space_needed_gb = 10 + ceil(input_size)

  call regtools {
    input:
      input_vcf = input_vcf,
      input_bam = input_bam,
      input_bam_bai = input_bam_bai,
      input_reference_dna_fasta = input_reference_dna_fasta,
      input_reference_gtf = input_reference_gtf,
      output_filename_tsv = output_filename_tsv,
      output_filename_vcf = output_filename_vcf,
      output_filename_bed = output_filename_bed,
      strand = strand,
      window_size = window_size,
      max_distance_exon = max_distance_exon,
      max_distance_intron = max_distance_intron,
      annotate_intronic_variant = annotate_intronic_variant,
      annotate_exonic_variant = annotate_exonic_variant,
      not_skipping_single_exon_transcripts = not_skipping_single_exon_transcripts,
      singecell_barcode = singecell_barcode,
      intron_motif_priority = intron_motif_priority,
      space_needed_gb = space_needed_gb
  }

  output {
    File output_splice_junction_tsv = regtools.output_splice_junction_tsv
    File? output_splice_variant_vcf = regtools.output_splice_variant_vcf
  }
}

task regtools {
  input {
    String output_filename_tsv
    String? output_filename_vcf
    String? output_filename_bed

    String strand
    Int? window_size
    Int? max_distance_exon
    Int? max_distance_intron

    Boolean annotate_intronic_variant
    Boolean annotate_exonic_variant
    Boolean not_skipping_single_exon_transcripts
    Boolean singecell_barcode
    Boolean intron_motif_priority

    File input_vcf
    File input_bam # indexed,aligned (and preferably sorted) bam or cram
    File input_bam_bai
    File input_reference_dna_fasta
    File input_reference_gtf

    Int space_needed_gb
  }

  Map[String, String] strandness = {
    "first": "RF",
    "second": "FR",
    "unstranded": "XS"
  }

  command <<<
    /regtools/build/regtools cis-splice-effects identify \
    -o ~{output_filename_tsv} \
    ~{if defined(output_filename_vcf) then "-v ~{output_filename_vcf}" else ""} \
    ~{if defined(output_filename_bed) then "-j ~{output_filename_bed}" else ""} \
    -s ~{strandness[strand]} \
    ~{if defined(window_size) then "-w ~{window_size}" else ""} \
    ~{if defined(max_distance_exon) then "-e ~{max_distance_exon}" else ""} \
    ~{if defined(max_distance_intron) then "-i ~{max_distance_intron}" else ""} \
    ~{if annotate_intronic_variant then "-I" else ""} \
    ~{if annotate_exonic_variant then "-E" else ""} \
    ~{if not_skipping_single_exon_transcripts then "-S" else ""} \
    ~{if singecell_barcode then "-b" else ""} \
    ~{if intron_motif_priority then "-C" else ""} \
    ~{input_vcf} ~{input_bam} ~{input_reference_dna_fasta} ~{input_reference_gtf}
  >>>

  runtime {
    memory: "44GB"
    preemptible: 1
    maxRetries: 2
    docker: "griffithlab/regtools:release-1.0.0"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  output {
    File output_splice_junction_tsv = output_filename_tsv
    File? output_splice_variant_vcf = output_filename_vcf
  }
}
