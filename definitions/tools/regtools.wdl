version 1.0

task regtools {
  input {
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

    File input_vcf
    File input_bam # indexed,aligned (and preferably sorted) bam or cram
    File input_reference_dna_fasta
    File input_reference_gtf

  }

  Float input_size = size([input_vcf,input_bam, input_reference_dna_fasta,input_reference_gtf], "GB")  
  Int space_needed_gb = 10 + round(input_size) 
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "44GB"
    docker: "griffithlab/regtools:release-1.0.0"
    disks: "local-disk ~{space_needed_gb} HDD"
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

  output {
    File output_splice_junction_tsv = output_filename_tsv
    File? output_splice_variant_vcf = output_filename_vcf
  }
}

