version 1.0

import "./subworkflows/bam_readcount.wdl" as br
import "./subworkflows/vcf_readcount_annotator.wdl" as vra
import "./tools/vcf_expression_annotator.wdl" as vea
import "./tools/index_vcf.wdl" as iv
import "./tools/regtools.wdl" as reg
import "./tools/pvacsplice.wdl" as pspl


workflow pvacsplice {
  input {
    File detect_variants_vcf 
    File detect_variants_vcf_tbi 
    String sample_name = "TUMOR" 
    String normal_sample_name = "NORMAL"
    File rnaseq_bam 
    File rnaseq_bam_bai 
    File reference 
    File reference_fai 
    File reference_dict 
    File? peptide_fasta 
    Int? readcount_minimum_base_quality 
    Int? readcount_minimum_mapping_quality 
    File gene_expression_file 
    File transcript_expression_file 
    String expression_tool = "kallisto" 
    
    #REGTOOLS inputs : 
    String output_filename_tsv
    String? output_filename_vcf
    String? output_filename_bed
    String? strand
    Int? window_size
    Int? max_distance_exon 
    Int? max_distance_intron
    Boolean annotate_intronic_variant 
    Boolean annotate_exonic_variant 
    Boolean not_skipping_single_exon_transcripts 
    Boolean singecell_barcode 
    Boolean intron_motif_priority
    File reference_annotation # gtf
    
    #PVACSPLICE inputs:
    Array[String] alleles
    Array[String] prediction_algorithms
    Array[Int]? epitope_lengths_class_i
    Array[Int]? epitope_lengths_class_ii
    Int? binding_threshold
    Int? percentile_threshold
    Int? iedb_retries
    String? top_score_metric  # enum [lowest, median]
    String? additional_report_columns  # enum [sample_name]
    Int? fasta_size
    Boolean? exclude_nas
    Int? maximum_transcript_support_level  # enum [1 2 3 4 5]
    Int? normal_cov
    Int? tdna_cov
    Int? trna_cov
    Float? normal_vaf
    Float? tdna_vaf
    Float? trna_vaf
    Float? expn_val
    String? net_chop_method  # enum [cterm 20s]
    Float? net_chop_threshold
    Boolean? netmhc_stab
    Boolean? run_reference_proteome_similarity
    Int? n_threads
    Float? tumor_purity
    Boolean? allele_specific_binding_thresholds
    Int? aggregate_inclusion_binding_threshold
    Array[String]? problematic_amino_acids
    Array[String]? biotypes
    Int? aggregate_inclusion_count_limit
    Int? junction_score
    Int? variant_distance
    Boolean? save_gtf 
    Array[String]? junction_anchor_types
    Boolean? keep_tmp_files 
  }

  call br.bamReadcount as tumorRnaBamReadcount {
    input:
    vcf=detect_variants_vcf,
    vcf_tbi=detect_variants_vcf_tbi,
    sample=sample_name,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=rnaseq_bam,
    bam_bai=rnaseq_bam_bai,
    min_base_quality=readcount_minimum_base_quality,
    min_mapping_quality=readcount_minimum_mapping_quality
  }

  call vra.vcfReadcountAnnotator as addTumorRnaBamReadcountToVcf {
    input:
    vcf=tumorRnaBamReadcount.normalized_vcf,
    snv_bam_readcount_tsv=tumorRnaBamReadcount.snv_bam_readcount_tsv,
    indel_bam_readcount_tsv=tumorRnaBamReadcount.indel_bam_readcount_tsv,
    data_type="RNA",
    sample_name=sample_name
  }

  call vea.vcfExpressionAnnotator as addGeneExpressionDataToVcf {
    input:
    vcf=addTumorRnaBamReadcountToVcf.annotated_bam_readcount_vcf,
    expression_file=gene_expression_file,
    expression_tool=expression_tool,
    data_type="gene",
    sample_name=sample_name
  }

  call vea.vcfExpressionAnnotator as addTranscriptExpressionDataToVcf {
    input:
    vcf=addGeneExpressionDataToVcf.annotated_expression_vcf,
    expression_file=transcript_expression_file,
    expression_tool=expression_tool,
    data_type="transcript",
    sample_name=sample_name
  }

  call iv.indexVcf as index {
    input: vcf=addTranscriptExpressionDataToVcf.annotated_expression_vcf
  }

  call reg.regtools as runregtools {
    input:
    output_filename_tsv=output_filename_tsv,
    output_filename_vcf=output_filename_vcf,
    output_filename_bed=output_filename_bed,
    strand=strand,
    window_size=window_size,
    max_distance_exon=max_distance_exon,
    max_distance_intron=max_distance_intron,
    annotate_intronic_variant=annotate_intronic_variant,
    annotate_exonic_variant=annotate_exonic_variant,
    not_skipping_single_exon_transcripts=not_skipping_single_exon_transcripts,
    singecell_barcode=singecell_barcode,
    intron_motif_priority=intron_motif_priority,
    input_vcf=index.indexed_vcf,
    input_bam=rnaseq_bam,
    input_reference_dna_fasta=reference,
    input_reference_gtf=reference_annotation
  }

  call pspl.pvacsplice as runpvacsplice {
    input:
    n_threads=n_threads,
    input_vcf=index.indexed_vcf,
    input_vcf_tbi=index.indexed_vcf_tbi,
    input_regtools_tsv=runregtools.output_splice_junction_tsv,
    input_reference_dna_fasta=reference,
    input_reference_gtf=reference_annotation, 
    sample_name=sample_name,
    alleles=alleles,
    prediction_algorithms=prediction_algorithms,
    peptide_fasta=peptide_fasta, 
    epitope_lengths_class_i=epitope_lengths_class_i, 
    epitope_lengths_class_ii=epitope_lengths_class_ii,
    binding_threshold=binding_threshold,
    percentile_threshold=percentile_threshold,
    iedb_retries=iedb_retries,
    normal_sample_name=normal_sample_name,
    net_chop_method=net_chop_method,
    top_score_metric=top_score_metric,  
    net_chop_threshold=net_chop_threshold,
    additional_report_columns=additional_report_columns,
    fasta_size=fasta_size,
    exclude_nas=exclude_nas,
    normal_cov=normal_cov,
    tdna_cov=tdna_cov,
    trna_cov=trna_cov,
    normal_vaf=normal_vaf,
    tdna_vaf=tdna_vaf,
    trna_vaf=trna_vaf,
    expn_val=expn_val,
    maximum_transcript_support_level=maximum_transcript_support_level,
    aggregate_inclusion_binding_threshold=aggregate_inclusion_binding_threshold,
    problematic_amino_acids=problematic_amino_acids,
    biotypes=biotypes,
    aggregate_inclusion_count_limit=aggregate_inclusion_count_limit,
    junction_score=junction_score,
    variant_distance=variant_distance,
    save_gtf=save_gtf,
    junction_anchor_types=junction_anchor_types,
    allele_specific_binding_thresholds=allele_specific_binding_thresholds,
    keep_tmp_files=keep_tmp_files,
    netmhc_stab=netmhc_stab,
    run_reference_proteome_similarity=run_reference_proteome_similarity,
    tumor_purity=tumor_purity
  }


  output {
    File annotated_vcf = index.indexed_vcf
    File annotated_vcf_tbi = index.indexed_vcf_tbi
    File regtools_tsv = runregtools.output_splice_junction_tsv
    File? regtools_vcf = runregtools.output_splice_variant_vcf
    Array[File] mhc_i = runpvacsplice.mhc_i
    Array[File] mhc_ii = runpvacsplice.mhc_ii
    Array[File] combined = runpvacsplice.combined
    File? splice_transcript_combined_report = runpvacsplice.splice_transcript_combined_report
    File? splice_fasta = runpvacsplice.splice_fasta
    File? splice_fasta_fai = runpvacsplice.splice_fasta_fai
  }
}