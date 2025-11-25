version 1.0

import "./tools/split_n_cigar_reads.wdl" as sncr
import "./tools/left_align_indels.wdl" as lai
import "./subworkflows/bam_readcount.wdl" as br
import "./subworkflows/vcf_readcount_annotator.wdl" as vra
import "./tools/vcf_expression_annotator.wdl" as vea
import "./tools/index_vcf.wdl" as iv
import "./tools/pvacseq.wdl" as p
import "./tools/pvacseq_aggregated_report_to_preferred_transcripts_list.wdl" as ptl
import "./tools/variants_to_table.wdl" as vtt
import "./tools/add_vep_fields_to_table.wdl" as avftt

workflow pvacseq {
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
    File? genes_of_interest_file
    Int? readcount_minimum_base_quality
    Int? readcount_minimum_mapping_quality
    File gene_expression_file
    File transcript_expression_file
    String expression_tool = "kallisto"
    Array[String] alleles
    Array[String] prediction_algorithms
    Array[Int]? epitope_lengths_class_i
    Array[Int]? epitope_lengths_class_ii
    Int? binding_threshold
    Int? percentile_threshold
    String? percentile_threshold_strategy
    Float? minimum_fold_change
    String? top_score_metric  # enum [lowest, median]
    String? top_score_metric2  # enum [ic50, percentile]
    String? additional_report_columns  # enum [sample_name]
    Int? fasta_size
    Int? downstream_sequence_length
    Boolean? exclude_nas
    File? phased_proximal_variants_vcf
    File? phased_proximal_variants_vcf_tbi
    Array[String]? transcript_prioritization_strategy # allowed values ['canonical', 'mane_select', 'tsl']
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
    Int? iedb_retries
    Array[String] variants_to_table_fields = ["CHROM", "POS", "ID", "REF", "ALT"]
    Array[String] variants_to_table_genotype_fields = ["GT", "AD", "AF", "DP", "RAD", "RAF", "RDP", "GX", "TX"]
    Array[String] vep_to_table_fields = ["HGVSc", "HGVSp"]
    Float? tumor_purity
    Boolean? allele_specific_binding_thresholds
    Int? aggregate_inclusion_binding_threshold
    Int? aggregate_inclusion_count_limit
    Array[String]? problematic_amino_acids
    Boolean? allele_specific_anchors
    Float? anchor_contribution_threshold
    String? prefix = "pvacseq"
    Array[String]? biotypes
    Boolean? allow_incomplete_transcripts
    String? netmhciipan_version # enum [4.3, 4.2, 4.1, 4.0]
  }

  call sncr.splitNCigarReads as tumorRnaSplitNCigarReads{
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=rnaseq_bam,
    bam_bai=rnaseq_bam_bai,
    vcf=detect_variants_vcf,
    vcf_tbi=detect_variants_vcf_tbi,
    output_bam_basename="split_n_cigar"
  }

  call lai.leftAlignIndels as tumorRnaLeftAlignIndels{
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=tumorRnaSplitNCigarReads.split_n_cigar_bam,
    bam_bai=tumorRnaSplitNCigarReads.split_n_cigar_bam_bai,
    output_bam_basename="left_align_indels"
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
    min_mapping_quality=readcount_minimum_mapping_quality,
    indel_counting_bam=tumorRnaLeftAlignIndels.left_align_indels_bam,
    indel_counting_bai=tumorRnaLeftAlignIndels.left_align_indels_bam_bai
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

  call p.pvacseq as ps {
    input:
    input_vcf=index.indexed_vcf,
    input_vcf_tbi=index.indexed_vcf_tbi,
    sample_name=sample_name,
    alleles=alleles,
    prediction_algorithms=prediction_algorithms,
    epitope_lengths_class_i=epitope_lengths_class_i,
    epitope_lengths_class_ii=epitope_lengths_class_ii,
    binding_threshold=binding_threshold,
    percentile_threshold=percentile_threshold,
    percentile_threshold_strategy=percentile_threshold_strategy,
    normal_sample_name=normal_sample_name,
    minimum_fold_change=minimum_fold_change,
    top_score_metric=top_score_metric,
    top_score_metric2=top_score_metric2,
    additional_report_columns=additional_report_columns,
    fasta_size=fasta_size,
    downstream_sequence_length=downstream_sequence_length,
    exclude_nas=exclude_nas,
    phased_proximal_variants_vcf=phased_proximal_variants_vcf,
    phased_proximal_variants_vcf_tbi=phased_proximal_variants_vcf_tbi,
    transcript_prioritization_strategy=transcript_prioritization_strategy,
    maximum_transcript_support_level=maximum_transcript_support_level,
    normal_cov=normal_cov,
    tdna_cov=tdna_cov,
    trna_cov=trna_cov,
    normal_vaf=normal_vaf,
    tdna_vaf=tdna_vaf,
    trna_vaf=trna_vaf,
    expn_val=expn_val,
    net_chop_method=net_chop_method,
    net_chop_threshold=net_chop_threshold,
    netmhc_stab=netmhc_stab,
    run_reference_proteome_similarity=run_reference_proteome_similarity,
    peptide_fasta=peptide_fasta,
    genes_of_interest_file=genes_of_interest_file,
    n_threads=n_threads,
    iedb_retries=iedb_retries,
    netmhciipan_version=netmhciipan_version,
    tumor_purity=tumor_purity,
    allele_specific_binding_thresholds=allele_specific_binding_thresholds,
    aggregate_inclusion_binding_threshold=aggregate_inclusion_binding_threshold,
    aggregate_inclusion_count_limit=aggregate_inclusion_count_limit,
    problematic_amino_acids=problematic_amino_acids,
    allele_specific_anchors=allele_specific_anchors,
    anchor_contribution_threshold=anchor_contribution_threshold,
    biotypes=biotypes,
    allow_incomplete_transcripts=allow_incomplete_transcripts
  }

  call vtt.variantsToTable {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    vcf=index.indexed_vcf,
    vcf_tbi=index.indexed_vcf_tbi,
    fields=variants_to_table_fields,
    genotype_fields=variants_to_table_genotype_fields
  }

  if (length(select_all([ps.mhc_i_aggregated_report, ps.mhc_ii_aggregated_report])) > 0 ) {
      call ptl.pvacseqAggregatedReportToPreferredTranscriptsList as pvacseqAggregatedReportToPreferredTranscriptsList {
        input:
        pvacseq_aggregated_report=select_first([ps.mhc_i_aggregated_report, ps.mhc_ii_aggregated_report])
      }

      call avftt.addVepFieldsToTable as addVepFieldsToTableWithPreferredTranscriptsTsv{
        input:
        vcf=index.indexed_vcf,
        vep_fields=vep_to_table_fields,
        tsv=variantsToTable.variants_tsv,
        prefix=prefix,
        preferred_transcripts_tsv=pvacseqAggregatedReportToPreferredTranscriptsList.preferred_transcripts_tsv
      }
  }
  if (length(select_all([ps.mhc_i_aggregated_report, ps.mhc_ii_aggregated_report])) == 0 ) {
      call avftt.addVepFieldsToTable as addVepFieldsToTableWithoutPreferredTranscriptsTsv{
        input:
        vcf=index.indexed_vcf,
        vep_fields=vep_to_table_fields,
        tsv=variantsToTable.variants_tsv,
        prefix=prefix,
      }
  }

  output {
    File annotated_vcf = index.indexed_vcf
    File annotated_vcf_tbi = index.indexed_vcf_tbi
    File annotated_tsv = select_first([addVepFieldsToTableWithPreferredTranscriptsTsv.annotated_variants_tsv, addVepFieldsToTableWithoutPreferredTranscriptsTsv.annotated_variants_tsv])
    Array[File] mhc_i = ps.mhc_i
    File? mhc_i_log = ps.mhc_i_log
    Array[File] mhc_ii = ps.mhc_ii
    File? mhc_ii_log = ps.mhc_ii_log
    Array[File] combined = ps.combined
    File indel_counting_bam = tumorRnaLeftAlignIndels.left_align_indels_bam
    File indel_counting_bai = tumorRnaLeftAlignIndels.left_align_indels_bam_bai
  }
}
