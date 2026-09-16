version 1.0

import "./tools/orfanage.wdl" as orfanage
import "./tools/pvacnc_prepare.wdl" as prep
import "./tools/pvacnc_transdecoder_prepare.wdl" as tdprep
import "./tools/pvacnc_annotate_pvacview.wdl" as view
import "./tools/pvacseq.wdl" as pvacseq
import "./tools/stringtie.wdl" as stringtie
import "./tools/transdecoder.wdl" as transdecoder

workflow pvacnc {
  input {
    String sample_name
    String tumor_sample_name
    String normal_sample_name
    File rnaseq_bam
    String? strand
    File with_e_stringtie_gtf
    File source_vcf
    File reference
    File reference_fai
    File reference_annotation
    Array[File] canonical_pvacseq_outputs
    Array[String] alleles
    Array[String] prediction_algorithms
    File? peptide_fasta
    File? genes_of_interest_file
    File? phased_proximal_variants_vcf
    File? phased_proximal_variants_vcf_tbi

    Array[Int]? epitope_lengths_class_i
    Array[Int]? epitope_lengths_class_ii
    Int? binding_threshold
    Float? binding_percentile_threshold
    Float? presentation_percentile_threshold
    Float? immunogenicity_percentile_threshold
    String? percentile_threshold_strategy
    Int? iedb_retries
    String? net_chop_method
    String? top_score_metric
    Array[String]? top_score_metric2
    Float? net_chop_threshold
    String? additional_report_columns
    Int? fasta_size
    Int? downstream_sequence_length
    Boolean exclude_nas = false
    Float? minimum_fold_change
    Int? normal_cov
    Int? tdna_cov
    Int? trna_cov
    Float? normal_vaf
    Float? tdna_vaf
    Float? trna_vaf
    Float? expn_val
    Int? maximum_transcript_support_level
    Array[String]? transcript_prioritization_strategy
    Int? aggregate_inclusion_binding_threshold
    Int? aggregate_inclusion_count_limit
    Array[String]? problematic_amino_acids
    Float? anchor_contribution_threshold
    Array[String]? biotypes
    String? netmhciipan_version
    File? reference_scores_zip
    Boolean allele_specific_binding_thresholds = false
    Boolean netmhc_stab = false
    Boolean run_reference_proteome_similarity = false
    Boolean allele_specific_anchors = false
    Boolean allow_incomplete_transcripts = false
    Boolean use_normalized_percentiles = false
    Float? tumor_purity
    Boolean run_ml_predictions = false
    Float? ml_threshold_accept
    Float? ml_threshold_reject
    Int? n_threads

    String orfanage_docker = "quay.io/biocontainers/orfanage:1.2.0--heaafb18_2"
    String transdecoder_docker = "quay.io/biocontainers/transdecoder:6.0.0--pl5321hdfd78af_0"
    Int transdecoder_minimum_orf_length = 30
    String preparation_docker = "jinglunli/pvacnc:0.2.2"
  }

  call stringtie.stringtie as stringtieNoE {
    input:
      strand = strand,
      reference_annotation = reference_annotation,
      sample_name = sample_name + ".no_e",
      bam = rnaseq_bam,
      reference_estimation_only = false
  }

  call orfanage.orfanageBest as withEOrfanage {
    input:
      stringtie_gtf = with_e_stringtie_gtf,
      reference_annotation = reference_annotation,
      reference = reference,
      output_prefix = sample_name + ".with_e",
      docker_image = orfanage_docker
  }

  call orfanage.orfanageBest as noEOrfanage {
    input:
      stringtie_gtf = stringtieNoE.transcript_gtf,
      reference_annotation = reference_annotation,
      reference = reference,
      output_prefix = sample_name + ".no_e",
      docker_image = orfanage_docker
  }

  call transdecoder.transdecoder as withETransdecoder {
    input:
      stringtie_gtf = with_e_stringtie_gtf,
      reference = reference,
      reference_fai = reference_fai,
      output_prefix = sample_name + ".with_e.transdecoder",
      minimum_orf_length = transdecoder_minimum_orf_length,
      docker_image = transdecoder_docker
  }

  call transdecoder.transdecoder as noETransdecoder {
    input:
      stringtie_gtf = stringtieNoE.transcript_gtf,
      reference = reference,
      reference_fai = reference_fai,
      output_prefix = sample_name + ".no_e.transdecoder",
      minimum_orf_length = transdecoder_minimum_orf_length,
      docker_image = transdecoder_docker
  }

  call prep.pvacncPrepare {
    input:
      sample_name = sample_name,
      tumor_sample_name = tumor_sample_name,
      source_vcf = source_vcf,
      reference = reference,
      reference_fai = reference_fai,
      reference_annotation = reference_annotation,
      with_e_orfanage_gtf = withEOrfanage.best_gtf,
      no_e_orfanage_gtf = noEOrfanage.best_gtf,
      canonical_pvacseq_outputs = canonical_pvacseq_outputs,
      docker_image = preparation_docker
  }

  call tdprep.pvacncTransdecoderPrepare {
    input:
      sample_name = sample_name,
      tumor_sample_name = tumor_sample_name,
      source_vcf = source_vcf,
      reference = reference,
      reference_fai = reference_fai,
      reference_annotation = reference_annotation,
      with_e_stringtie_gtf = with_e_stringtie_gtf,
      with_e_transdecoder_genome_gff3 = withETransdecoder.genome_gff3,
      with_e_transdecoder_transcript_gff3 = withETransdecoder.transcript_gff3,
      with_e_transdecoder_cdna_fasta = withETransdecoder.cdna_fasta,
      with_e_transdecoder_cds_fasta = withETransdecoder.cds_fasta,
      with_e_transdecoder_peptide_fasta = withETransdecoder.peptide_fasta,
      no_e_stringtie_gtf = stringtieNoE.transcript_gtf,
      no_e_transdecoder_genome_gff3 = noETransdecoder.genome_gff3,
      no_e_transdecoder_transcript_gff3 = noETransdecoder.transcript_gff3,
      no_e_transdecoder_cdna_fasta = noETransdecoder.cdna_fasta,
      no_e_transdecoder_cds_fasta = noETransdecoder.cds_fasta,
      no_e_transdecoder_peptide_fasta = noETransdecoder.peptide_fasta,
      canonical_pvacseq_outputs = canonical_pvacseq_outputs,
      docker_image = preparation_docker
  }

  call pvacseq.pvacseq as orfanagePredict {
    input:
      n_threads = n_threads,
      input_vcf = pvacncPrepare.pvacnc_input_vcf_gz,
      input_vcf_tbi = pvacncPrepare.pvacnc_input_vcf_tbi,
      sample_name = tumor_sample_name,
      alleles = alleles,
      prediction_algorithms = prediction_algorithms,
      peptide_fasta = peptide_fasta,
      genes_of_interest_file = genes_of_interest_file,
      epitope_lengths_class_i = epitope_lengths_class_i,
      epitope_lengths_class_ii = epitope_lengths_class_ii,
      binding_threshold = binding_threshold,
      binding_percentile_threshold = binding_percentile_threshold,
      presentation_percentile_threshold = presentation_percentile_threshold,
      immunogenicity_percentile_threshold = immunogenicity_percentile_threshold,
      percentile_threshold_strategy = percentile_threshold_strategy,
      iedb_retries = iedb_retries,
      normal_sample_name = normal_sample_name,
      net_chop_method = net_chop_method,
      top_score_metric = top_score_metric,
      top_score_metric2 = top_score_metric2,
      net_chop_threshold = net_chop_threshold,
      additional_report_columns = additional_report_columns,
      fasta_size = fasta_size,
      downstream_sequence_length = downstream_sequence_length,
      exclude_nas = exclude_nas,
      phased_proximal_variants_vcf = phased_proximal_variants_vcf,
      phased_proximal_variants_vcf_tbi = phased_proximal_variants_vcf_tbi,
      minimum_fold_change = minimum_fold_change,
      normal_cov = normal_cov,
      tdna_cov = tdna_cov,
      trna_cov = trna_cov,
      normal_vaf = normal_vaf,
      tdna_vaf = tdna_vaf,
      trna_vaf = trna_vaf,
      expn_val = expn_val,
      maximum_transcript_support_level = maximum_transcript_support_level,
      transcript_prioritization_strategy = transcript_prioritization_strategy,
      aggregate_inclusion_binding_threshold = aggregate_inclusion_binding_threshold,
      aggregate_inclusion_count_limit = aggregate_inclusion_count_limit,
      problematic_amino_acids = problematic_amino_acids,
      anchor_contribution_threshold = anchor_contribution_threshold,
      biotypes = biotypes,
      netmhciipan_version = netmhciipan_version,
      reference_scores_zip = reference_scores_zip,
      allele_specific_binding_thresholds = allele_specific_binding_thresholds,
      netmhc_stab = netmhc_stab,
      run_reference_proteome_similarity = run_reference_proteome_similarity,
      allele_specific_anchors = allele_specific_anchors,
      allow_incomplete_transcripts = allow_incomplete_transcripts,
      use_normalized_percentiles = use_normalized_percentiles,
      tumor_purity = tumor_purity,
      run_ml_predictions = run_ml_predictions,
      ml_threshold_accept = ml_threshold_accept,
      ml_threshold_reject = ml_threshold_reject
  }

  call pvacseq.pvacseq as transdecoderPredict {
    input:
      n_threads = n_threads,
      input_vcf = pvacncTransdecoderPrepare.pvacnc_input_vcf_gz,
      input_vcf_tbi = pvacncTransdecoderPrepare.pvacnc_input_vcf_tbi,
      sample_name = tumor_sample_name,
      alleles = alleles,
      prediction_algorithms = prediction_algorithms,
      peptide_fasta = peptide_fasta,
      genes_of_interest_file = genes_of_interest_file,
      epitope_lengths_class_i = epitope_lengths_class_i,
      epitope_lengths_class_ii = epitope_lengths_class_ii,
      binding_threshold = binding_threshold,
      binding_percentile_threshold = binding_percentile_threshold,
      presentation_percentile_threshold = presentation_percentile_threshold,
      immunogenicity_percentile_threshold = immunogenicity_percentile_threshold,
      percentile_threshold_strategy = percentile_threshold_strategy,
      iedb_retries = iedb_retries,
      normal_sample_name = normal_sample_name,
      net_chop_method = net_chop_method,
      top_score_metric = top_score_metric,
      top_score_metric2 = top_score_metric2,
      net_chop_threshold = net_chop_threshold,
      additional_report_columns = additional_report_columns,
      fasta_size = fasta_size,
      downstream_sequence_length = downstream_sequence_length,
      exclude_nas = exclude_nas,
      phased_proximal_variants_vcf = phased_proximal_variants_vcf,
      phased_proximal_variants_vcf_tbi = phased_proximal_variants_vcf_tbi,
      minimum_fold_change = minimum_fold_change,
      normal_cov = normal_cov,
      tdna_cov = tdna_cov,
      trna_cov = trna_cov,
      normal_vaf = normal_vaf,
      tdna_vaf = tdna_vaf,
      trna_vaf = trna_vaf,
      expn_val = expn_val,
      maximum_transcript_support_level = maximum_transcript_support_level,
      transcript_prioritization_strategy = transcript_prioritization_strategy,
      aggregate_inclusion_binding_threshold = aggregate_inclusion_binding_threshold,
      aggregate_inclusion_count_limit = aggregate_inclusion_count_limit,
      problematic_amino_acids = problematic_amino_acids,
      anchor_contribution_threshold = anchor_contribution_threshold,
      biotypes = biotypes,
      netmhciipan_version = netmhciipan_version,
      reference_scores_zip = reference_scores_zip,
      allele_specific_binding_thresholds = allele_specific_binding_thresholds,
      netmhc_stab = netmhc_stab,
      run_reference_proteome_similarity = run_reference_proteome_similarity,
      allele_specific_anchors = allele_specific_anchors,
      allow_incomplete_transcripts = allow_incomplete_transcripts,
      use_normalized_percentiles = use_normalized_percentiles,
      tumor_purity = tumor_purity,
      run_ml_predictions = run_ml_predictions,
      ml_threshold_accept = ml_threshold_accept,
      ml_threshold_reject = ml_threshold_reject
  }

  call view.pvacncAnnotatePvacview as orfanageAnnotatePvacview {
    input:
      mapping_tsv = pvacncPrepare.mapping_tsv,
      source_mhc_i = orfanagePredict.mhc_i,
      source_mhc_ii = orfanagePredict.mhc_ii,
      source_combined = orfanagePredict.combined,
      docker_image = preparation_docker
  }

  call view.pvacncAnnotatePvacview as transdecoderAnnotatePvacview {
    input:
      mapping_tsv = pvacncTransdecoderPrepare.mapping_tsv,
      source_mhc_i = transdecoderPredict.mhc_i,
      source_mhc_ii = transdecoderPredict.mhc_ii,
      source_combined = transdecoderPredict.combined,
      docker_image = preparation_docker
  }

  output {
    File no_e_stringtie_gtf = stringtieNoE.transcript_gtf

    File with_e_orfanage_gtf = withEOrfanage.best_gtf
    File no_e_orfanage_gtf = noEOrfanage.best_gtf
    File orfanage_input_vcf_gz = pvacncPrepare.pvacnc_input_vcf_gz
    File orfanage_input_vcf_tbi = pvacncPrepare.pvacnc_input_vcf_tbi
    File orfanage_mapping_tsv = pvacncPrepare.mapping_tsv
    File orfanage_routing_tsv = pvacncPrepare.routing_tsv
    File orfanage_effects_tsv = pvacncPrepare.effects_tsv
    Array[File] orfanage_preparation_outputs = pvacncPrepare.preparation_outputs
    Array[File] orfanage_mhc_i = orfanageAnnotatePvacview.mhc_i
    File? orfanage_mhc_i_log = orfanagePredict.mhc_i_log
    Array[File] orfanage_mhc_ii = orfanageAnnotatePvacview.mhc_ii
    File? orfanage_mhc_ii_log = orfanagePredict.mhc_ii_log
    Array[File] orfanage_combined = orfanageAnnotatePvacview.combined
    File orfanage_pvacview_annotation_log = orfanageAnnotatePvacview.annotation_log

    Array[File] with_e_transdecoder_outputs = [
      withETransdecoder.cdna_fasta,
      withETransdecoder.peptide_fasta,
      withETransdecoder.cds_fasta,
      withETransdecoder.transcript_gff3,
      withETransdecoder.genome_gff3,
      withETransdecoder.summary
    ]
    Array[File] no_e_transdecoder_outputs = [
      noETransdecoder.cdna_fasta,
      noETransdecoder.peptide_fasta,
      noETransdecoder.cds_fasta,
      noETransdecoder.transcript_gff3,
      noETransdecoder.genome_gff3,
      noETransdecoder.summary
    ]
    File transdecoder_input_vcf_gz = pvacncTransdecoderPrepare.pvacnc_input_vcf_gz
    File transdecoder_input_vcf_tbi = pvacncTransdecoderPrepare.pvacnc_input_vcf_tbi
    File transdecoder_mapping_tsv = pvacncTransdecoderPrepare.mapping_tsv
    File transdecoder_routing_tsv = pvacncTransdecoderPrepare.routing_tsv
    File transdecoder_effects_tsv = pvacncTransdecoderPrepare.effects_tsv
    Array[File] transdecoder_preparation_outputs = pvacncTransdecoderPrepare.preparation_outputs
    Array[File] transdecoder_mhc_i = transdecoderAnnotatePvacview.mhc_i
    File? transdecoder_mhc_i_log = transdecoderPredict.mhc_i_log
    Array[File] transdecoder_mhc_ii = transdecoderAnnotatePvacview.mhc_ii
    File? transdecoder_mhc_ii_log = transdecoderPredict.mhc_ii_log
    Array[File] transdecoder_combined = transdecoderAnnotatePvacview.combined
    File transdecoder_pvacview_annotation_log = transdecoderAnnotatePvacview.annotation_log
  }
}
