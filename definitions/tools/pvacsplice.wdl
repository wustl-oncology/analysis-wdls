version 1.0

workflow pvacsplice_workflow {
  input {
    Int n_threads = 8
    File input_vcf
    File input_vcf_tbi
    File input_regtools_tsv
    File input_reference_dna_fasta
    File input_reference_gtf
    String sample_name
    Array[String] alleles
    Array[String] prediction_algorithms
    File? peptide_fasta

    Array[Int]? epitope_lengths_class_i
    Array[Int]? epitope_lengths_class_ii
    Int? binding_threshold
    Int? percentile_threshold
    Int? iedb_retries

    String? normal_sample_name
    String? net_chop_method  # enum [cterm , 20s]
    String? top_score_metric  # enum [lowest, median]
    Float? net_chop_threshold
    String? additional_report_columns  # enum [sample_name]
    Int? fasta_size
    Boolean exclude_nas = false
    Int? normal_cov
    Int? tdna_cov
    Int? trna_cov
    Float? normal_vaf
    Float? tdna_vaf
    Float? trna_vaf
    Float? expn_val
    Int? maximum_transcript_support_level  # enum [1, 2, 3, 4, 5]
    Int? aggregate_inclusion_binding_threshold
    Array[String]? problematic_amino_acids
    Array[String]? biotypes
    Int? aggregate_inclusion_count_limit

    Int? junction_score
    Int? variant_distance
    Boolean save_gtf = false
    Array[String]? junction_anchor_types

    Boolean allele_specific_binding_thresholds = false
    Boolean keep_tmp_files = false
    Boolean netmhc_stab = false
    Boolean run_reference_proteome_similarity = false

    Float? tumor_purity
  }

  Float input_size = size([input_vcf, input_vcf_tbi, input_regtools_tsv, input_reference_dna_fasta, input_reference_gtf], "GB")
  Int space_needed_gb = 10 + round(input_size)

  call pvacsplice {
    input:
      n_threads = n_threads,
      input_vcf = input_vcf,
      input_vcf_tbi = input_vcf_tbi,
      input_regtools_tsv = input_regtools_tsv,
      input_reference_dna_fasta = input_reference_dna_fasta,
      input_reference_gtf = input_reference_gtf,
      sample_name = sample_name,
      alleles = alleles,
      prediction_algorithms = prediction_algorithms,
      peptide_fasta = peptide_fasta,
      epitope_lengths_class_i = epitope_lengths_class_i,
      epitope_lengths_class_ii = epitope_lengths_class_ii,
      binding_threshold = binding_threshold,
      percentile_threshold = percentile_threshold,
      iedb_retries = iedb_retries,
      normal_sample_name = normal_sample_name,
      net_chop_method = net_chop_method,
      top_score_metric = top_score_metric,
      net_chop_threshold = net_chop_threshold,
      additional_report_columns = additional_report_columns,
      fasta_size = fasta_size,
      exclude_nas = exclude_nas,
      normal_cov = normal_cov,
      tdna_cov = tdna_cov,
      trna_cov = trna_cov,
      normal_vaf = normal_vaf,
      tdna_vaf = tdna_vaf,
      trna_vaf = trna_vaf,
      expn_val = expn_val,
      maximum_transcript_support_level = maximum_transcript_support_level,
      aggregate_inclusion_binding_threshold = aggregate_inclusion_binding_threshold,
      problematic_amino_acids = problematic_amino_acids,
      biotypes = biotypes,
      aggregate_inclusion_count_limit = aggregate_inclusion_count_limit,
      junction_score = junction_score,
      variant_distance = variant_distance,
      save_gtf = save_gtf,
      junction_anchor_types = junction_anchor_types,
      allele_specific_binding_thresholds = allele_specific_binding_thresholds,
      keep_tmp_files = keep_tmp_files,
      netmhc_stab = netmhc_stab,
      run_reference_proteome_similarity = run_reference_proteome_similarity,
      tumor_purity = tumor_purity,
      space_needed_gb = space_needed_gb
  }

  output {
    Array[File] mhc_i = pvacsplice.mhc_i
    Array[File] mhc_ii = pvacsplice.mhc_ii
    Array[File] combined = pvacsplice.combined
    File? splice_transcript_combined_report = pvacsplice.splice_transcript_combined_report
    File? splice_fasta = pvacsplice.splice_fasta
    File? splice_fasta_fai = pvacsplice.splice_fasta_fai
  }
}
task pvacsplice {
  input {
    Int n_threads = 8
    File input_vcf
    File input_vcf_tbi
    File input_regtools_tsv
    File input_reference_dna_fasta
    File input_reference_gtf 
    String sample_name
    Array[String] alleles
    Array[String] prediction_algorithms
    File? peptide_fasta 

    Array[Int]? epitope_lengths_class_i
    Array[Int]? epitope_lengths_class_ii
    Int? binding_threshold
    Int? percentile_threshold
    Int? iedb_retries

    String? normal_sample_name
    String? net_chop_method  # enum [cterm , 20s]
    String? top_score_metric  # enum [lowest, median]
    Float? net_chop_threshold
    String? additional_report_columns  # enum [sample_name]
    Int? fasta_size
    Boolean exclude_nas = false
    Int? normal_cov
    Int? tdna_cov
    Int? trna_cov
    Float? normal_vaf
    Float? tdna_vaf
    Float? trna_vaf
    Float? expn_val
    Int? maximum_transcript_support_level  # enum [1, 2, 3, 4, 5]
    Int? aggregate_inclusion_binding_threshold
    Array[String]? problematic_amino_acids
    Array[String]? biotypes
    Int? aggregate_inclusion_count_limit
    
    Int? junction_score
    Int? variant_distance
    Boolean save_gtf = false
    Array[String]? junction_anchor_types

    Boolean allele_specific_binding_thresholds = false
    Boolean keep_tmp_files = false
    Boolean netmhc_stab = false
    Boolean run_reference_proteome_similarity = false

    Float? tumor_purity
    Int space_needed_gb
  }

  # Float input_size = size([input_vcf, input_vcf_tbi, input_regtools_tsv, input_reference_dna_fasta,input_reference_gtf], "GB") #input files: annotated vcf, regtools tsv, reference fasta, reference gtf 
  # Int space_needed_gb = 10 + round(input_size) 
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "32GB"
    cpu: n_threads
    docker: "griffithlab/pvactools:5.0.1"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  # explicit typing required, don't inline
  Array[Int] epitope_i = select_first([epitope_lengths_class_i, []])
  Array[Int] epitope_ii = select_first([epitope_lengths_class_ii, []])
  Array[String] problematic_aa = select_first([problematic_amino_acids, []])
  Array[String] biotypes_list = select_first([biotypes, []])
  Array[String] junction_anchor_types_list = select_first([junction_anchor_types,[]])
#   command <<<

#     # touch each tbi to ensure they have a timestamp after the vcf
#     touch ~{input_vcf_tbi}

#     ln -s "$TMPDIR" /tmp/pvacsplice && export TMPDIR=/tmp/pvacsplice && \
#     /usr/local/bin/pvacsplice run \
#     --iedb-install-directory /opt/iedb \
#     --pass-only \
#     ~{if defined(tumor_purity) then "--tumor-purity " + select_first([tumor_purity]) else ""} \
#     ~{if length(epitope_i) > 0 then "-e1 " + sep(",", epitope_i) else ""} \
#     ~{if length(epitope_ii) > 0 then "-e2 " + sep(",", epitope_ii) else ""} \
#     ~{if defined(binding_threshold) then "-b " + binding_threshold else ""} \
#     ~{if defined(percentile_threshold) then "--percentile-threshold " + percentile_threshold else ""} \
#     ~{if allele_specific_binding_thresholds then "--allele-specific-binding-thresholds" else ""} \
#     ~{if defined(aggregate_inclusion_binding_threshold) then "--aggregate-inclusion-binding-threshold " + aggregate_inclusion_binding_threshold else ""} \
#     ~{if defined(aggregate_inclusion_count_limit) then "--aggregate-inclusion-count-limit " + aggregate_inclusion_count_limit else ""} \
#     ~{if defined(iedb_retries) then "-r " + iedb_retries else ""} \
#     ~{if keep_tmp_files then "-k" else ""} \
#     ~{if defined(normal_sample_name) then "--normal-sample-name " + normal_sample_name else ""} \
#     ~{if defined(net_chop_method) then "--net-chop-method " + net_chop_method else ""} \
#     ~{if netmhc_stab then "--netmhc-stab" else ""} \
#     ~{if run_reference_proteome_similarity then "--run-reference-proteome-similarity" else ""} \
#     ~{if defined(peptide_fasta) then "--peptide-fasta " + peptide_fasta else ""} \
#     ~{if defined(top_score_metric) then "-m " + top_score_metric else ""} \
#     ~{if defined(net_chop_threshold) then "--net-chop-threshold " + net_chop_threshold else ""} \
#     ~{if defined(additional_report_columns) then "-a " + additional_report_columns else ""} \
#     ~{if defined(fasta_size) then "-s " + fasta_size else ""} \
#     ~{if exclude_nas then "--exclude-NAs" else ""} \
#     ~{if defined(normal_cov) then "--normal-cov " + normal_cov else ""} \
#     ~{if defined(tdna_cov) then "--tdna-cov " + tdna_cov else ""} \
#     ~{if defined(trna_cov) then "--trna-cov " + trna_cov else ""} \
#     ~{if defined(normal_vaf) then "--normal-vaf " + normal_vaf else ""} \
#     ~{if defined(tdna_vaf) then "--tdna-vaf " + tdna_vaf else ""} \
#     ~{if defined(trna_vaf) then "--trna-vaf " + trna_vaf else ""} \
#     ~{if defined(expn_val) then "--expn-val " + expn_val else ""} \
#     ~{if defined(maximum_transcript_support_level) then "--maximum-transcript-support-level " + maximum_transcript_support_level else ""} \
#     ~{if length(problematic_aa) > 0 then "--problematic-amino-acids " + sep(",", problematic_aa) else ""} \
#     ~{if length(biotypes_list) > 0 then "--biotypes " + sep(",", biotypes_list) else ""} \
#     ~{if defined(junction_score) then "--junction-score " + junction_score else ""} \
#     ~{if defined(variant_distance) then "--variant-distance " + variant_distance else ""} \
#     ~{if save_gtf then "-g" else ""} \
#     ~{if length(junction_anchor_types_list) > 0 then "--anchor-types " + sep(",", junction_anchor_types_list) else ""} \
#     -t ~{n_threads} \
#     ~{input_regtools_tsv} \
#     ~{sample_name} \
#     ~{sep="," alleles} \
#     ~{sep=" " prediction_algorithms} \
#     pvacsplice_predictions \
#     ~{input_vcf} \
#     ~{input_reference_dna_fasta} \
#     ~{input_reference_gtf}

# >>>

  command <<<
    # touch each tbi to ensure they have a timestamp after the vcf
    touch ~{input_vcf_tbi}

    ln -s "$TMPDIR" /tmp/pvacsplice && export TMPDIR=/tmp/pvacsplice && \
    /usr/local/bin/pvacsplice run --iedb-install-directory /opt/iedb \
    --pass-only \
    ~{if defined(tumor_purity) then "--tumor-purity " + select_first([tumor_purity]) else ""} \
    ~{if length(epitope_i ) > 0 then "-e1 " else ""} ~{sep="," epitope_i} \
    ~{if length(epitope_ii) > 0 then "-e2 " else ""} ~{sep="," epitope_ii} \
    ~{if defined(binding_threshold) then "-b ~{binding_threshold}" else ""} \
    ~{if defined(percentile_threshold) then "--percentile-threshold ~{percentile_threshold}" else ""} \
    ~{if allele_specific_binding_thresholds then "--allele-specific-binding-thresholds" else ""} \
    ~{if defined(aggregate_inclusion_binding_threshold) then "--aggregate-inclusion-binding-threshold ~{aggregate_inclusion_binding_threshold}" else ""} \
    ~{if defined(aggregate_inclusion_count_limit) then "--aggregate-inclusion-count-limit ~{aggregate_inclusion_count_limit}" else ""} \
    ~{if defined(iedb_retries) then "-r ~{iedb_retries}" else ""} \
    ~{if keep_tmp_files then "-k" else ""} \
    ~{if defined(normal_sample_name) then "--normal-sample-name ~{normal_sample_name}" else ""} \
    ~{if defined(net_chop_method) then "--net-chop-method ~{net_chop_method}" else ""} \
    ~{if netmhc_stab then "--netmhc-stab" else ""} \
    ~{if run_reference_proteome_similarity then "--run-reference-proteome-similarity" else ""} \
    ~{if defined(peptide_fasta) then "--peptide-fasta ~{peptide_fasta}" else ""} \
    ~{if defined(top_score_metric) then "-m ~{top_score_metric}" else ""} \
    ~{if defined(net_chop_threshold) then "--net-chop-threshold ~{net_chop_threshold}" else ""} \
    ~{if defined(additional_report_columns) then "-a ~{additional_report_columns}" else ""} \
    ~{if defined(fasta_size) then "-s ~{fasta_size}" else ""} \
    ~{if exclude_nas then "--exclude-NAs" else ""} \
    ~{if defined(normal_cov) then "--normal-cov ~{normal_cov}" else ""} \
    ~{if defined(tdna_cov) then "--tdna-cov ~{tdna_cov}" else ""} \
    ~{if defined(trna_cov) then "--trna-cov ~{trna_cov}" else ""} \
    ~{if defined(normal_vaf) then "--normal-vaf ~{normal_vaf}" else ""} \
    ~{if defined(tdna_vaf) then "--tdna-vaf ~{tdna_vaf}" else ""} \
    ~{if defined(trna_vaf) then "--trna-vaf ~{trna_vaf}" else ""} \
    ~{if defined(expn_val) then "--expn-val ~{expn_val}" else ""} \
    ~{if defined(maximum_transcript_support_level) then "--maximum-transcript-support-level ~{maximum_transcript_support_level}" else ""} \
    ~{if length(problematic_aa) > 0 then "--problematic-amino-acids" else ""} ~{sep="," problematic_aa} \
    ~{if length(biotypes_list) > 0 then "--biotypes" else ""} ~{sep="," biotypes_list} \
    ~{if defined(junction_score) then "--junction-score ~{junction_score}" else ""} \
    ~{if defined(variant_distance) then "--variant-distance ~{variant_distance}" else ""} \
    ~{if save_gtf then "-g" else ""} \
    ~{if length(junction_anchor_types_list) > 0 then "--anchor-types" else ""} ~{sep="," junction_anchor_types_list} \
    -t ~{n_threads} \
    ~{input_regtools_tsv} ~{sample_name} ~{sep="," alleles} ~{sep=" " prediction_algorithms} \
    pvacsplice_predictions ~{input_vcf} ~{input_reference_dna_fasta} ~{input_reference_gtf} 
  >>>

  output {
    File? mhc_i_all_epitopes = "pvacsplice_predictions/MHC_Class_I/~{sample_name}.all_epitopes.tsv"
    File? mhc_i_aggregated_report = "pvacsplice_predictions/MHC_Class_I/~{sample_name}.all_epitopes.aggregated.tsv"
    File? mhc_i_filtered_epitopes = "pvacsplice_predictions/MHC_Class_I/~{sample_name}.filtered.tsv"
    File? mhc_ii_all_epitopes = "pvacsplice_predictions/MHC_Class_II/~{sample_name}.all_epitopes.tsv"
    File? mhc_ii_aggregated_report = "pvacsplice_predictions/MHC_Class_II/~{sample_name}.all_epitopes.aggregated.tsv"
    File? mhc_ii_filtered_epitopes = "pvacsplice_predictions/MHC_Class_II/~{sample_name}.filtered.tsv"
    File? combined_all_epitopes = "pvacsplice_predictions/combined/~{sample_name}.all_epitopes.tsv"
    File? combined_aggregated_report = "pvacsplice_predictions/combined/~{sample_name}.all_epitopes.aggregated.tsv"
    File? combined_filtered_epitopes = "pvacsplice_predictions/combined/~{sample_name}.filtered.tsv"
    File? splice_transcript_combined_report = "pvacsplice_predictions/~{sample_name}_combined.tsv"
    File? splice_gtf = "pvacsplice_predictions/~{sample_name}_gtf.tsv"
    File? splice_fasta = "pvacsplice_predictions/~{sample_name}_.transcripts.fa"
    File? splice_fasta_fai = "pvacsplice_predictions/~{sample_name}_.transcripts.fa.fai"
    
    # glob documentations
    # https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#globs
    Array[File] mhc_i = glob("pvacsplice_predictions/MHC_Class_I/*")
    Array[File] mhc_ii = glob("pvacsplice_predictions/MHC_Class_II/*")
    Array[File] combined = glob("pvacsplice_predictions/combined/*")

  }
}
