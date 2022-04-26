version 1.0

task pvacfuse {
  input {
    File input_fusions_zip
    String sample_name
    Array[String] alleles
    Array[String] prediction_algorithms
    Array[Int]? epitope_lengths_class_i
    Array[Int]? epitope_lengths_class_ii
    Int? binding_threshold
    Int? percentile_threshold
    Int? iedb_retries
    Boolean keep_tmp_files = false
    String? net_chop_method  # enum [cterm 20s]
    Boolean netmhc_stab = false
    String? top_score_metric  # enum [lowest, median]
    Float? net_chop_threshold
    Boolean run_reference_proteome_similarity = false
    String? blastp_db  # enum [refseq_select_prot, refseq_protein]
    String? additional_report_columns  # enum [sample_name]
    Int? fasta_size
    Int? downstream_sequence_length
    Boolean exclude_nas = false
    Int n_threads = 8
  }

  Int space_needed_gb = 10 + round(size([input_fusions_zip], "GB") * 3)
  runtime {
    docker: "griffithlab/pvactools:3.0.0"
    memory: "16GB"
    cpu: n_threads
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  # explicit typing required, don't inline
  Array[Int] epitope_i = select_first([epitope_lengths_class_i, []])
  Array[Int] epitope_ii = select_first([epitope_lengths_class_ii, []])
  command <<<
    mkdir agfusion_dir && unzip -qq ~{input_fusions_zip} -d agfusion_dir

    ln -s "$TMPDIR" /tmp/pvacseq && export TMPDIR=/tmp/pvacseq && \
    /usr/local/bin/pvacfuse run --iedb-install-directory /opt/iedb \
    --blastp-path /opt/ncbi-lbast-2.12.0+/bin/blastp \
    ~{if defined(blastp_db) then "--blastp-db " + select_first([blastp_db]) else ""} \
    agfusion_dir ~{sample_name} \
    ~{sep="," alleles} \
    ~{sep=" " prediction_algorithms} \
    pvacfuse_predictions \
    ~{if defined(epitope_lengths_class_i ) then "-e1 " else ""} ~{sep="," epitope_i} \
    ~{if defined(epitope_lengths_class_ii) then "-e2 " else ""} ~{sep="," epitope_ii} \
    ~{if defined(binding_threshold) then "-b ~{binding_threshold}" else ""} \
    ~{if defined(percentile_threshold) then "--percentile-threshold ~{percentile_threshold}" else ""} \
    ~{if defined(iedb_retries) then "-r ~{iedb_retries}" else ""} \
    ~{if keep_tmp_files then "-k" else ""} \
    ~{if defined(net_chop_method) then "--net-chop-method ~{net_chop_method}" else ""} \
    ~{if netmhc_stab then "--netmhc-stab" else ""} \
    ~{if defined(top_score_metric) then "-m ~{top_score_metric}" else ""} \
    ~{if defined(top_score_metric) then "-m ~{top_score_metric}" else ""} \
    ~{if defined(net_chop_threshold) then "--net-chop-threshold ~{net_chop_threshold}" else ""} \
    ~{if run_reference_proteome_similarity then "--run-reference-proteome-similarity" else ""} \
    ~{if defined(additional_report_columns) then "-m ~{additional_report_columns}" else ""} \
    ~{if defined(fasta_size) then "-s ~{fasta_size}" else ""} \
    ~{if defined(downstream_sequence_length) then "-d ~{downstream_sequence_length}" else ""} \
    ~{if exclude_nas then "--exclude-NAs" else ""} \
    --n-threads ~{n_threads}
  >>>

  output {
    File? mhc_i_all_epitopes = "pvacseq_predictions/MHC_Class_I/~{sample_name}.all_epitopes.tsv"
    File? mhc_i_aggregated_report = "pvacseq_predictions/MHC_Class_I/~{sample_name}.all_epitopes.aggregated.tsv"
    File? mhc_i_filtered_epitopes = "pvacseq_predictions/MHC_Class_I/~{sample_name}.filtered.tsv"
    File? mhc_ii_all_epitopes = "pvacseq_predictions/MHC_Class_II/~{sample_name}.all_epitopes.tsv"
    File? mhc_ii_aggregated_report = "pvacseq_predictions/MHC_Class_II/~{sample_name}.all_epitopes.aggregated.tsv"
    File? mhc_ii_filtered_epitopes = "pvacseq_predictions/MHC_Class_II/~{sample_name}.filtered.tsv"
    File? combined_all_epitopes = "pvacseq_predictions/combined/~{sample_name}.all_epitopes.tsv"
    File? combined_aggregated_report = "pvacseq_predictions/combined/~{sample_name}.all_epitopes.aggregated.tsv"
    File? combined_filtered_epitopes = "pvacseq_predictions/combined/~{sample_name}.filtered.tsv"
    Array[File] pvacfuse_predictions = glob("pvacseq_predictions/**/*")
  }
}
