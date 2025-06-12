version 1.0
import "pvacseq.wdl" as pvacseq
import "rnaseq.wdl" as rnaseq
import "somatic_exome_nonhuman.wdl" as sen

workflow nonhuman_immuno {
    input {
        # Reference Files
        File reference
        File reference_fai
        File reference_dict
        File reference_alt
        File reference_amb
        File reference_ann
        File reference_bwt
        File reference_pac
        File reference_0123
        File reference_index
        File reference_index_1ht2
        File reference_index_2ht2
        File reference_index_3ht2
        File reference_index_4ht2
        File reference_index_5ht2
        File reference_index_6ht2
        File reference_index_7ht2
        File reference_index_8ht2
        File reference_annotation

        # Rnaseq Inputs
        Array[String] read_group_id
        Array[SequenceData] rna_sequence
        Array[Array[String]] read_group_fields
        String? strand  # [first, second, unstranded]
        String sample_name = "TUMOR"
        String normal_sample_name = "NORMAL"

        File trimming_adapters
        String trimming_adapter_trim_end
        Int trimming_adapter_min_overlap
        Int trimming_max_uncalled
        Int trimming_min_readlength

        File kallisto_index
        File gene_transcript_lookup_table
        File refFlat
        File? ribosomal_intervals
        # Pvacseq Inputs
        # File detect_variants_vcf
        # File detect_variants_vcf_tbi
        # String sample_name = "TUMOR"
        # File rnaseq_bam
        # File rnaseq_bam_bai
        # File reference
        # File reference_fai
        # File reference_dict
        File? peptide_fasta
        Int? readcount_minimum_base_quality
        Int? readcount_minimum_mapping_quality
        # File gene_expression_file
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
        String? additional_report_columns  # enum [sample_name]
        Int? fasta_size
        Int? downstream_sequence_length
        Boolean? exclude_nas
        File? phased_proximal_variants_vcf
        File? phased_proximal_variants_vcf_tbi
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
        Array[String] pvacseq_variants_to_table_fields = ["CHROM", "POS", "ID", "REF", "ALT"]
        Array[String] pvacseq_variants_to_table_genotype_fields = ["GT", "AD", "AF", "DP", "RAD", "RAF", "RDP", "GX", "TX"]
        Array[String] pvacseq_vep_to_table_fields = ["HGVSc", "HGVSp"]
        Float? tumor_purity
        Boolean? allele_specific_binding_thresholds
        Int? aggregate_inclusion_binding_threshold
        Int? aggregate_inclusion_count_limit
        Array[String]? problematic_amino_acids
        Boolean? allele_specific_anchors
        Float? anchor_contribution_threshold
        String? prefix = "pvacseq"
        Array[String]? biotypes
        # Somatic Exome Nonhuman Inputs:
        # File reference
        # File reference_fai
        # File reference_dict
        

        # String tumor_name = "tumor"
        Array[SequenceData] tumor_sequence

        # String normal_name = "normal"
        Array[SequenceData] normal_sequence

        TrimmingOptions? trimming

        File bait_intervals
        File target_intervals
        Int target_interval_padding = 100
        Array[LabelledFile] per_base_intervals
        Array[LabelledFile] per_target_intervals
        Array[LabelledFile] summary_intervals

        String picard_metric_accumulation_level
        Int qc_minimum_mapping_quality = 0
        Int qc_minimum_base_quality = 0

        Int strelka_cpu_reserved = 8
        Int scatter_count = 50

        Int? varscan_strand_filter
        Int? varscan_min_coverage
        Float? varscan_min_var_freq
        Float? varscan_p_value
        Float? varscan_max_normal_freq

        Float? fp_min_var_freq

        File vep_cache_dir_zip
        String vep_ensembl_assembly
        String vep_ensembl_version
        String vep_ensembl_species
        File? synonyms_file
        Boolean annotate_coding_only = false
        # one of [pick, flag_pick, pick-allele, per_gene, pick_allele_gene, flag_pick_allele, flag_pick_allele_gene]
        String? vep_pick
        Boolean cle_vcf_filter = false

        Float? filter_somatic_llr_threshold
        Float? filter_somatic_llr_tumor_purity
        Float? filter_somatic_llr_normal_contamination_rate

        # Array[String] vep_to_table_fields = ["Consequence", "SYMBOL", "Feature"]
        # Array[String] variants_to_table_genotype_fields = ["GT", "AD"]
        # Array[String] variants_to_table_fields = ["CHROM", "POS", "ID", "REF", "ALT", "set", "AC", "AF"]

        # String tumor_sample_name
        # String normal_sample_name

        Int? max_mm_qualsum_diff
        Int? max_var_mm_qualsum
        Array[String] somatic_exome_variants_to_table_fields = ["CHROM", "POS", "ID", "REF", "ALT"]
        Array[String] somatic_exome_variants_to_table_genotype_fields = ["GT", "AD", "AF", "DP", "RAD", "RAF", "RDP", "GX", "TX"]
        Array[String] somatic_exome_vep_to_table_fields = ["HGVSc", "HGVSp"]
    }

    String rna_sample_name = "RNA_" + sample_name

    
    call rnaseq.rnaseq as rns {
        input:
        reference=reference,
        reference_fai=reference_fai,
        reference_dict=reference_dict,
        reference_index=reference_index,
        reference_index_1ht2=reference_index_1ht2,
        reference_index_2ht2=reference_index_2ht2,
        reference_index_3ht2=reference_index_3ht2,
        reference_index_4ht2=reference_index_4ht2,
        reference_index_5ht2=reference_index_5ht2,
        reference_index_6ht2=reference_index_6ht2,
        reference_index_7ht2=reference_index_7ht2,
        reference_index_8ht2=reference_index_8ht2,
        reference_annotation=reference_annotation,
        read_group_id=read_group_id,
        rna_sequence=rna_sequence,
        read_group_fields=read_group_fields,
        strand=strand,
        sample_name=rna_sample_name,
        trimming_adapters=trimming_adapters,
        trimming_adapter_trim_end=trimming_adapter_trim_end,
        trimming_adapter_min_overlap=trimming_adapter_min_overlap,
        trimming_max_uncalled=trimming_max_uncalled,
        trimming_min_readlength=trimming_min_readlength,
        kallisto_index=kallisto_index,
        gene_transcript_lookup_table=gene_transcript_lookup_table,
        refFlat=refFlat,
        ribosomal_intervals=ribosomal_intervals
    }
    call sen.somaticExomeNonhuman as somat {
        input:
        reference=reference,
        reference_fai=reference_fai,
        reference_dict=reference_dict,
        reference_alt=reference_alt,
        reference_amb=reference_amb,
        reference_ann=reference_ann,
        reference_bwt=reference_bwt,
        reference_pac=reference_pac,
        reference_0123=reference_0123,
        tumor_name=sample_name,
        tumor_sequence=tumor_sequence,
        normal_name=normal_sample_name,
        normal_sequence=normal_sequence,
        trimming=trimming,
        bait_intervals=bait_intervals,
        target_intervals=target_intervals,
        target_interval_padding=target_interval_padding,
        per_base_intervals=per_base_intervals,
        per_target_intervals=per_target_intervals,
        summary_intervals=summary_intervals,
        picard_metric_accumulation_level=picard_metric_accumulation_level,
        qc_minimum_mapping_quality=qc_minimum_mapping_quality,
        qc_minimum_base_quality=qc_minimum_base_quality,
        strelka_cpu_reserved=strelka_cpu_reserved,
        scatter_count=scatter_count,
        varscan_strand_filter=varscan_strand_filter,
        varscan_min_coverage=varscan_min_coverage,
        varscan_min_var_freq=varscan_min_var_freq,
        varscan_p_value=varscan_p_value,
        varscan_max_normal_freq=varscan_max_normal_freq,
        fp_min_var_freq=fp_min_var_freq,
        vep_cache_dir_zip=vep_cache_dir_zip,
        vep_ensembl_assembly=vep_ensembl_assembly,
        vep_ensembl_version=vep_ensembl_version,
        vep_ensembl_species=vep_ensembl_species,
        synonyms_file=synonyms_file,
        annotate_coding_only=annotate_coding_only,
        vep_pick=vep_pick,
        cle_vcf_filter=cle_vcf_filter,
        filter_somatic_llr_threshold=filter_somatic_llr_threshold,
        filter_somatic_llr_tumor_purity=filter_somatic_llr_tumor_purity,
        filter_somatic_llr_normal_contamination_rate=filter_somatic_llr_normal_contamination_rate,
        vep_to_table_fields=somatic_exome_vep_to_table_fields,
        variants_to_table_genotype_fields=somatic_exome_variants_to_table_genotype_fields,
        variants_to_table_fields=somatic_exome_variants_to_table_fields,
        tumor_sample_name=sample_name,
        normal_sample_name=normal_sample_name,
        max_mm_qualsum_diff=max_mm_qualsum_diff,
        max_var_mm_qualsum=max_var_mm_qualsum
    }
    call pvacseq.pvacseq as pvs {
        input:
        detect_variants_vcf=somat.final_vcf,
        detect_variants_vcf_tbi=somat.final_vcf_tbi,
        sample_name=sample_name,
        normal_sample_name=normal_sample_name,
        rnaseq_bam=rns.final_bam,
        rnaseq_bam_bai=rns.final_bam_bai,
        reference=reference,
        reference_fai=reference_fai,
        reference_dict=reference_dict,
        peptide_fasta=peptide_fasta,
        readcount_minimum_base_quality=readcount_minimum_base_quality,
        readcount_minimum_mapping_quality=readcount_minimum_mapping_quality,
        gene_expression_file=rns.kallisto_gene_abundance,
        transcript_expression_file=rns.kallisto_transcript_abundance_tsv, # If this is wrong, rns.stringtie_transcript_expression_tsv
        expression_tool=expression_tool,
        alleles=alleles,
        prediction_algorithms=prediction_algorithms,
        epitope_lengths_class_i=epitope_lengths_class_i,
        epitope_lengths_class_ii=epitope_lengths_class_ii,
        binding_threshold=binding_threshold,
        percentile_threshold=percentile_threshold,
        percentile_threshold_strategy=percentile_threshold_strategy,
        minimum_fold_change=minimum_fold_change,
        top_score_metric=top_score_metric,
        additional_report_columns=additional_report_columns,
        fasta_size=fasta_size,
        downstream_sequence_length=downstream_sequence_length,
        exclude_nas=exclude_nas,
        phased_proximal_variants_vcf=phased_proximal_variants_vcf,
        phased_proximal_variants_vcf_tbi=phased_proximal_variants_vcf_tbi,
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
        n_threads=n_threads,
        iedb_retries=iedb_retries,
        variants_to_table_fields=pvacseq_variants_to_table_fields,
        variants_to_table_genotype_fields=pvacseq_variants_to_table_genotype_fields,
        vep_to_table_fields=pvacseq_vep_to_table_fields,
        tumor_purity=tumor_purity,
        allele_specific_binding_thresholds=allele_specific_binding_thresholds,
        aggregate_inclusion_binding_threshold=aggregate_inclusion_binding_threshold,
        aggregate_inclusion_count_limit=aggregate_inclusion_count_limit,
        problematic_amino_acids=problematic_amino_acids,
        allele_specific_anchors=allele_specific_anchors,
        anchor_contribution_threshold=anchor_contribution_threshold,
        prefix=prefix,
        biotypes=biotypes
    }
    output {
        #PVacseq:
        File annotated_vcf = pvs.annotated_vcf
        File annotated_vcf_tbi = pvs.annotated_vcf_tbi
        File annotated_tsv = pvs.annotated_tsv
        Array[File] mhc_i = pvs.mhc_i
        File? mhc_i_log = pvs.mhc_i_log
        Array[File] mhc_ii = pvs.mhc_ii
        File? mhc_ii_log = pvs.mhc_ii_log
        Array[File] combined = pvs.combined
        File indel_counting_bam = pvs.indel_counting_bam
        File indel_counting_bai = pvs.indel_counting_bai

        # RnaSeq:
        File final_bam = rns.final_bam
        File final_bam_bai = rns.final_bam_bai 
        File stringtie_transcript_gtf = rns.stringtie_transcript_gtf
        File stringtie_gene_expression_tsv = rns.stringtie_gene_expression_tsv
        File kallisto_transcript_abundance_tsv = rns.kallisto_transcript_abundance_tsv
        File kallisto_transcript_abundance_h5 = rns.kallisto_transcript_abundance_h5
        File kallisto_gene_abundance = rns.kallisto_gene_abundance
        File metrics = rns.metrics
        File? chart = rns.chart
        File kallisto_fusion_evidence = rns.kallisto_fusion_evidence
        File bamcoverage_bigwig = rns.bamcoverage_bigwig

        # Somatic Exome:
        File tumor_cram = somat.tumor_cram
        File tumor_cram_crai = somat.tumor_cram_crai
        File tumor_mark_duplicates_metrics = somat.tumor_mark_duplicates_metrics
        QCMetrics tumor_qc_metrics = somat.tumor_qc_metrics
        File normal_cram = somat.normal_cram
        File normal_cram_crai = somat.normal_cram_crai
        File normal_mark_duplicates_metrics = somat.normal_mark_duplicates_metrics
        QCMetrics normal_qc_metrics = somat.normal_qc_metrics
        File mutect_unfiltered_vcf = somat.mutect_unfiltered_vcf
        File mutect_unfiltered_vcf_tbi = somat.mutect_unfiltered_vcf_tbi
        File mutect_filtered_vcf = somat.mutect_filtered_vcf
        File mutect_filtered_vcf_tbi = somat.mutect_filtered_vcf_tbi
        File strelka_unfiltered_vcf = somat.strelka_unfiltered_vcf
        File strelka_unfiltered_vcf_tbi = somat.strelka_unfiltered_vcf_tbi
        File strelka_filtered_vcf = somat.strelka_filtered_vcf
        File strelka_filtered_vcf_tbi = somat.strelka_filtered_vcf_tbi
        File varscan_unfiltered_vcf = somat.varscan_unfiltered_vcf
        File varscan_unfiltered_vcf_tbi = somat.varscan_unfiltered_vcf_tbi
        File varscan_filtered_vcf = somat.varscan_filtered_vcf
        File varscan_filtered_vcf_tbi = somat.varscan_filtered_vcf_tbi
        File final_vcf = somat.final_vcf
        File final_vcf_tbi = somat.final_vcf_tbi
        File final_filtered_vcf = somat.final_filtered_vcf
        File final_filtered_vcf_tbi = somat.final_filtered_vcf_tbi
        File final_tsv = somat.final_tsv
        File vep_summary = somat.vep_summary
        File tumor_snv_bam_readcount_tsv = somat.tumor_snv_bam_readcount_tsv
        File tumor_indel_bam_readcount_tsv = somat.tumor_indel_bam_readcount_tsv
        File normal_snv_bam_readcount_tsv = somat.normal_snv_bam_readcount_tsv
        File normal_indel_bam_readcount_tsv = somat.normal_indel_bam_readcount_tsv

    }
}