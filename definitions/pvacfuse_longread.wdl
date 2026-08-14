version 1.0

import "tools/ctat_lr_fusion.wdl" as clf
import "tools/agfusion.wdl" as af
import "tools/pvacfuse.wdl" as pf
import "tools/zip_directory.wdl" as zd

# Starting point of this workflow is a long-read RNA BAM file from PacBio or Nanopore
# The workflow will:
# 1. Detect fusions using CTAT-LR-Fusion
# 2. Annotate fusions using AGfusion
# 3. Predict neoepitopes using pVACfuse
#
# TWO ENTRY POINTS
#
#   (a) Default -- supply lr_bam. Runs the full chain above.
#
#   (b) Resume -- supply agfusion_dir instead. Steps 1 and 2 are skipped and
#       pVACfuse runs straight off an AGFusion output directory you already
#       have. CTAT-LR-Fusion is by far the most expensive step here, so this is
#       the path to use when you want to re-run pVACfuse against different HLA
#       alleles, epitope lengths or binding thresholds.
#
#       lr_bam is still required by WDL's type system when using (b); it is
#       simply never read. Point it at the same bam or any placeholder.
#
# On the zip: pVACfuse itself always consumes a DIRECTORY. tools/pvacfuse.wdl
# unzips into `agfusion_dir` before invoking it -- the zip is only a transport
# wrapper for crossing task boundaries, since WDL 1.0 has no portable Directory
# type. Entry point (b) therefore re-zips your directory rather than adding a
# second pVACfuse task, which keeps tools/pvacfuse.wdl unmodified and leaves
# every output type below unchanged.

workflow pvacfuse_longread {
  input {
    # Long-read RNA BAM inputs
    File lr_bam
    File lr_bam_bai
    File star_fusion_genome_lib_zip
    Int ctat_cpu = 30
    Boolean examine_coding_effect = true
    Boolean vis = true

    # Entry point (b): an existing AGFusion output directory. When set,
    # CTAT-LR-Fusion and AGFusion are both skipped.
    Directory? agfusion_dir

    # AGfusion inputs
    File agfusion_database
    Boolean agfusion_annotate_noncanonical = true

    # HLA typing inputs
    Array[String] alleles

    # pVACfuse inputs
    String sample_name = "TUMOR"
    Array[String] prediction_algorithms
    Array[Int]? epitope_lengths_class_i
    Array[Int]? epitope_lengths_class_ii
    Int? binding_threshold
    Int? percentile_threshold
    String? percentile_threshold_strategy
    Int? iedb_retries
    Boolean pvacfuse_keep_tmp_files = false
    String? net_chop_method  # enum [cterm 20s]
    Boolean netmhc_stab = false
    String? top_score_metric  # enum [lowest, median]
    String? top_score_metric2  # enum [ic50, percentile]
    Float? net_chop_threshold
    Boolean run_reference_proteome_similarity = false
    String? additional_report_columns  # enum [sample_name]
    Int? fasta_size
    Int? downstream_sequence_length
    Boolean exclude_nas = false
    Int pvacfuse_threads = 8
    Boolean allele_specific_binding_thresholds = false
    Int? aggregate_inclusion_binding_threshold
    Int? aggregate_inclusion_count_limit
    Array[String]? problematic_amino_acids
    File? peptide_fasta
    File? genes_of_interest_file
    Int? read_support
    Float? expn_val
    String? netmhciipan_version # enum [4.3, 4.2, 4.1, 4.0]
  }

  # Steps 1 and 2 run only on entry point (a).
  if (!defined(agfusion_dir)) {
    # Step 1: Detect fusions using CTAT-LR-Fusion
    call clf.ctat_lr_fusion {
      input:
      lr_bam=lr_bam,
      star_fusion_genome_lib_zip=star_fusion_genome_lib_zip,
      cpu=ctat_cpu,
      examine_coding_effect=examine_coding_effect,
      vis=vis
    }

    # Step 2: Annotate fusions using AGfusion
    call af.agfusion {
      input:
      fusion_predictions=ctat_lr_fusion.fusion_predictions_abridged,
      agfusion_database=agfusion_database,
      annotate_noncanonical=agfusion_annotate_noncanonical
    }
  }

  # Entry point (b): package the supplied directory the same way agfusion.wdl
  # would have, so the single pVACfuse call below does not care which path ran.
  if (defined(agfusion_dir)) {
    call zd.zipDirectory as zipAgfusionDir {
      input:
      dir=select_first([agfusion_dir]),
      output_basename="agfusion_results"
    }
  }

  File fusions_zip_to_use = select_first([agfusion.annotated_fusion_predictions_zip, zipAgfusionDir.zipped])

  # Step 3: Predict neoepitopes using pVACfuse
  call pf.pvacfuse {
    input:
    input_fusions_zip=fusions_zip_to_use,
    sample_name=sample_name,
    alleles=alleles,
    prediction_algorithms=prediction_algorithms,
    epitope_lengths_class_i=epitope_lengths_class_i,
    epitope_lengths_class_ii=epitope_lengths_class_ii,
    binding_threshold=binding_threshold,
    percentile_threshold=percentile_threshold,
    percentile_threshold_strategy=percentile_threshold_strategy,
    iedb_retries=iedb_retries,
    keep_tmp_files=pvacfuse_keep_tmp_files,
    net_chop_method=net_chop_method,
    netmhc_stab=netmhc_stab,
    top_score_metric=top_score_metric,
    top_score_metric2=top_score_metric2,
    net_chop_threshold=net_chop_threshold,
    run_reference_proteome_similarity=run_reference_proteome_similarity,
    additional_report_columns=additional_report_columns,
    fasta_size=fasta_size,
    downstream_sequence_length=downstream_sequence_length,
    exclude_nas=exclude_nas,
    n_threads=pvacfuse_threads,
    allele_specific_binding_thresholds=allele_specific_binding_thresholds,
    aggregate_inclusion_binding_threshold=aggregate_inclusion_binding_threshold,
    aggregate_inclusion_count_limit=aggregate_inclusion_count_limit,
    problematic_amino_acids=problematic_amino_acids,
    peptide_fasta=peptide_fasta,
    genes_of_interest_file=genes_of_interest_file,
    read_support=read_support,
    expn_val=expn_val,
    netmhciipan_version=netmhciipan_version
  }

  output {
    # Fusion detection outputs.
    # Optional: absent on entry point (b), where CTAT-LR-Fusion never ran.
    File? ctat_fusion_predictions = ctat_lr_fusion.fusion_predictions
    File? ctat_fusion_predictions_abridged = ctat_lr_fusion.fusion_predictions_abridged
    Directory? ctat_lr_fusion_outdir = ctat_lr_fusion.ctat_lr_fusion_outdir

    # AGfusion outputs. On entry point (b) this is the re-zipped copy of the
    # directory you supplied, not a fresh AGFusion run.
    File agfusion_annotated_predictions_zip = fusions_zip_to_use

    # pVACfuse outputs - MHC Class I
    File? mhc_i_all_epitopes = pvacfuse.mhc_i_all_epitopes
    File? mhc_i_aggregated_report = pvacfuse.mhc_i_aggregated_report
    File? mhc_i_filtered_epitopes = pvacfuse.mhc_i_filtered_epitopes
    File? mhc_i_log = pvacfuse.mhc_i_log

    # pVACfuse outputs - MHC Class II
    File? mhc_ii_all_epitopes = pvacfuse.mhc_ii_all_epitopes
    File? mhc_ii_aggregated_report = pvacfuse.mhc_ii_aggregated_report
    File? mhc_ii_filtered_epitopes = pvacfuse.mhc_ii_filtered_epitopes
    File? mhc_ii_log = pvacfuse.mhc_ii_log

    # pVACfuse outputs - Combined
    File? combined_all_epitopes = pvacfuse.combined_all_epitopes
    File? combined_aggregated_report = pvacfuse.combined_aggregated_report
    File? combined_filtered_epitopes = pvacfuse.combined_filtered_epitopes

    # pVACfuse globbed outputs for all files
    Array[File] mhc_i_all = pvacfuse.mhc_i
    Array[File] mhc_ii_all = pvacfuse.mhc_ii
    Array[File] combined_all = pvacfuse.combined
  }
}
