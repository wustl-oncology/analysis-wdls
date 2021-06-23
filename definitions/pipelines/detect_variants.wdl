version 1.0

import "subworkflows/pindel.wdl" as p
import "subworkflows/mutect.wdl" as m
import "subworkflows/strelka_and_post_processing.wdl" as sapp
import "subworkflows/varscan_pre_and_post_processing.wdl" as vpapp
import "tools/combine_variants.wdl" as cv
import "tools/docm_add_variants.wdl" as dav
import "tools/vt_decompose.wdl" as vd

workflow DetectVariants {
  input {
    File reference
    File reference_fai
    File reference_dict

    String tumor_sample_name
    File tumor_bam
    File tumor_bam_bai

    String normal_sample_name
    File normal_bam
    File normal_bam_bai

    File roi_intervals
    Int scatter_count = 50

    Boolean strelka_exome_mode
    Int strelka_cpu_reserved = 8

    Int varscan_strand_filter = 0
    Int varscan_min_coverage = 8
    Float varscan_min_var_freq = 0.1
    Float varscan_p_value = 0.99
    Float? varscan_max_normal_freq

    Boolean filter_docm_variants = true
  }

  call m.mutect {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    tumor_sample_name=tumor_sample_name,
    tumor_bam=tumor_bam,
    tumor_bam_bai=tumor_bam_bai,
    normal_bam=normal_bam,
    normal_bam_bai=normal_bam_bai,
    interval_list=roi_intervals,
    scatter_count=scatter_count
  }

  call sapp.strelkaAndPostProcessing as strelka {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    tumor_bam=tumor_bam,
    tumor_bam_bai=tumor_bam_bai,
    normal_bam=normal_bam,
    normal_bam_bai=normal_bam_bai,
    interval_list=roi_intervals,
    exome_mode=strelka_exome_mode,
    cpu_reserved=strelka_cpu_reserved,
    normal_sample_name=normal_sample_name,
    tumor_sample_name=normal_sample_name
  }

  call vpapp.varscanPreAndPostProcessing as varscan {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    tumor_bam=tumor_bam,
    tumor_bam_bai=tumor_bam_bai,
    normal_bam=normal_bam,
    normal_bam_bai=normal_bam_bai,
    interval_list=roi_intervals,
    scatter_count=scatter_count,
    strand_filter=varscan_strand_filter,
    min_coverage=varscan_min_coverage,
    min_var_freq=varscan_min_var_freq,
    p_value=varscan_p_value,
    max_normal_freq=varscan_max_normal_freq,
    normal_sample_name=normal_sample_name,
    tumor_sample_name=tumor_sample_name
  }

  call p.pindel {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    tumor_bam=tumor_bam,
    tumor_bam_bai=tumor_bam_bai,
    normal_bam=normal_bam,
    normal_bam_bai=normal_bam_bai,
    interval_list=roi_intervals,
    scatter_count=scatter_count,
    insert_size=pindel_insert_size,
    tumor_sample_name=tumor_sample_name,
    normal_sample_name=normal_sample_name,
  }

  call dc.docmCle as docm {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    tumor_bam=tumor_bam,
    tumor_bam_bai=tumor_bam_bai,
    normal_bam=normal_bam,
    normal_bam_bai=normal_bam_bai,
    docm_vcf=docm_vcf,
    docm_vcf_tbi=docm_vcf_tbi,
    interval_list=roi_intervals,
    filter_docm_variants=filter_docm_variants
  }

  call cv.combineVariants as combine {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,

    mutect_vcf=mutect.filtered_vcf,
    mutect_vcf_tbi=mutect.filtered_vcf_tbi,

    strelka_vcf=strelka.filtered_vcf,
    strelka_vcf_tbi=strelka.filtered_vcf_tbi,

    varscan_vcf=varscan.filtered_vcf,
    varscan_vcf_tbi=varscan.filtered_vcf_tbi,

    pindel_vcf=pindel.filtered_vcf,
    pindel_vcf_tbi=pindel.filtered_vcf_tbi
  }

  call dav.docmAddVariants as addDocmVariants {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    docm_vcf=docm.docm_variants_vcf,
    docm_vcf_tbi=docm.docm_variants_vcf_tbi,
    callers_vcf=combine.combined_vcf,
    callers_vcf_tbi=combine.combined_vcf_tbi
  }

  call vd.vtDecompose as decompose {
    input:
    vcf=addDocmVariants.merged_vcf,
    vcf_tbi=addDocmVariants.merged_vcf_tbi
  }

  call iv.indexVcf as decomposeIndex {
    input: vcf=decompose.decomposed_vcf
  }

  # ---- TODO: I am here
}
