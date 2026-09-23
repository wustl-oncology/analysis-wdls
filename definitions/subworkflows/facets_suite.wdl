version 1.0

import "../tools/snp_pileup.wdl" as sp
import "../tools/facets.wdl" as f
import "../tools/vcf_to_maf.wdl" as vtm
import "../tools/annotate_maf.wdl" as am

# Estimate tumor purity/ploidy and allele-specific copy number with FACETS, then annotate
# the somatic variants with local copy number and cancer-cell fraction (CCF).
#
# The MAF branch needs nothing from the VCF but Chromosome/Start_Position/End_Position and
# t_alt_count/t_ref_count -- ccf_annotate_maf() reads no other column and passes the rest
# through untouched. t_alt_count/t_ref_count come from the tumor genotype's AD field, which
# the callers do not write uniformly (strelka emits AU/CU/GU/TU and TAR/TIR instead), so
# `variants_vcf` must be a VCF that has already been through the bam-readcount annotation
# in detect_variants.wdl. `final_filtered_vcf` is the natural choice.
workflow facetsSuite {
  input {
    File tumor_bam
    File tumor_bam_bai
    File normal_bam
    File normal_bam_bai

    # Common SNP positions to pile up at, e.g. dbSNP
    File snps_vcf
    File snps_vcf_tbi

    File variants_vcf
    String tumor_sample_name
    String normal_sample_name

    File reference
    File reference_fai

    String genome = "hg38"  # enum [hg18, hg19, hg38]
    String ncbi_build = "GRCh38"

    Int purity_cval = 1000
    Int cval = 500
    Int seed = 918
    Int? min_nhet
    Int? purity_min_nhet
    Int? snp_window_size
    Int? normal_depth
    Int? pseudo_snps
    Int? max_depth
    Boolean facets_everything = false

    String facets_algorithm = "em"  # enum [em, cncf]
    # Annotate CCF from the high-sensitivity run by default; set true to use the purity run.
    Boolean annotate_with_purity_run = false
  }

  call sp.snpPileup {
    input:
    normal_bam=normal_bam,
    normal_bam_bai=normal_bam_bai,
    tumor_bam=tumor_bam,
    tumor_bam_bai=tumor_bam_bai,
    snps_vcf=snps_vcf,
    snps_vcf_tbi=snps_vcf_tbi,
    pseudo_snps=pseudo_snps,
    max_depth=max_depth
  }

  call f.runFacets {
    input:
    counts_file=snpPileup.snp_pileup_counts,
    sample_name=tumor_sample_name,
    genome=genome,
    purity_cval=purity_cval,
    cval=cval,
    seed=seed,
    min_nhet=min_nhet,
    purity_min_nhet=purity_min_nhet,
    snp_window_size=snp_window_size,
    normal_depth=normal_depth,
    everything=facets_everything
  }

  call vtm.vcfToMaf {
    input:
    vcf=variants_vcf,
    tumor_sample_name=tumor_sample_name,
    normal_sample_name=normal_sample_name,
    reference=reference,
    reference_fai=reference_fai,
    ncbi_build=ncbi_build
  }

  File rds_to_annotate = if annotate_with_purity_run then runFacets.purity_rds else runFacets.hisens_rds

  call am.annotateMaf {
    input:
    maf=vcfToMaf.maf,
    facets_rds=rds_to_annotate,
    facets_algorithm=facets_algorithm
  }

  output {
    File facets_purity_rds = runFacets.purity_rds
    File facets_purity_seg = runFacets.purity_seg
    File facets_purity_png = runFacets.purity_png
    File facets_hisens_rds = runFacets.hisens_rds
    File facets_hisens_seg = runFacets.hisens_seg
    File facets_hisens_png = runFacets.hisens_png
    File facets_run_details = runFacets.run_details
    File? facets_qc = runFacets.qc
    File? facets_gene_level = runFacets.gene_level
    File? facets_arm_level = runFacets.arm_level
    File facets_ccf_maf = annotateMaf.ccf_maf
  }
}
