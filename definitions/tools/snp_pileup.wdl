version 1.0

task snpPileup {
  input {
    File normal_bam
    File normal_bam_bai
    File tumor_bam
    File tumor_bam_bai
    # VCF of common SNP positions to pile up at, e.g. dbSNP
    File snps_vcf
    File snps_vcf_tbi
    Int? pseudo_snps  # wrapper default 50
    Int? max_depth    # wrapper default 4000
  }

  Int space_needed_gb = 10 + round(2*size([normal_bam, normal_bam_bai, tumor_bam, tumor_bam_bai, snps_vcf, snps_vcf_tbi], "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "mskcc/facets-suite:2.0.9"
    memory: "32GB"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  # The wrapper appends ".snp_pileup.gz" to --output-prefix. The name is arbitrary: it is
  # only ever passed on as --counts-file, and cromwell gives each task its own directory.
  String output_prefix = "facets"

  command <<<
    set -eo pipefail

    # snp-pileup index-seeks both bams and the vcf, so each index has to sit beside its
    # file; cromwell may localize them into separate input directories. Skip the move when
    # the file is already in the working directory so this stays safe under `set -e`.
    for f in ~{normal_bam} ~{normal_bam_bai} ~{tumor_bam} ~{tumor_bam_bai} ~{snps_vcf} ~{snps_vcf_tbi}; do
      [ -e "$(basename "$f")" ] || mv "$f" .
    done

    # normal is passed before tumor: the wrapper emits them in that order, so the counts
    # file has the normal as File1 and the tumor as File2, which is what run_facets expects.
    Rscript /usr/bin/facets-suite/snp-pileup-wrapper.R \
    --snp-pileup-path /usr/bin/snp-pileup \
    --vcf-file ~{basename(snps_vcf)} \
    --normal-bam ~{basename(normal_bam)} \
    --tumor-bam ~{basename(tumor_bam)} \
    --output-prefix ~{output_prefix} \
    ~{if defined(pseudo_snps) then "--pseudo-snps ~{pseudo_snps}" else ""} \
    ~{if defined(max_depth) then "--max-depth ~{max_depth}" else ""}
  >>>

  output {
    File snp_pileup_counts = output_prefix + ".snp_pileup.gz"
  }
}

workflow wf {
  input {
    File normal_bam
    File normal_bam_bai
    File tumor_bam
    File tumor_bam_bai
    File snps_vcf
    File snps_vcf_tbi
    Int? pseudo_snps
    Int? max_depth
  }
  call snpPileup {
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
}
