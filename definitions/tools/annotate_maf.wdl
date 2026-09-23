version 1.0

task annotateMaf {
  input {
    File maf
    # A single-sample facetsSuite .rds (purity or hisens run) from runFacets
    File facets_rds
    String facets_algorithm = "em"  # enum [em, cncf]
  }

  Int space_needed_gb = 10 + round(3*size(maf, "GB") + size(facets_rds, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "mskcc/facets-suite:2.0.9"
    memory: "32GB"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String outfile = basename(maf, ".maf") + ".ccf.maf"

  command <<<
    set -eo pipefail

    # ccf_annotate_maf() joins the MAF to the FACETS segments with foverlaps() on
    # Chromosome/Start_Position/End_Position, and the segment chromosomes are unprefixed
    # (1..22, X). Strip "chr" from the MAF's Chromosome column (field 5) so the join hits.
    # Line 1 is vcf2maf's "#version" comment and line 2 is the header, hence NR>2.
    awk -F'\t' -v OFS='\t' 'NR>2 {gsub(/^chr/, "", $5)} 1' ~{maf} > no_chr.maf

    # The wrapper requires exactly one distinct Tumor_Sample_Barcode when given a
    # single-sample --facets-output, which vcfToMaf's --tumor-id guarantees.
    Rscript /usr/bin/facets-suite/annotate-maf-wrapper.R \
    --maf-file no_chr.maf \
    --facets-output ~{facets_rds} \
    --facets-algorithm ~{facets_algorithm} \
    --output ~{outfile}
  >>>

  output {
    File ccf_maf = outfile
  }
}

workflow wf {
  input {
    File maf
    File facets_rds
    String facets_algorithm = "em"
  }
  call annotateMaf {
    input:
    maf=maf,
    facets_rds=facets_rds,
    facets_algorithm=facets_algorithm
  }
}
