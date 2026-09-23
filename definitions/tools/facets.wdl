version 1.0

task runFacets {
  input {
    File counts_file  # <sample>.snp_pileup.gz from snpPileup
    String sample_name
    String genome = "hg38"  # enum [hg18, hg19, hg38]
    # Two-pass run: purity_cval drives the first (coarse) pass, cval the second (hisens) pass.
    Int purity_cval = 1000
    Int cval = 500
    Int seed = 918
    Int? min_nhet
    Int? purity_min_nhet
    Int? snp_window_size
    Int? normal_depth
    Float? diplogr
    # --everything additionally emits gene-level, arm-level and fit-QC tables
    Boolean everything = false
    String facets_lib_path = "/usr/local/lib/R/site-library"
  }

  Int space_needed_gb = 10 + round(4*size(counts_file, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "mskcc/facets-suite:2.0.9"
    memory: "32GB"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    set -eo pipefail

    # The 2.0.9 wrapper takes --directory as a plain path and writes flat into it, so "."
    # puts every output in the task working directory. --facets-lib-path is required.
    Rscript /usr/bin/facets-suite/run-facets-wrapper.R \
    --counts-file ~{counts_file} \
    --sample-id ~{sample_name} \
    --directory . \
    --facets-lib-path ~{facets_lib_path} \
    --genome ~{genome} \
    --purity-cval ~{purity_cval} \
    --cval ~{cval} \
    --seed ~{seed} \
    ~{if everything then "--everything" else ""} \
    ~{if defined(min_nhet) then "--min-nhet ~{min_nhet}" else ""} \
    ~{if defined(purity_min_nhet) then "--purity-min-nhet ~{purity_min_nhet}" else ""} \
    ~{if defined(snp_window_size) then "--snp-window-size ~{snp_window_size}" else ""} \
    ~{if defined(normal_depth) then "--normal-depth ~{normal_depth}" else ""} \
    ~{if defined(diplogr) then "--dipLogR ~{diplogr}" else ""}
  >>>

  # Output names are fixed by the wrapper's two-pass branch: <sample>_{purity,hisens}.{rds,seg,png}
  # plus a <sample>.txt run-details table. Verified against mskcc/facets-suite:2.0.9 -- note that
  # later versions rename the .seg files to _diplogR.{un,}adjusted.seg, so this list is tied to the
  # pinned tag above.
  output {
    File purity_rds = sample_name + "_purity.rds"
    File purity_seg = sample_name + "_purity.seg"
    File purity_png = sample_name + "_purity.png"
    File hisens_rds = sample_name + "_hisens.rds"
    File hisens_seg = sample_name + "_hisens.seg"
    File hisens_png = sample_name + "_hisens.png"
    File run_details = sample_name + ".txt"
    File? qc = sample_name + ".qc.txt"
    File? gene_level = sample_name + ".gene_level.txt"
    File? arm_level = sample_name + ".arm_level.txt"
  }
}

workflow wf {
  input {
    File counts_file
    String sample_name
    String genome = "hg38"
    Int purity_cval = 1000
    Int cval = 500
    Int seed = 918
    Boolean everything = false
  }
  call runFacets {
    input:
    counts_file=counts_file,
    sample_name=sample_name,
    genome=genome,
    purity_cval=purity_cval,
    cval=cval,
    seed=seed,
    everything=everything
  }
}
