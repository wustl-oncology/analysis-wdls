version 1.0

# Build the --phased-proximal-variants-vcf (-p) input for pVACseq.
#
# pVACseq wants ONE vcf containing both germline and somatic calls, under the
# TUMOR sample name, with phase set (PS) tags. In the short-read arm that took
# five steps (subworkflows/phase_vcf.wdl: reheader -> SelectVariants ->
# GATK3 CombineVariants -> SortVcf -> ReadBackedPhasing).
#
# Here HiPhase has already done the phasing, and it phased both vcfs against the
# same bam in a single invocation, so they already share PS ids. All that is
# left is bookkeeping: subset the somatic vcf to the tumor column, rename the
# germline sample to match, concatenate, sort.
#
# Note this vcf is NOT VEP annotated, and does not need to be -- pVACseq reads
# only genotypes and phase sets from it. (The short-read arm's equivalent file
# is likewise not re-annotated after the merge.)
task combinePhasedProximalVcf {
  input {
    File germline_vcf
    File germline_vcf_tbi
    File somatic_vcf
    File somatic_vcf_tbi
    String tumor_sample_name
    String output_basename = "phased_proximal_variants"
    Int cores = 4
  }

  Int space_needed_gb = 10 + round(6*size([germline_vcf, somatic_vcf], "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"
    memory: "8GB"
    cpu: cores
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    set -euxo pipefail
    bcftools --version

    # 1. Somatic vcf may carry both tumor and normal columns (DeepSomatic emits
    #    both). Keep only the tumor column so the two files have identical
    #    sample sets before concat.
    if bcftools query -l ~{somatic_vcf} | grep -qx "~{tumor_sample_name}"; then
      bcftools view -s "~{tumor_sample_name}" -Oz -o somatic.tumor.vcf.gz ~{somatic_vcf}
    else
      echo "WARNING: '~{tumor_sample_name}' not found in somatic vcf; using first sample" >&2
      first=$(bcftools query -l ~{somatic_vcf} | head -n1)
      bcftools view -s "${first}" -Oz -o somatic.tumor.vcf.gz ~{somatic_vcf}
      printf "%s\t%s\n" "${first}" "~{tumor_sample_name}" > somatic_rename.txt
      bcftools reheader -s somatic_rename.txt -o somatic.tumor.renamed.vcf.gz somatic.tumor.vcf.gz
      mv somatic.tumor.renamed.vcf.gz somatic.tumor.vcf.gz
    fi
    bcftools index -t somatic.tumor.vcf.gz

    # 2. Rename the germline sample to the tumor sample name. The germline calls
    #    came from the tumor bam (Clair3 on tumor), so this is a relabel, not a
    #    claim about provenance -- pVACseq just needs the columns to line up.
    old=$(bcftools query -l ~{germline_vcf} | head -n1)
    printf "%s\t%s\n" "${old}" "~{tumor_sample_name}" > germline_rename.txt
    bcftools reheader -s germline_rename.txt -o germline.renamed.vcf.gz ~{germline_vcf}
    bcftools index -t germline.renamed.vcf.gz

    # 3. Concatenate and sort. -a allows overlapping records; -D drops exact
    #    duplicates in case a site was called by both Clair3 and DeepSomatic.
    bcftools concat -a -D -Oz -o combined.unsorted.vcf.gz \
      germline.renamed.vcf.gz somatic.tumor.vcf.gz

    bcftools sort -Oz -o ~{output_basename}.vcf.gz combined.unsorted.vcf.gz
    bcftools index -t ~{output_basename}.vcf.gz

    # Sanity check: how many records carry a phase set?
    echo "records with PS tag:"
    bcftools query -f '[%PS\n]' ~{output_basename}.vcf.gz 2>/dev/null | grep -cv '^\.$' || true
  >>>

  output {
    File phased_vcf = output_basename + ".vcf.gz"
    File phased_vcf_tbi = output_basename + ".vcf.gz.tbi"
  }
}

workflow wf {
  input {
    File germline_vcf
    File germline_vcf_tbi
    File somatic_vcf
    File somatic_vcf_tbi
    String tumor_sample_name
  }
  call combinePhasedProximalVcf {
    input:
    germline_vcf=germline_vcf,
    germline_vcf_tbi=germline_vcf_tbi,
    somatic_vcf=somatic_vcf,
    somatic_vcf_tbi=somatic_vcf_tbi,
    tumor_sample_name=tumor_sample_name
  }
}
