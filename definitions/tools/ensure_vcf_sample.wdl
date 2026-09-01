version 1.0

# Guarantee a named sample column exists in a VCF, adding an all-missing one if
# it does not.
#
# DeepSomatic in tumor/normal mode does not reliably emit a normal genotype
# column, even though it takes --sample_name_normal. Both bam-readcount and
# vcf-readcount-annotator look the sample up by name and hard-fail when it is
# absent:
#
#     sample_index = vcf_file.samples.index(sample)   # ValueError
#
# This task ensures the normal column exists by adding it when necessary.
# If the sample is already present, the VCF passes through unchanged (just
# bgzipped and indexed), so it is safe to leave in the pipeline once
# DeepSomatic's behaviour is confirmed either way.
task ensureVcfSample {
  input {
    File vcf
    String sample_name
    String output_basename = "with_sample"
    Int cores = 4
  }

  Int space_needed_gb = 10 + round(6*size(vcf, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"
    memory: "8GB"
    cpu: cores
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    set -euo pipefail
    bcftools --version

    echo "samples currently in vcf:"
    bcftools query -l ~{vcf}

    if bcftools query -l ~{vcf} | grep -qx "~{sample_name}"; then
      echo "'~{sample_name}' already present -- passing through"
      bcftools view -Oz -o ~{output_basename}.vcf.gz ~{vcf}
      echo "false" > added.txt
    else
      echo "'~{sample_name}' absent -- appending an all-missing genotype column"
      # Build a per-record genotype string shaped like that record's FORMAT
      # field: "./." for GT, "." for everything else. vcf-readcount-annotator
      # overwrites these with real values downstream.
      bcftools view ~{vcf} \
        | awk -v s="~{sample_name}" 'BEGIN{OFS="\t"}
            /^##/      { print; next }
            /^#CHROM/  { print $0, s; next }
            {
              n = split($9, f, ":");
              g = "";
              for (i = 1; i <= n; i++) {
                v = (f[i] == "GT") ? "./." : ".";
                g = (i == 1) ? v : g ":" v;
              }
              print $0, g
            }' \
        | bcftools view -Oz -o ~{output_basename}.vcf.gz
      echo "true" > added.txt
    fi

    bcftools index --threads ~{cores - 1} -t ~{output_basename}.vcf.gz

    echo "samples after:"
    bcftools query -l ~{output_basename}.vcf.gz
  >>>

  output {
    File vcf_with_sample = output_basename + ".vcf.gz"
    File vcf_with_sample_tbi = output_basename + ".vcf.gz.tbi"
    # "true" if a column had to be created. Worth checking on run 1 -- it tells
    # you definitively whether DeepSomatic emitted a normal column.
    String sample_was_added = read_string("added.txt")
  }
}

workflow wf {
  input {
    File vcf
    String sample_name
  }
  call ensureVcfSample {
    input: vcf=vcf, sample_name=sample_name
  }
}
