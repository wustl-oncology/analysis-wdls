version 1.0

# Decompose multi-allelic records and left-align indels.
#
# This is the HiFi equivalent of tools/vt_decompose.wdl in the short-read arm.
# pVACseq requires one ALT allele per record; a multi-allelic site will either
# be skipped or mis-annotated.
#
# DeepSomatic writes `##FORMAT=<ID=AD,Number=.>`, but
# with Number=. bcftools does not know AD is per-allele. Using Number=R
# ensures it subsets correctly (instead of copying the whole AD array to both output records). Downstream this is what pVACseq reads for --tdna-cov.
task bcftoolsNorm {
  input {
    File input_vcf
    File reference
    File reference_fai
    Int cores = 4
  }

  String outbase = basename(basename(input_vcf, ".gz"), ".vcf")
  Int space_needed_gb = 10 + round(4*size(input_vcf, "GB") + size(reference, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"
    memory: "16GB"
    cpu: cores
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    set -euxo pipefail
    bcftools --version

    bcftools view ~{input_vcf} \
      | sed -e 's/ID=AD,Number=\./ID=AD,Number=R/' \
      | bcftools norm \
          --threads ~{cores - 1} \
          --multiallelics - \
          --fasta-ref ~{reference} \
          --output-type u \
      | bcftools sort -Oz -o ~{outbase}.norm.vcf.gz

    bcftools index --threads ~{cores - 1} -t ~{outbase}.norm.vcf.gz
  >>>

  output {
    File normalized_vcf = outbase + ".norm.vcf.gz"
    File normalized_vcf_tbi = outbase + ".norm.vcf.gz.tbi"
  }
}

workflow wf {
  input {
    File input_vcf
    File reference
    File reference_fai
  }
  call bcftoolsNorm {
    input: input_vcf=input_vcf, reference=reference, reference_fai=reference_fai
  }
}
