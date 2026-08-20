version 1.0

# Goal: Normalize/Process vcf , to make them suitable for downstream analysis (eg. by pVACseq)
# Decomposes multiple-allelic records (variants with multiple alternate alleles) into single-allele records
# Left align idels
# Fix formating issue

# Steps:
# Approach: read input VCF file (bcftools view)
# fix AD field format (with sed) (change header from AD,Number=. to AD, Number=R)
# normalize variants (bcftools norm) (split multip-allelic sites into separate records)(and left-align indels using the reference genome)
# sort and index


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
