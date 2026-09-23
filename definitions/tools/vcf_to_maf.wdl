version 1.0

task vcfToMaf {
  input {
    File vcf
    String tumor_sample_name
    String normal_sample_name
    File reference
    File reference_fai
    String ncbi_build = "GRCh38"
  }

  Int space_needed_gb = 10 + round(3*size(vcf, "GB") + size([reference, reference_fai], "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    # TODO: replace with an org-owned image before merging. Every other tool in this repo
    # pulls from mgibio/, griffithlab/ or quay.io/biocontainers/; a personal Docker Hub
    # namespace is a single point of failure for a shared pipeline.
    docker: "myhoang04/vcf2maf:1.6.22"
    memory: "32GB"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String outfile = tumor_sample_name + ".maf"

  command <<<
    set -eo pipefail

    # vcf2maf reads a plain-text VCF, so decompress first if needed.
    case "~{basename(vcf)}" in
      *.gz) gunzip -c ~{vcf} > input.vcf ;;
      *)    cp ~{vcf} input.vcf ;;
    esac

    # --inhibit-vep: the input VCF already carries VEP CSQ annotations from the
    # detect_variants step, so re-running VEP here would be wasted work.
    perl /opt/vcf2maf.pl \
    --input-vcf input.vcf \
    --output-maf ~{outfile} \
    --tumor-id ~{tumor_sample_name} \
    --normal-id ~{normal_sample_name} \
    --vcf-tumor-id ~{tumor_sample_name} \
    --vcf-normal-id ~{normal_sample_name} \
    --ref-fasta ~{reference} \
    --inhibit-vep \
    --ncbi-build ~{ncbi_build}
  >>>

  output {
    File maf = outfile
  }
}

workflow wf {
  input {
    File vcf
    String tumor_sample_name
    String normal_sample_name
    File reference
    File reference_fai
    String ncbi_build = "GRCh38"
  }
  call vcfToMaf {
    input:
    vcf=vcf,
    tumor_sample_name=tumor_sample_name,
    normal_sample_name=normal_sample_name,
    reference=reference,
    reference_fai=reference_fai,
    ncbi_build=ncbi_build
  }
}
