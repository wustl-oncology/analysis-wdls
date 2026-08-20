version 1.0

# Clair3 germline small-variant calling on HiFi reads.

task clair3 {
  input {
    File bam
    File bam_bai
    File reference
    File reference_fai
    String sample_name
    # Model shipped inside the container. hifi_revio for Revio, hifi_sequel2 for
    # older instruments.
    String clair_model = "hifi_revio"
    Int cores = 16
  }

  Int space_needed_gb = 20 + round(2*size(bam, "GB") + size(reference, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "hkubal/clair3@sha256:1430f7b520674bb282b5a1ba165e8555b012df67b9cdb22ff3d48f26cf1474a3"
    memory: "~{cores * 4}GB"
    cpu: cores
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    set -euxo pipefail
    /opt/bin/run_clair3.sh --version
    mkdir -p clair3_out

    /opt/bin/run_clair3.sh \
      --bam_fn=~{bam} \
      --ref_fn=~{reference} \
      --threads=~{cores} \
      --platform="hifi" \
      --model_path="/opt/models/~{clair_model}" \
      --output=clair3_out \
      --sample_name=~{sample_name}

    # Drop LowQual; HiPhase should only see confident heterozygous sites.
    gunzip -c clair3_out/merge_output.vcf.gz \
      | awk '{if($7!="LowQual"){print $0}}' OFS=$'\t' \
      | bgzip -c > ~{sample_name}.clair3.vcf.gz
    tabix -p vcf ~{sample_name}.clair3.vcf.gz
  >>>

  output {
    File vcf = sample_name + ".clair3.vcf.gz"
    File vcf_tbi = sample_name + ".clair3.vcf.gz.tbi"
  }
}

workflow wf {
  input {
    File bam
    File bam_bai
    File reference
    File reference_fai
    String sample_name
    String? clair_model
  }
  call clair3 {
    input:
    bam=bam,
    bam_bai=bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    sample_name=sample_name,
    clair_model=select_first([clair_model, "hifi_revio"])
  }
}
