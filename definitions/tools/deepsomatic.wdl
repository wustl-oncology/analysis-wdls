version 1.0

# DeepSomatic: the somatic SNV/indel caller for PacBio HiFi tumor/normal pairs.
# Replaces the Mutect2 + Strelka2 + VarScan2 + DoCM ensemble used in the
# short-read arm -- there is no equivalent consensus caller for HiFi, and
# DeepSomatic's PACBIO model is trained on exactly this data.

task deepSomatic {
  input {
    File tumor_bam
    File tumor_bam_bai
    File normal_bam
    File normal_bam_bai
    File reference
    File reference_fai
    File? regions_bed
    # These become the vcf column headers, and downstream
    # tools/combine_phased_proximal_vcf.wdl looks the tumor one up by name.
    # Keep them identical to the workflow-level sample names.
    String tumor_sample_name
    String normal_sample_name
    String output_basename
    Int cores = 16
  }

  Int space_needed_gb = 30 + round(2*size([tumor_bam, normal_bam], "GB") + size(reference, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "google/deepsomatic@sha256:d9797b8950bf615ec7010d1336b7ee0a2f12ea09323dc3585f7e9fe39b082bde"
    memory: "~{cores * 8}GB"
    cpu: cores
    bootDiskSizeGb: 30
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    set -euxo pipefail
    mkdir -p deepsomatic_out/logs
    /opt/deepvariant/bin/deepsomatic/run_deepsomatic --version

    /opt/deepvariant/bin/deepsomatic/run_deepsomatic \
      --model_type=PACBIO \
      --ref=~{reference} \
      --reads_normal=~{normal_bam} \
      --reads_tumor=~{tumor_bam} \
      --output_vcf=deepsomatic_out/~{output_basename}.vcf.gz \
      --output_gvcf=deepsomatic_out/~{output_basename}.g.vcf.gz \
      --sample_name_tumor=~{tumor_sample_name} \
      --sample_name_normal=~{normal_sample_name} \
      --num_shards=~{cores} \
      --postprocess_variants_extra_args="--cpus=~{cores / 2},--num_partitions=~{cores / 2}" \
      --logging_dir=deepsomatic_out/logs \
      ~{"--regions=" + regions_bed}

    # pVACseq is run with --pass-only, so drop non-PASS here and keep the
    # downstream files small.
    bcftools view -f PASS -Oz \
      -o deepsomatic_out/~{output_basename}.PASS.vcf.gz \
      deepsomatic_out/~{output_basename}.vcf.gz
    tabix -p vcf deepsomatic_out/~{output_basename}.PASS.vcf.gz
  >>>

  output {
    File vcf = "deepsomatic_out/" + output_basename + ".PASS.vcf.gz"
    File vcf_tbi = "deepsomatic_out/" + output_basename + ".PASS.vcf.gz.tbi"
    File gvcf = "deepsomatic_out/" + output_basename + ".g.vcf.gz"
  }
}

task gatherDeepSomatic {
  input {
    Array[File] vcfs
    String output_basename
    Int cores = 8
  }

  Int space_needed_gb = 10 + round(4*size(vcfs, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"
    memory: "~{cores * 4}GB"
    cpu: cores
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    set -euxo pipefail
    bcftools --version

    # Re-localize and re-index defensively; scattered indices sometimes arrive
    # with a timestamp older than their vcf, which bcftools rejects.
    count=0
    for file in ~{sep=' ' vcfs}; do
      cp "$file" "${count}.vcf.gz"
      tabix -p vcf "${count}.vcf.gz"
      count=$((count+1))
    done
    ls -1 ./*.vcf.gz > vcf.list

    # DeepSomatic emits somatic calls as 1/1. HiPhase will not phase a
    # homozygous genotype, so force everything to 0/1 before phasing.
    bcftools concat -a -f vcf.list \
      | bcftools +setGT -- -t q -n c:"0/1" -i 'FMT/GT="1/1"' \
      | bcftools sort -Oz -o ~{output_basename}.vcf.gz
    tabix -p vcf ~{output_basename}.vcf.gz
  >>>

  output {
    File vcf = output_basename + ".vcf.gz"
    File vcf_tbi = output_basename + ".vcf.gz.tbi"
  }
}

workflow wf {
  input {
    File tumor_bam
    File tumor_bam_bai
    File normal_bam
    File normal_bam_bai
    File reference
    File reference_fai
    Array[File] contigs
    String tumor_sample_name
    String normal_sample_name
    String output_basename = "somatic"
  }

  scatter (ctg in contigs) {
    call deepSomatic {
      input:
      tumor_bam=tumor_bam,
      tumor_bam_bai=tumor_bam_bai,
      normal_bam=normal_bam,
      normal_bam_bai=normal_bam_bai,
      reference=reference,
      reference_fai=reference_fai,
      regions_bed=ctg,
      tumor_sample_name=tumor_sample_name,
      normal_sample_name=normal_sample_name,
      output_basename=output_basename + "." + basename(ctg, ".bed")
    }
  }

  call gatherDeepSomatic {
    input: vcfs=deepSomatic.vcf, output_basename=output_basename
  }

  output {
    File vcf = gatherDeepSomatic.vcf
    File vcf_tbi = gatherDeepSomatic.vcf_tbi
  }
}
