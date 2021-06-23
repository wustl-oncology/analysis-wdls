version 1.0

task strelka {
  input {
    File tumor_bam
    File tumor_bam_bai
    File normal_bam
    File normal_bam_bai
    File reference
    File reference_fai
    File reference_dict
    Boolean exome_mode
    Int? cpu_reserved
  }

  Int space_needed_gb = 10 + round(size([tumor_bam, tumor_bam_bai, normal_bam, normal_bam_bai, reference, reference_fai, reference_dict], "GB"))
  runtime {
    memory: "4GB"
    cpu: 4
    docker: "mgibio/strelka-cwl:2.9.9"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String outdir = "/cromwell_root/"  # TODO: output dir
  command <<<
    mv ~{tumor_bam} ~{basename(tumor_bam)}; mv ~{tumor_bam_bai} ~{basename(tumor_bam_bai)}
    mv ~{normal_bam} ~{basename(normal_bam)}; mv ~{normal_bam_bai} ~{basename(normal_bam_bai)}
    /usr/bin/perl /usr/bin/docker_helper.pl \
    ~{if defined(cpu_reserved) then cpu_reserved else ""} \
    ~{outdir} --tumorBam=~{basename(tumor_bam)} --normalBam=~{basename(normal_bam)} \
    --referenceFasta=~{reference} \
    ~{if exome_mode then "--exome" else ""}
  >>>

  output {
    File indels = "results/variants/somatic.indels.vcf.gz"
    File indels_tbi = "results/variants/somatic.indels.vcf.gz.tbi"
    File snvs = "results/variants/somatic.snvs.vcf.gz"
    File snvs_tbi = "results/variants/somatic.snvs.vcf.gz.tbi"
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
    File reference_dict
    Boolean exome_mode
    Int? cpu_reserved
  }

  call strelka {
    input:
    tumor_bam=tumor_bam,
    tumor_bam_bai=tumor_bam_bai,
    normal_bam=normal_bam,
    normal_bam_bai=normal_bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    exome_mode=exome_mode,
    cpu_reserved=cpu_reserved,
  }
}
