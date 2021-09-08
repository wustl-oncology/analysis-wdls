version 1.0

task docmGatkHaplotypeCaller {
  input {
    File reference
    File reference_fai
    File reference_dict
    File normal_bam
    File normal_bam_bai
    File bam
    File bam_bai
    File docm_vcf
    File docm_vcf_tbi
    File interval_list
  }

  Float copied_size = size([docm_vcf, interval_list], "GB")
  Int space_needed_gb = 10 + round(copied_size*3 + size([reference, reference_fai, reference_dict, normal_bam, normal_bam_bai, bam, bam_bai, docm_vcf_tbi], "GB"))
  runtime {
    memory: "9GB"
    docker: "broadinstitute/gatk:4.1.2.0"
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  command <<<
    ~{if defined(normal_bam) then "mv ~{normal_bam} ~{basename(normal_bam)}" else ""}
    ~{if defined(normal_bam_bai) then "mv ~{normal_bam_bai} ~{basename(normal_bam_bai)}" else ""}
    mv ~{bam} ~{basename(bam)}; mv ~{bam_bai} ~{basename(bam_bai)}


    # Extracting the header from the interval_list
    cat "~{interval_list}" | grep '^@' > docm.interval_list
    # Extracting the docm regions with a 100bp flanking region on both directions
    zcat "~{docm_vcf}" | grep ^chr | awk '{FS = "\t";OFS = "\t";print $1,$2-100,$2+100,"+",$1"_"$2-100"_"$2+100}' >> docm.interval_list

    /gatk/gatk HaplotypeCaller --java-options "-Xmx8g" -R "~{reference}" \
    -I "~{basename(bam)}" ~{if defined(normal_bam) then "-I ~{basename(normal_bam)}" else ""} \
    --alleles "~{docm_vcf}" \
    -L docm.interval_list \
    --genotyping-mode GENOTYPE_GIVEN_ALLELES \
    -O docm_raw_variants.vcf
  >>>

  output {
    File docm_raw_variants = "docm_raw_variants.vcf"
  }
}

workflow wf {
  input {
    File reference
    File reference_fai
    File reference_dict
    File normal_bam
    File normal_bam_bai
    File bam
    File bam_bai
    File docm_vcf
    File docm_vcf_tbi
    File interval_list
  }
  call docmGatkHaplotypeCaller {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    normal_bam=normal_bam,
    normal_bam_bai=normal_bam_bai,
    bam=bam,
    bam_bai=bam_bai,
    docm_vcf=docm_vcf,
    docm_vcf_tbi=docm_vcf_tbi,
    interval_list=interval_list
  }
}
