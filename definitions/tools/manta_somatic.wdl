version 1.0

task mantaSomatic {
  input {
    File? normal_bam
    File? normal_bam_bai
    File tumor_bam
    File tumor_bam_bai
    File reference
    File reference_fai
    File reference_dict
    File? call_regions
    File? call_regions_tbi
    Boolean non_wgs = false
    Boolean output_contigs = false
  }

  Int cores = 12
  Float ref_size = size([reference, reference_fai, reference_dict], "GB")
  Float tumor_size = size([tumor_bam, tumor_bam_bai], "GB")
  Float normal_size = size([normal_bam, normal_bam_bai], "GB")
  Float regions_size = size([call_regions, call_regions_tbi], "GB")
  Int size_needed_gb = 10 + round(ref_size + tumor_size + normal_size + regions_size)
  runtime {
    docker: "mgibio/manta_somatic-cwl:1.6.0"
    cpu: cores
    ram: "24GB"
    bootDiskSizeGb: 10
    disks: "local-disk ~{size_needed_gb} HDD"
  }

  String outdir = "/cromwell_root"
  command <<<
    /usr/bin/python /usr/bin/manta/bin/configManta.py \
    ~{if non_wgs then "--exome" else ""} \
    ~{if output_contigs then "--outputContig" else ""} \
    ~{if defined(call_regions) then "--callRegions ~{call_regions}" else ""} \  # -5
    --referenceFasta ~{reference} \  # -4
    --tumorBam ~{tumor_bam} \  # -3
    ~{if defined(normal_bam) then "--normalBam ~{normal_bam}" else ""} \  # -2
    --runDir ~{outdir} \  # -1
    && /usr/bin/python runWorkflow.py -m local \  # 0
    -j ~{cores} \  # 1
  >>>

  output {
    File? diploid_variants = "results/variants/diploidSV.vcf.gz"
    File? diploid_variants_tbi = "results/variants/diploidSV.vcf.gz.tbi"

    File? somatic_variants = "results/variants/somaticSV.vcf.gz"
    File? somatic_variants_tbi = "results/variants/somaticSV.vcf.gz.tbi"

    File all_candidates = "results/variants/candidateSV.vcf.gz"
    File all_candidates_tbi = "results/variants/candidateSV.vcf.gz.tbi"

    File small_candidates = "results/variants/candidateSmallIndels.vcf.gz"
    File small_candidates_tbi = "results/variants/candidateSmallIndels.vcf.gz.tbi"

    File? tumor_only_variants = "results/variants/tumorSV.vcf.gz"
    File? tumor_only_variants_tbi = "results/variants/tumorSV.vcf.gz.tbi"
  }
}
