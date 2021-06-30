version 1.0

task varscanSomatic {
  input {
    File reference
    File reference_fai
    File reference_dict

    File tumor_bam
    File tumor_bam_bai

    File normal_bam
    File normal_bam_bai

    Int strand_filter = 0
    Int min_coverage = 8
    Float min_var_freq = 0.1
    Float p_value = 0.99
    File? roi_bed
  }

  Int space_needed_gb = 10 + round(size([reference, reference_fai, reference_dict, tumor_bam, tumor_bam_bai, normal_bam, normal_bam_bai], "GB"))
  runtime {
    memory: "12GB"
    cpu: 2
    docker: "mgibio/cle:v1.3.1"
    disks: "local-disk ~{space_needed_gb} HDD"
  }


  String outdir = "/cromwell_root"
  command <<<
    set -o errexit
    set -o nounset

    java -jar /opt/varscan/VarScan.jar somatic \
    <(/opt/samtools/bin/samtools mpileup --no-baq ~{if defined(roi_bed) then "-l ~{roi_bed}" else ""} -f "~{reference}" "~{normal_bam}" "~{tumor_bam}") \
    "~{outdir}/output" \
    --strand-filter "~{strand_filter}" \
    --min-coverage "~{min_coverage}" \
    --min-var-freq "~{min_var_freq}" \
    --p-value "~{p_value}" \
    --mpileup 1 \
    --output-vcf
  >>>

  output {
    File snvs = "~{outdir}/output.snp.vcf"
    File indels = "~{outdir}/output.indel.vcf"
  }
}

workflow wf {
  input {
    File reference
    File reference_fai
    File reference_dict

    File tumor_bam
    File tumor_bam_bai

    File normal_bam
    File normal_bam_bai

    Int strand_filter = 0
    Int min_coverage = 8
    Float min_var_freq = 0.1
    Float p_value = 0.99
    File? roi_bed
  }

  call varscanSomatic {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    tumor_bam=tumor_bam,
    tumor_bam_bai=tumor_bam_bai,
    normal_bam=normal_bam,
    normal_bam_bai=normal_bam_bai,
    strand_filter=strand_filter,
    min_coverage=min_coverage,
    min_var_freq=min_var_freq,
    p_value=p_value,
    roi_bed=roi_bed
  }
}
