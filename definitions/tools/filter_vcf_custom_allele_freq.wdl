version 1.0

task filterVcfCustomAlleleFreq {
  input {
    File vcf
    Float maximum_population_allele_frequency
    String field_name
  }

  Int space_needed_gb = 10 + round(size(vcf, "GB")*2)
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "mgibio/vep_helper-cwl:vep_113.3_v1"
    memory: "4GB"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String outfile = "annotated.af_filtered.vcf"
  String intermediate_file = "soft_filtered_af_intermediate.vcf"
  command <<<
    /usr/bin/perl /usr/bin/vcf_check.pl ~{vcf} ~{intermediate_file} \
    /usr/bin/perl /opt/vep/src/ensembl-vep/filter_vep --format vcf -o ~{intermediate_file} -i ~{vcf} \
    --filter "~{field_name} < ~{maximum_population_allele_frequency} or not ~{field_name}" \
    --soft_filter

    # Added step: Parse "PASS;filter_vep_pass" to "PASS", "PASS;filter_vep_fail" to "gnomade_af_0.001"
    awk 'BEGIN {FS=OFS="\t"} 
        /^#/ {print $0} 
        !/^#/ {
            if ($7 == "PASS;filter_vep_pass") $7 = "PASS";
            else if ($7 == "PASS;filter_vep_fail") $7 = "gnomade_af_0.001";
            print $0
        }' ~{intermediate_file} > ~{outfile}
  >>>

  output {
    File filtered_vcf = outfile
  }
}

workflow wf {
  input {
    File vcf
    Float maximum_population_allele_frequency
    String field_name
  }

  call filterVcfCustomAlleleFreq {
    input:
    vcf=vcf,
    maximum_population_allele_frequency=maximum_population_allele_frequency,
    field_name=field_name
  }
}
