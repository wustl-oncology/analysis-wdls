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
    --filter "~{field_name} < ~{maximum_population_allele_frequency} or ~{field_name} = '.' or ~{field_name} = ''" \
    --soft_filter

    # Compute a dynamic FILTER ID that reflects the threshold, e.g., gnomade_af_0.001
    FILTER_ID=$(awk -v v="~{maximum_population_allele_frequency}" 'BEGIN{printf "gnomade_af_%g", v+0}')

    # Process the VCF to add custom filter definition in header and update filter names (Cromwell processes WDL files differently than bash scripts, and multi-line awk commands can cause syntax parsing issues)
    awk -v filter_id="$FILTER_ID" -v thresh="~{maximum_population_allele_frequency}" 'BEGIN {FS=OFS="\t"; header_processed=0} /^##FILTER=/ {print $0; next} /^#CHROM/ {if (!header_processed) {print "##FILTER=<ID=" filter_id ",Description=\"Variant failed gnomAD allele frequency filter (AF >= " thresh ")\">"; header_processed = 1} print $0; next} /^#/ {print $0; next} !/^#/ {if ($7 ~ /filter_vep_pass/) {gsub(/;filter_vep_pass/, "", $7); gsub(/^filter_vep_pass;/, "", $7); gsub(/^filter_vep_pass$/, "PASS", $7); if ($7 == "") $7 = "PASS"} else if ($7 ~ /filter_vep_fail/) {gsub(/;filter_vep_fail/, "", $7); gsub(/^filter_vep_fail;/, "", $7); gsub(/^filter_vep_fail$/, filter_id, $7); if ($7 == "") $7 = filter_id; if ($7 !~ filter_id) {if ($7 == "PASS") $7 = filter_id; else $7 = $7 ";" filter_id}} print $0}' ~{intermediate_file} > ~{outfile}
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
