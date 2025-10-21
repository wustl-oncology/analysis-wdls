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
  String step1_soft_filtered_file = "1_soft_filtered_af_intermediate.vcf"
  String step2_header_processed_file = "2_header_processed_intermediate.vcf"
  String step3_filter_vep_pass_processed_file = "3_filter_vep_pass_processed_intermediate.vcf"
  String step4_filter_vep_fail_processed_file = "4_filter_vep_fail_processed_intermediate.vcf"
  command <<<
    /usr/bin/perl /usr/bin/vcf_check.pl ~{vcf} ~{step1_soft_filtered_file} \
    /usr/bin/perl /opt/vep/src/ensembl-vep/filter_vep --format vcf -o ~{step1_soft_filtered_file} -i ~{vcf} \
    --filter "~{field_name} < ~{maximum_population_allele_frequency} or not ~{field_name}" \
    --soft_filter

    # Compute a dynamic FILTER ID that reflects the threshold, e.g., gnomade_af_0.001
    FILTER_ID=$(awk -v v="~{maximum_population_allele_frequency}" 'BEGIN{printf "gnomade_af_%g", v+0}')

    # Step 1: Remove unwanted filter headers and add custom filter definition
    echo "Step 1: Removing unwanted filter headers and adding custom filter definition..."
    awk -v filter_id="$FILTER_ID" -v thresh="~{maximum_population_allele_frequency}" '
    BEGIN {FS=OFS="\t"; header_processed=0} 
    /^##FILTER=<ID=filter_vep_pass/ {next}
    /^##FILTER=<ID=filter_vep_fail/ {next}
    /^##FILTER=/ {print $0; next} 
    /^#CHROM/ {
        if (!header_processed) {
            print "##FILTER=<ID=" filter_id ",Description=\"Variant failed gnomAD allele frequency filter (AF >= " thresh ")\">"
            header_processed = 1
        }
        print $0
        next
    }
    /^#/ {print $0; next} 
    !/^#/ {print $0}
    ' ~{step1_soft_filtered_file} > ~{step2_header_processed_file}

    # Step 2: Process filter_vep_pass entries
    # Handles three cases: 
      # "OTHER_FILTER;filter_vep_pass" -> "OTHER_FILTER"
      # "filter_vep_pass;OTHER_FILTER" -> "OTHER_FILTER"
      # "PASS;filter_vep_pass" -> "PASS"
    echo "Step 2: Processing filter_vep_pass entries..."
    awk '
    BEGIN {FS=OFS="\t"} 
    /^#/ {print $0; next} 
    !/^#/ {
        if ($7 ~ /filter_vep_pass/) {
            # Remove all occurrences of filter_vep_pass
            gsub(/filter_vep_pass/, "", $7)
            # Clean up any resulting double semicolons or leading/trailing semicolons
            gsub(/;;+/, ";", $7)
            gsub(/^;/, "", $7)
            gsub(/;$/, "", $7)
            # If the result is empty, set to PASS
            if ($7 == "") $7 = "PASS"
        }
        print $0
    }
    ' ~{step2_header_processed_file} > ~{step3_filter_vep_pass_processed_file}

    # Step 3: Process filter_vep_fail entries and add custom filter
    # Handles four cases: 
      # "OTHER_FILTER;filter_vep_fail" -> "OTHER_FILTER;gnomade_af_0.001"
      # "filter_vep_fail;OTHER_FILTER" -> "gnomade_af_0.001;OTHER_FILTER"
      # "PASS;filter_vep_fail" -> "gnomade_af_0.001"
    echo "Step 3: Processing filter_vep_fail entries and adding custom filter: replace all filter_vep_fail label to gnomad_af label"
    awk -v filter_id="$FILTER_ID" '
    BEGIN {FS=OFS="\t"} 
    /^#/ {print $0; next} 
    !/^#/ {
        if ($7 ~ /filter_vep_fail/) {
            # Replace all occurrences of filter_vep_fail with the custom filter ID
            gsub(/filter_vep_fail/, filter_id, $7)
            # Remove PASS if it appears with the custom filter ID (e.g., "PASS;gnomade_af_0.001" -> "gnomade_af_0.001")
            gsub(/PASS;/, "", $7)
            gsub(/;PASS/, "", $7)
        }
        print $0
    }
    ' ~{step3_filter_vep_pass_processed_file} > ~{step4_filter_vep_fail_processed_file}

    # Step 4: Copy final result to output file
    echo "Step 4: Copying final result to output file..."
    cp ~{step4_filter_vep_fail_processed_file} ~{outfile}
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
