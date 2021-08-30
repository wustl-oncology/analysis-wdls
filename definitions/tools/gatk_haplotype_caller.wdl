version 1.0

task gatkHaplotypeCaller {
  input {
    File reference
    File reference_fai
    File reference_dict
    File bam
    File bai
    String emit_reference_confidence  # enum [NONE, BP_RESOLUTION, GVCF]
    Array[String] gvcf_gq_bands
    Array[String] intervals
    File? dbsnp_vcf
    File? dbsnp_vcf_tbi
    File verify_bam_id_metrics
    Int? max_alternate_alleles
    Int? ploidy
    String? read_filter
    String output_file_name
  }

  Float reference_size = size([reference, reference_fai, reference_dict], "GB")
  Float bam_size = size([bam, bai], "GB")
  Float vcf_size = size([dbsnp_vcf, dbsnp_vcf_tbi], "GB")
  Int space_needed_gb = 10 + round(reference_size + 2*bam_size + vcf_size)
  runtime {
    memory: "18GB"
    docker: "broadinstitute/gatk:4.1.8.1"
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  # TODO: also check alphanumeric, /^0-9A-Za-z]+$/
  command <<<
    # requires .bai not .bam.bai
    mv ~{bam} ~{basename(bam)}
    mv ~{bai} ~{basename(basename(bai, ".bai"), ".bam") + ".bai"}

    # get FREEMIX value for contamination fraction
    python <<CODE
    with open("~{verify_bam_id_metrics}", "r") as f:
        header = f.readline().split("\t")
        if len(header) >= 7 and header[6] == "FREEMIX":
            with open("freemix.txt", "w") as w:
                w.write(f.readline().split("\t")[6])
    CODE
    CONTAMINATION_FRACTION=`cat freemix.txt 2> /dev/null`
    if [ -z CONTAMINATION_FRACTION ]; then
        CONTAMINATION_ARG="--contamination $CONTAMINATION_FRACTION"
    fi

    # do the task itself
    /gatk/gatk --java-options -Xmx16g HaplotypeCaller \
    -R ~{reference} \
    -I ~{basename(bam)} \
    -ERC ~{emit_reference_confidence} \
    ~{sep=" " prefix("-GQB", gvcf_gq_bands)} \
    -L ~{sep="," intervals} \
    ~{if(defined(dbsnp_vcf)) then "--dbsnp " + dbsnp_vcf else ""} \
    $CONTAMINATION_ARG \
    ~{if(defined(max_alternate_alleles)) then "--max_alternate_alleles " + max_alternate_alleles else ""} \
    ~{if(defined(ploidy)) then "-ploidy " + ploidy else ""} \
    ~{if(defined(read_filter)) then "--read_filter " + read_filter else ""} \
    -O ~{output_file_name}
  >>>

  output {
    File gvcf = output_file_name
    File gvcf_tbi = output_file_name + ".tbi"
  }
}

workflow wf { call gatkHaplotypeCaller { input: } }
