version 1.0

task filterVcfCodingVariant {
  input {
    File vcf
    String? coding_filter
  }

  Int space_needed_gb = 10 + round(2*size(vcf, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "4GB"
    docker: "mgibio/vep_helper-cwl:vep_113.3_v1"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String outfile = "annotated.coding_variant_filtered"

  #Note on the '--ontology' option and how this tool filters variants using ensembl filter_vep
  #If the '--ontology' option is used with ensembl filter_vep, the tool will query a SequenceOntology web service to construct the --filter query
  #This service will be used to find all SO terms that are children of coding_sequence_variant
  #This approach assumes that the variants will be annotated by VEP using Sequence Ontology terms (e.g. --terms SO) which is default for VEP
  #If this tool is being run in a secure environment that prevents access to external services (e.g. due to firewall settings, private cloud network, etc.) it will fail.
  #To prevent use of the external web service, the user can define a specific query directly using the 'coding_filter' WDL input
  #If the WDL input 'coding_filter' is defined, it will be used verbatim as the --filter query for ensembl filter_vep
  #Refer to the ensembl filter_vep documentation for instructions on how to construct a valid --filter query
  #The following coding_filter string should be equivalent to use of the --ontology approach:

  #coding_filter: "Consequence in coding_sequence_variant,initiator_codon_variant,start_lost,start_retained_variant,protein_altering_variant,frameshift_variant,frame_restoring_variant,frameshift_elongation,frameshift_truncation,frameshift_variant_NMD_escaping,frameshift_variant_NMD_triggering,minus_1_frameshift_variant,minus_2_frameshift_variant,plus_1_frameshift_variant,plus_2_frameshift_variant,inframe_variant,incomplete_terminal_codon_variant,inframe_indel,inframe_deletion,conservative_inframe_deletion,disruptive_inframe_deletion,inframe_insertion,conservative_inframe_insertion,disruptive_inframe_insertion,nonsynonymous_variant,missense_variant,conservative_missense_variant,non_conservative_missense_variant,rare_amino_acid_variant,pyrrolysine_loss,selenocysteine_gain,selenocysteine_loss,redundant_inserted_stop_gained,start_lost,stop_gained,stop_gained_NMD_escaping,stop_gained_NMD_triggering,stop_lost,synonymous_variant,start_retained_variant,stop_retained_variant,terminator_codon_variant,incomplete_terminal_codon_variant,stop_lost,stop_retained_variant or Consequence match coding_sequence_variant" 

  # Prepare filter options in advance
  String filter_opts = if defined(coding_filter) then "--filter '~{coding_filter}'" else "--filter 'Consequence is coding_sequence_variant' --ontology"

  command <<<
    /usr/bin/perl /usr/bin/vcf_check.pl ~{vcf} ~{outfile} \
    /usr/bin/perl /opt/vep/src/ensembl-vep/filter_vep \
    --format vcf \
    -o ~{outfile} \
    ~{filter_opts} \
    -i ~{vcf}
  >>>

  output {
    File filtered_vcf = outfile
  }
}

workflow wf {
  input {
    File vcf
    String? coding_filter
  }

  call filterVcfCodingVariant {
    input:
    vcf=vcf,
    coding_filter=coding_filter
  }
}

