version 1.0

import "types.wdl"

task vep {
  input {
    File vcf
    File cache_dir_zip
    File reference
    File reference_fai
    File reference_dict
    String ensembl_assembly
    String ensembl_version
    String ensembl_species
    Array[String] plugins
    Boolean coding_only = false
    Array[VepCustomAnnotation] custom_annotations = []
    Boolean everything = true
    # one of [pick, flag_pick, pick-allele, per_gene, pick_allele_gene, flag_pick_allele, flag_pick_allele_gene]
    String pick = "flag_pick"
    File? synonyms_file
  }

  Int space_needed_gb = 10 + round(size([vcf, reference, reference_fai, reference_dict, synonyms_file], "GB") + 2*size(cache_dir_zip, "GB"))
  runtime {
    memory: "64GB"
    bootDiskSizeGb: 30
    cpu: 4
    docker: "mgibio/vep_helper-cwl:vep_101.0_v2"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String annotated_vcf = basename(vcf, ".vcf") + "_annotated.vcf"
  String runtime_dir = "/cromwell_root"
  String cache_dir = runtime_dir + "/" + basename(cache_dir_zip, ".zip")
  # TODO: custom annotations
  command <<<
    mkdir ~{cache_dir} && unzip -qq ~{cache_dir_zip} -d ~{cache_dir}

    /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /usr/bin/variant_effect_predictor.pl \
    --format vcf \
    --vcf \
    --fork 4 \
    --term SO \
    --transcript_version \
    --offline \
    --cache \
    --symbol \
    -o ~{annotated_vcf} \
    -i ~{vcf} \
    ~{if defined(synonyms_file) then "--synonyms ~{synonyms_file}" else ""} \
    --coding_only ~{coding_only} \
    --~{pick} \
    --dir ~{cache_dir} \
    --fasta ~{reference} \
    ~{sep=" " prefix("--plugin ", plugins)}  \
    ~{if everything then "--everything" else ""} \
    --assembly ~{ensembl_assembly} \
    --cache_version ~{ensembl_version} \
    --species ~{ensembl_species}
  >>>

  output {
    File annotated_vcf = annotated_vcf
    File vep_summary = basename(vcf, ".vcf") + "_annotated.vcf_summary.html"
  }
}

workflow wf {
  input {
    File vcf
    File cache_dir_zip
    File reference
    File reference_fai
    File reference_dict
    Array[String] plugins
    String ensembl_assembly
    String ensembl_version
    String ensembl_species
    File? synonyms_file
    # TODO custom_annotations
    Array[VepCustomAnnotation] custom_annotations = []
    Boolean coding_only = false
    Boolean everything = true
    # one of [pick, flag_pick, pick-allele, per_gene, pick_allele_gene, flag_pick_allele, flag_pick_allele_gene]
    String pick = "flag_pick"
  }

  call vep {
    input:
    vcf=vcf,
    cache_dir_zip=cache_dir_zip,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    plugins=plugins,
    ensembl_assembly=ensembl_assembly,
    ensembl_version=ensembl_version,
    ensembl_species=ensembl_species,
    synonyms_file=synonyms_file,
    custom_annotations=custom_annotations,
    coding_only=coding_only,
    everything=everything,
    pick=pick,
  }
}
