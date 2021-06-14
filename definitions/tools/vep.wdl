version 1.0

import "types.wdl"

task vep {
  input {
    File vcf
    File cache_dir_zip
    File? synonyms_file
    Boolean coding_only = false
    # one of [pick, flag_pick, pick-allele, per_gene, pick_allele_gene, flag_pick_allele, flag_pick_allele_gene]
    String pick = "flag_pick"
    # TODO custom_annotations
    File reference
    File reference_fai
    File reference_dict
    Array[String] plugins
    Boolean everything = true
    String ensembl_assembly
    String ensembl_version
    String ensembl_species
  }
  runtime {
    memory: "64GB"
    bootDiskSizeGb: 30
    cpu: 4
    docker: "mgibio/vep_helper-cwl:vep_101.0_v2"
  }
  String annotated_vcf = basename(vcf, ".vcf") + "_annotated.vcf"
  command <<<
    mkdir VEP_cache && unzip -qq ~{cache_dir_zip} -d VEP_cache

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
    # TODO: custom annotations
    --fasta ~{reference}
    ~{sep(" ", prefix("--plugin", plugins))}  \ # TODO: prefix on all or just collection?
    ~{if everything then "--everything" else ""} \
    --assembly ~{ensembl_assembly} \
    --cache_version ~{ensembl_version} \
    --species ~{ensemble_species}
  >>>
  output {
    File annotated_vcf = annotated_vcf
    File vep_summary = basename(vcf, ".vcf") + "_annotated.vcf_summary.html"
  }
}
