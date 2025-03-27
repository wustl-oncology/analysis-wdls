version 1.0

task addVepFieldsToTable {
  input {
    File vcf
    Array[String] vep_fields = ["Consequence", "SYMBOL", "Feature", "HGVSc", "HGVSp"]
    File? tsv
    String prefix = "variants"
    File? preferred_transcripts_tsv
  }

  Int space_needed_gb = 10 + round(size([vcf, tsv], "GB")*2)
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "4GB"
    docker: "griffithlab/vatools:5.2.0"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    vep-annotation-reporter -o ~{prefix}.annotated.tsv \
    ~{vcf} ~{sep=" " vep_fields} \
    ~{if defined(tsv) then "-t ~{tsv}" else ""} \
    ~{if defined(preferred_transcripts_tsv) then "-p ~{preferred_transcripts_tsv}" else ""}
  >>>

  output {
    File annotated_variants_tsv = "~{prefix}.annotated.tsv"
  }
}

workflow wf {
  input {
    File vcf
    Array[String]? vep_fields
    File? tsv
    String? prefix
    File? preferred_transcripts_tsv
  }

  call addVepFieldsToTable {
    input:
    vcf=vcf,
    vep_fields=vep_fields,
    tsv=tsv,
    prefix=prefix,
    preferred_transcripts_tsv=preferred_transcripts_tsv
  }
}
