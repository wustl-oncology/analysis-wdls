version 1.0

task arriba {
  input {
    File aligned_bam
    File reference_annotation
    File reference
    File reference_fai
    File reference_dict
  }

  Float bam_size_gb = size([aligned_bam, reference_annotation], "GB")
  Float reference_size_gb = size([reference, reference_fai, reference_dict], "GB")
  Int space_needed_gb = 10 + round(3*bam_size_gb + reference_size_gb)

  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "64GB"
    docker: "uhrigs/arriba:2.4.0"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  # explicit typing required, don't inline
  command <<<
    /arriba_v2.4.0/arriba \
    -b /arriba_v2.4.0/database/blacklist_hg38_GRCh38_v2.4.0.tsv.gz \
    -k /arriba_v2.4.0/database/known_fusions_hg38_GRCh38_v2.4.0.tsv.gz \
    -t /arriba_v2.4.0/database/known_fusions_hg38_GRCh38_v2.4.0.tsv.gz \
    -p /arriba_v2.4.0/database/protein_domains_hg38_GRCh38_v2.4.0.gff3 \
    -o arriba_fusions.tsv \
    -O arriba_fusions.discarded.tsv \
    -x ~{aligned_bam} \
    -g ~{reference_annotation} \
    -a ~{reference}
  >>>

  output {
    File fusion_predictions = "arriba_fusions.tsv"
    File discarded_fusion_predictions = "arriba_fusions.discarded.tsv"
  }
}

workflow wf { call arriba { input: } }
