version 1.0

# Rename FILTER / INFO / FORMAT keys in place.
#
# The reason this exists: DeepSomatic writes the variant allele fraction as
# FORMAT/VAF, but pVACtools reads FORMAT/AF for its --tdna-vaf and --normal-vaf
# coverage filters. Without this rename those filters silently see nothing and
# every candidate passes the VAF gate.
#
# Each entry in `renames` is a "TYPE/OLD<TAB>NEW" line as understood by
# `bcftools annotate --rename-annots`, e.g. "FORMAT/VAF\tAF". Pass them as plain
# strings with a literal tab, or use the default below.
task renameVcfAnnotations {
  input {
    File input_vcf
    Array[String] renames = ["FORMAT/VAF\tAF"]
    Int cores = 4
  }

  String outbase = basename(basename(input_vcf, ".gz"), ".vcf")
  Int space_needed_gb = 10 + round(4*size(input_vcf, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"
    memory: "8GB"
    cpu: cores
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    set -euxo pipefail
    bcftools --version

    printf '%b\n' "~{sep='" "' renames}" > rename_annots.txt
    cat rename_annots.txt

    bcftools annotate \
      --threads ~{cores - 1} \
      --rename-annots rename_annots.txt \
      -Oz -o ~{outbase}.renamed.vcf.gz \
      ~{input_vcf}

    bcftools index --threads ~{cores - 1} -t ~{outbase}.renamed.vcf.gz
  >>>

  output {
    File renamed_vcf = outbase + ".renamed.vcf.gz"
    File renamed_vcf_tbi = outbase + ".renamed.vcf.gz.tbi"
  }
}

workflow wf {
  input {
    File input_vcf
    Array[String]? renames
  }
  call renameVcfAnnotations {
    input:
    input_vcf=input_vcf,
    renames=select_first([renames, ["FORMAT/VAF\tAF"]])
  }
}
