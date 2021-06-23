version 1.0

task pindelSomaticFilter {
  input {
    File reference
    File reference_fai
    File reference_dict
    File pindel_output_summary
  }

  runtime {
    memory: "16GB"
    docker: "mgibio/cle:v1.3.1"
  }

  String outdir = "/cromwell_root"
  command <<<
    /usr/bin/perl /usr/bin/write_pindel_filter_config.pl ~{pindel_output_summary} ~{reference} ~{outdir}
    /usr/bin/perl /usr/bin/somatic_indelfilter.pl filter.config
  >>>

  output {
    File vcf = "pindel.out.vcf"
  }
}
