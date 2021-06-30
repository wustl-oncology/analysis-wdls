version 1.0

task filterVcfDepth {
  input {
    File vcf
    Int minimum_depth
    Array[String] sample_names
  }

  runtime {
    docker: "mgibio/depth-filter:0.1.2"
    memory: "4GB"
  }

  String outfile = "depth_filtered.vcf"
  command <<<
    /opt/conda/bin/python3 /usr/bin/depth_filter.py \
    --minimum_depth ~{minimum_depth} \
    ~{vcf} ~{sep="," sample_names} \
    ~{outfile}
  >>>

  output {
    File depth_filtered_vcf = outfile
  }
}
