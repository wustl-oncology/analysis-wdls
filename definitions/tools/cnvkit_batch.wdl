version 1.0

task cnvkitBatch {
  input {
    File tumor_bam
    File? bait_intervals
    File? access
    File? normal_bam
    File reference
    String method = "hybrid"  # enum [hybrid, amplicon, wgs]
    Boolean diagram = false
    Boolean scatter_plot = false
    Boolean drop_low_coverage = false
    Boolean male_reference = false
    Int? target_average_size
  }


  Int size_needed_gb = 10 + round(size([tumor_bam, bait_intervals, access, normal_bam, reference], "GB") * 2)
  runtime {
    bootDiskSizeGb: 10
    memory: "4GB"
    cpu: 1
    docker: "etal/cnvkit:0.9.5"
    disks: "local-disk ~{size_needed_gb} SSD"
  }

  command <<<
    /usr/bin/python /usr/local/bin/cnvkit.py batch \
    ~{tumor_bam} \
    ~{if defined(normal_bam) then "--normal ~{normal_bam}" else ""} \
    --fasta ~{reference} \
    ~{if defined(bait_intervals) then "--targets ~{bait_intervals}" else ""} \
    ~{if defined(access) then "--access ~{access}" else ""} \
    --method ~{method} \
    ~{if diagram then "--diagram" else ""} \
    ~{if scatter_plot then "--scatter" else ""} \
    ~{if drop_low_coverage then "--drop-low-coverage" else ""} \
    ~{if male_reference then "--male-reference" else ""} \
    ~{if defined(target_average_size) then "--target-avg-size ~{target_average_size}" else ""}
  >>>

  String intervals_base = basename(if defined(bait_intervals) then "~{bait_intervals}" else "", ".interval_list")
  String normal_base = basename(if defined(normal_bam) then "~{normal_bam}" else "", ".bam")
  output {
    File? intervals_antitarget = intervals_base + ".antitarget.bed"
    File? intervals_target = intervals_base + ".target.bed"
    File? normal_antitarget_coverage = normal_base + ".antitarget.bed"
    File? normal_target_coverage = normal_base + ".target.bed"
    File? reference_coverage = "reference.cnn"
    File? cn_diagram = basename(tumor_bam, ".bam") + "-diagram.pdf"
    File? cn_scatter_plot = basename(tumor_bam, ".bam") + "-scatter.pdf"
    File tumor_antitarget_coverage = basename(tumor_bam, ".bam") + ".antitargetcoverage.cnn"
    File tumor_target_coverage = basename(tumor_bam, ".bam") + ".targetcoverage.cnn"
    File tumor_bin_level_ratios = basename(tumor_bam, ".bam") + ".cnr"
    File tumor_segmented_ratios = basename(tumor_bam, ".bam") + ".cns"
  }
}

workflow wf { call cnvkitBatch { input: } }
