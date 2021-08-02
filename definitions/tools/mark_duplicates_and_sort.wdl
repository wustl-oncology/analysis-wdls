version 1.0

task markDuplicatesAndSort {
  input {
    File bam
    String input_sort_order = "queryname"
    String output_name = "MarkedSorted.bam"
  }
  String metrics_file_name = sub(output_name, "\.bam$", ".mark_dups_metrics.txt")

  Float bam_size_gb = size(bam, "GB")
  runtime {
    docker: "mgibio/mark_duplicates-cwl:1.0.1"
    memory: "40GB"
    cpu: 8
    # add space to shift bam around via stdin/stdout and a bit more
    bootDiskSizeGb: 10 + round(bam_size_gb * 3)
    # add space for input bam, output bam, and a bit more
    disks: "local-disk ~{10 + round(bam_size_gb * 3)} SSD"
  }

  command <<<
    set -o pipefail
    set -o errexit
    /usr/bin/java -Xmx16g -jar /opt/picard/picard.jar MarkDuplicates I=~{bam} O=/dev/stdout ASSUME_SORT_ORDER=~{input_sort_order} METRICS_FILE=~{metrics_file_name} QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT \
        | /usr/bin/sambamba sort -t 8 -m 18G -o ~{output_name} /dev/stdin
  >>>

  output {
    File sorted_bam = output_name
    File sorted_bam_bai = output_name + ".bai"
    File metrics_file = metrics_file_name
  }
}

workflow wf {
  input { File bam }
  call markDuplicatesAndSort { input: bam=bam }
}
