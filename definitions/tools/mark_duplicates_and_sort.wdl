version 1.0

task markDuplicatesAndSort {
  input {
    File bam
    String input_sort_order = "queryname"
    String output_name = "MarkedSorted.bam"
  }
  String metrics_file_name = sub(output_name, "\.bam$", ".mark_dups_metrics.txt")
  Int space_needed_gb = 10 + round(5*size(bam, "GB"))
  #estimate 15M reads per Gb size of bam
  #markdup is listed as 2Gb per 100M reads
  Int mem_needed_gb = round(((size(bam, "GB")*15)/100)*2)+20
  runtime {
    docker: "mgibio/mark_duplicates-cwl:1.0.1"
    memory: "40GB"
    cpu: 8
    # add space to shift bam around via stdin/stdout and a bit more
    bootDiskSizeGb: space_needed_gb
    # add space for input bam, output bam, and a bit more
    disks: "local-disk ~{space_needed_gb} SSD"
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
  input {
    File bam
    String? input_sort_order
    String? output_name
  }
  call markDuplicatesAndSort {
    input:
    bam=bam,
    input_sort_order=input_sort_order,
    output_name=output_name
  }
}
