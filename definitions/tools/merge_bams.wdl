version 1.0

task mergeBams {
  input {
    Array[File] bams
    Boolean sorted = false
    String name = "merged"
  }

  Int cores = 4
  Int space_needed_gb = 10 + round(4*size(bams, "GB"))
  runtime {
    docker: "mgibio/bam-merge:0.1"
    memory: "8GB"
    cpu: cores
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  String outname = name + ".bam"
  command <<< #!/bin/bash
    set -o pipefail
    set -o errexit
    set -o nounset

    NUM_BAMS=~{length(bams)}
    #if there is only one bam, just copy it and index it
    if [[ $NUM_BAMS -eq 1 ]]; then
        cp "~{bams[0]}" "~{outname}";
    else
        if [[ "~{sorted}" == "true" ]];then
            /usr/bin/sambamba merge -t "~{cores}" "~{outname}" "~{S}{BAMS[@]}"
        else #unsorted bams, use picard
            java -jar -Xmx6g /opt/picard/picard.jar MergeSamFiles \
                OUTPUT="~{outname}" ASSUME_SORTED=true USE_THREADING=true \
                SORT_ORDER=unsorted VALIDATION_STRINGENCY=LENIENT \
                ~{sep=" " prefix("INPUT=", bams)}
        fi
    fi
    if [[ ~{sorted} == true ]];then
        /usr/bin/sambamba index "~{outname}"
    fi
  >>>

  output {
    File merged_bam = outname
  }
}

workflow wf {
  input { Array[File] bams }
  call mergeBams { input: bams=bams }
}
