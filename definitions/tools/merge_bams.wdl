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
  String S = "$"  # https://github.com/broadinstitute/cromwell/issues/1819
  command <<< #!/bin/bash
    set -o pipefail
    set -o errexit
    set -o nounset

    BAMS=(~{sep=" " bams})
    NUM_BAMS=~{length(bams)}
    #if there is only one bam, just copy it and index it
    if [[ $NUM_BAMS -eq 1 ]]; then
        cp "$BAMS" "~{outname}";
    else
        if [[ "~{sorted}" == "true" ]];then
            /usr/bin/sambamba merge -t "~{cores}" "~{outname}" "~{S}{BAMS[@]}"
        else #unsorted bams, use picard
            args=(OUTPUT="~{outname}" ASSUME_SORTED=true USE_THREADING=true SORT_ORDER=unsorted VALIDATION_STRINGENCY=LENIENT)
            for i in "~{S}{BAMS[@]}";do
                args+=("INPUT=$i")
            done
            java -jar -Xmx6g /opt/picard/picard.jar MergeSamFiles "~{S}{args[@]}"
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
