version 1.0

# TODO: add trimming inputs

task sequenceAlignAndTag {
  input {
    File? bam
    File? fastq1
    File? fastq2
    String readgroup
    File reference
    # secondary files. Must be separate in WDL
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_sa
  }
  Int cores = 8

  runtime {
    docker: "mgibio/alignment_helper-cwl:1.1.0"
    memory: "20GB"
    cpu: cores
    # 1 + just for a buffer
    # size(bam, "GB")*10 because bam uncompresses and streams to /dev/stdout and /dev/stdin, could have a couple flying at once
    bootDiskSizeGb: 1 + round(size(bam, "GB")*10 + size([reference, reference_amb, reference_ann, reference_bwt, reference_pac, reference_sa], "GB"))
  }

  command <<<
    set -o pipefail
    set -o errexit
    set -o nounset

    [[ ! -z "~{bam}" ]] && MODE="bam"
    [[ ! -z "~{fastq1}" ]] && [[ ! -z "~{fastq2}" ]] && MODE="fastq"
    RUN_TRIMMING="false"

    if [[ "$MODE" == "fastq" ]]; then
        if [[ "$RUN_TRIMMING" == 'false' ]]; then
            /usr/local/bin/bwa mem -K 100000000 -t ~{cores} -Y -R "~{readgroup}" "~{reference}" "~{fastq1}" "~{fastq2}" \
              | /usr/local/bin/samblaster -a --addMateTags | /opt/samtools/bin/samtools view -b -S /dev/stdin
        else
            /opt/flexbar/flexbar --adapters "$TRIMMING_ADAPTERS" --reads "~{fastq1}" --reads2 "~{fastq2}" --adapter-trim-end LTAIL --adapter-min-overlap "$TRIMMING_ADAPTER_MIN_OVERLAP" --adapter-error-rate 0.1 --max-uncalled 300 --stdout-reads \
              | /usr/local/bin/bwa mem -K 100000000 -t ~{cores} -Y -p -R "~{readgroup}" "~{reference}" /dev/stdin \
              | /usr/local/bin/samblaster -a --addMateTags \
              | /opt/samtools/bin/samtools view -b -S /dev/stdin
        fi
    fi
    if [[ "$MODE" == "bam" ]]; then
        if [[ $RUN_TRIMMING == "false" ]]; then
            /usr/bin/java -Xmx4g -jar /opt/picard/picard.jar SamToFastq I="~{bam}" INTERLEAVE=true INCLUDE_NON_PF_READS=true FASTQ=/dev/stdout \
              | /usr/local/bin/bwa mem -K 100000000 -t ~{cores} -Y -p -R "~{readgroup}" "~{reference}" /dev/stdin \
              | /usr/local/bin/samblaster -a --addMateTags \
              | /opt/samtools/bin/samtools view -b -S /dev/stdin
        else
            /usr/bin/java -Xmx4g -jar /opt/picard/picard.jar SamToFastq I="~{bam}" INTERLEAVE=true INCLUDE_NON_PF_READS=true FASTQ=/dev/stdout \
              | /opt/flexbar/flexbar --adapters "$TRIMMING_ADAPTERS" --reads - --interleaved --adapter-trim-end LTAIL --adapter-min-overlap "$TRIMMING_ADAPTER_MIN_OVERLAP" --adapter-error-rate 0.1 --max-uncalled 300 --stdout-reads \
              | /usr/local/bin/bwa mem -K 100000000 -t ~{cores} -Y -p -R "~{readgroup}" "~{reference}" /dev/stdin \
              | /usr/local/bin/samblaster -a --addMateTags \
              | /opt/samtools/bin/samtools view -b -S /dev/stdin
        fi
    fi
  >>>

  output {
    File aligned_bam = stdout()
  }
}


workflow wf {
  input {
    File? bam
    File? fastq1
    File? fastq2
    String readgroup
    File reference
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_sa
  }
  call doTask {
    input:
    bam=bam,
    readgroup=readgroup,
    reference=reference,
    reference_amb=reference_amb,
    reference_ann=reference_ann,
    reference_bwt=reference_bwt,
    reference_pac=reference_pac,
    reference_sa=reference_sa
  }
}
