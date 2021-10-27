version 1.0

import "../tools/sequence_to_fastq_rna.wdl" as stfr
import "../tools/trim_fastq.wdl" as tf
import "../types.wdl"

workflow sequenceToTrimmedFastq {
  input {
    SequenceData unaligned
    File adapters
    String adapter_trim_end
    Int adapter_min_overlap
    Int max_uncalled
    Int min_readlength
    Boolean? unzip_fastqs
  }

  call stfr.sequenceToFastqRna as sequenceToFastq {
    input:
    bam=unaligned.sequence.bam,
    fastq1=unaligned.sequence.fastq1,
    fastq2=unaligned.sequence.fastq2,
    unzip_fastqs=unzip_fastqs
  }

  call tf.trimFastq {
    input:
    reads1=sequenceToFastq.fastq1,
    reads2=sequenceToFastq.fastq2,
    adapters=adapters,
    adapter_trim_end=adapter_trim_end,
    adapter_min_overlap=adapter_min_overlap,
    max_uncalled=max_uncalled,
    min_readlength=min_readlength,
  }

  output {
    Array[File] fastqs = trimFastq.fastqs
    File fastq1 = trimFastq.fastq1
    File fastq2 = trimFastq.fastq2
  }
}
