version 1.0

import "../types.wdl"  # !UnusedImport
import "../tools/sequence_align_and_tag.wdl" as saat
import "../tools/merge_bams.wdl" as mb
import "../tools/name_sort.wdl" as ns
import "../tools/mark_duplicates_and_sort.wdl" as mdas
import "../tools/index_bam.wdl" as ib
import "../tools/bqsr.wdl" as b
import "../tools/apply_bqsr.wdl" as ab

workflow sequenceToBqsr {
  input {
    Array[SequenceData] unaligned

    Array[File] bqsr_known_sites
    Array[File] bqsr_known_sites_tbi
    Array[String]? bqsr_intervals

    TrimmingOptions? trimming

    File reference
    File reference_fai
    File reference_dict
    File reference_alt
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_sa

    String final_name = "final"
  }

  scatter(seq_data in unaligned) {
    call saat.sequenceAlignAndTag {
      input:
      unaligned=seq_data,
      trimming=trimming,
      reference=reference,
      reference_alt=reference_alt,
      reference_amb=reference_amb,
      reference_ann=reference_ann,
      reference_bwt=reference_bwt,
      reference_pac=reference_pac,
      reference_sa=reference_sa
    }
  }

  call mb.mergeBams {
    input:
    bams=sequenceAlignAndTag.aligned_bam,
    name=final_name
  }

  call ns.nameSort {
    input: bam=mergeBams.merged_bam
  }

  call mdas.markDuplicatesAndSort {
    input: bam=nameSort.name_sorted_bam
  }

  call b.bqsr {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,

    bam=markDuplicatesAndSort.sorted_bam,
    bam_bai=markDuplicatesAndSort.sorted_bam_bai,

    intervals=bqsr_intervals,

    known_sites=bqsr_known_sites,
    known_sites_tbi=bqsr_known_sites_tbi
  }

  call ab.applyBqsr {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,

    bam=markDuplicatesAndSort.sorted_bam,
    bam_bai=markDuplicatesAndSort.sorted_bam_bai,

    bqsr_table=bqsr.bqsr_table,
    output_name=final_name
  }

  call ib.indexBam {
    input: bam = applyBqsr.bqsr_bam
  }

  output {
    File final_bam = indexBam.indexed_bam
    File final_bam_bai = indexBam.indexed_bam_bai
    File final_bai = indexBam.indexed_bai
    File mark_duplicates_metrics_file = markDuplicatesAndSort.metrics_file
  }
}
