version 1.0

import "tools/minimap2_rnaseq_pacbioHifi.wdl" as mm
import "tools/stringtie3.wdl" as st3
import "tools/samtools_flagstat.wdl" as sf

# PacBio Iso-Seq RNA arm for the long-read immuno pipeline.
#
# This is the long-read counterpart of rnaseq_star_fusion.wdl. Deliberately much
# smaller: no adapter trimming (FLNC is already primer/polyA trimmed by
# `isoseq refine`), no duplicate marking (PCR duplicates are not a meaningful
# concept for FLNC), no STAR two-pass.
#
#   flnc.bam --minimap2 splice:hq--> aligned bam --stringtie3 -L -e--> expression
#
# The two expression files are consumed by tools/vcf_expression_annotator.wdl
# with expression_tool="stringtie":
#   gene mode       <- gene_expression_tsv (writes FORMAT/GX)
#   transcript mode <- transcript_gtf      (writes FORMAT/TX)
#
# aligned_bam is also a valid --LR_bam input for pvacfuse_longread.wdl, so the
# fusion arm can reuse this alignment instead of doing its own.
workflow rnaseqPacbioHifi {
  input {
    # One or more unaligned FLNC bams from `isoseq refine`. Multiple entries are
    # replicates of ONE biological sample: each is aligned with its own @RG and
    # the shards are merged into a single sorted, indexed bam. That single bam
    # then feeds StringTie, regtools and ctat-LR-fusion alike.
    Array[File] flnc_bams
    # Optional per-replicate @RG IDs, positionally matched to flnc_bams.
    # e.g. ["NTR004_rna_rep1", "NTR004_rna_rep2", "NTR004_rna_rep3"]
    Array[String] flnc_read_group_ids = []
    String sample_name

    # Either the genome fasta or a prebuilt .mmi (faster: skips genome indexing
    # on every run). If you pass an .mmi, confirm it was built with the same
    # preset -- see the caveat in tools/minimap2_rnaseq_pacbioHifi.wdl.
    File reference
    File reference_fai
    # Gene model GTF. Must be the same Ensembl release as the VEP cache used on
    # the DNA arm, otherwise transcript ids will not join and TX stays empty.
    File reference_annotation

    String minimap2_preset = "splice:hq"
    # Extra minimap2 flags. Default empty -- the base command already matches a
    # configuration proven against ctat-LR-fusion.
    String minimap2_additional_args = ""
    Int minimap2_cores = 16
    Int samtools_cores = 8
    String? minimap2_docker_image

    # FLNC is strand resolved and minimap2 is run with -uf, so the alignments
    # are already oriented; "unstranded" is correct here.
    String stringtie_strand = "unstranded"
    # true => quantify reference transcripts only (no MSTRG.* novel ids).
    # Required for the ids to match VEP's CSQ. Set false only for discovery runs.
    Boolean stringtie_reference_only = true
    Int stringtie_cores = 12
    String? stringtie_docker_image
  }

  call mm.minimap2Align as align {
    input:
    flnc_bams=flnc_bams,
    read_group_ids=flnc_read_group_ids,
    sample_name=sample_name,
    reference=reference,
    reference_fai=reference_fai,
    output_basename=sample_name + ".flnc.aligned",
    preset=minimap2_preset,
    additional_args=minimap2_additional_args,
    minimap2_cores=minimap2_cores,
    samtools_cores=samtools_cores,
    minimap2_docker_image=minimap2_docker_image
  }

  call st3.stringtie3 as quantify {
    input:
    bam=align.aligned_bam,
    bam_bai=align.aligned_bam_bai,
    reference_annotation=reference_annotation,
    sample_name=sample_name,
    strand=stringtie_strand,
    reference_only=stringtie_reference_only,
    cores=stringtie_cores,
    docker_image=select_first([stringtie_docker_image, "quay.io/biocontainers/stringtie:3.0.3--h29c0135_0"])
  }

  call sf.samtoolsFlagstat as flagstat {
    input:
    bam=align.aligned_bam,
    bam_bai=align.aligned_bam_bai
  }

  output {
    File final_bam = align.aligned_bam
    File final_bam_bai = align.aligned_bam_bai

    # -> vcf-expression-annotator ... stringtie transcript  (FORMAT/TX)
    File stringtie_transcript_gtf = quantify.transcript_gtf
    # -> vcf-expression-annotator ... stringtie gene        (FORMAT/GX)
    File stringtie_gene_expression_tsv = quantify.gene_expression_tsv

    File flagstats = flagstat.flagstats
  }
}
