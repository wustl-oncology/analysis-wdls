version 1.0

import "types.wdl"  # !UnusedImport

import "tools/pbmm2.wdl" as pb
import "tools/mosdepth.wdl" as md
import "tools/split_contigs.wdl" as sc
import "tools/deepsomatic.wdl" as ds
import "tools/clair3.wdl" as c3
import "tools/hiphase.wdl" as hp
import "tools/rename_vcf_annotations.wdl" as rn
import "tools/bcftools_norm.wdl" as bn
import "tools/combine_phased_proximal_vcf.wdl" as cp
import "tools/vep.wdl" as v
import "tools/bgzip.wdl" as bg
import "tools/index_vcf.wdl" as iv
import "tools/ensure_vcf_sample.wdl" as evs
import "subworkflows/vcf_readcount.wdl" as vr

# PacBio HiFi WGS somatic arm: tumor/normal bams -> VEP-annotated somatic vcf,
# annotated the way pVACseq needs it.
#
# This is the long-read counterpart of somatic_exome.wdl + detect_variants.wdl.
# Structural differences worth knowing before you read the calls:
#
#   * No BQSR. HiFi base qualities are already well calibrated and GATK's
#     recalibration model assumes short-read error profiles.
#   * No duplicate marking. HiFi library prep does not produce optical/PCR
#     duplicates in a way markdup can identify.
#   * One caller, not four. DeepSomatic's PACBIO model replaces the
#     Mutect2/Strelka2/VarScan2/DoCM ensemble; there is no HiFi consensus caller
#     to build an equivalent of CombineVariants around, and fpfilter's
#     assumptions (mismatch qualsum, strand bias) do not transfer.
#   * Phasing happens BEFORE annotation, not after. HiPhase phases the bam and
#     both vcfs in one pass, which is why phased_proximal_variants comes out of
#     the same step rather than needing a separate ReadBackedPhasing subworkflow.
#
# The VEP settings below are deliberately NOT the same as hifi-somatic-wdl's
# annotation.vep_annotate. See the comments on the vep call.
workflow somaticPacbioHifi {
  input {
    # ---------- Sequence inputs ----------
    # Unaligned (or already aligned, with skip_align=true) HiFi bams.
    # Supply exactly ONE of bams / fastqs per sample.
    # FASTQ mode exists to save disk -- a HiFi uBAM with kinetics is several
    # times the size of the equivalent fastq.gz. Nothing in this workflow reads
    # the PacBio tags that FASTQ discards (kinetics, MM/ML base modifications),
    # so the only real cost is that a methylation arm could never be added on
    # top of a FASTQ-derived run.
    Array[File] tumor_bams = []
    Array[File] tumor_fastqs = []
    Array[File] normal_bams = []
    Array[File] normal_fastqs = []
    String tumor_sample_name
    String normal_sample_name
    # Optional per-file @RG IDs, positionally matched to the fastq/bam arrays.
    # Keeps replicates distinguishable inside the merged bam (SM stays the
    # sample name; ID records which sequencing run each read came from).
    Array[String] tumor_read_group_ids = []
    Array[String] normal_read_group_ids = []
    # BAM mode only: treat inputs as already aligned and just merge/index.
    Boolean skip_align = false
    Boolean strip_kinetics = false

    # ---------- Reference ----------
    File reference
    File reference_fai
    File reference_dict
    # Optional pre-built pbmm2 index. Used for alignment only; every other step
    # takes the plain fasta.
    File? reference_mmi

    # ---------- Resources ----------
    Int pbmm2_cores = 24
    Int samtools_cores = 8
    Int deepsomatic_cores = 16
    Int clair3_cores = 16
    Int hiphase_cores = 24
    Int chunk_size = 75000000
    String clair3_model = "hifi_revio"
    Boolean run_mosdepth = true

    # ---------- VEP ----------
    File vep_cache_dir_zip
    String vep_ensembl_assembly = "GRCh38"
    String vep_ensembl_version
    String vep_ensembl_species = "homo_sapiens"
    File? vep_synonyms_file
    Array[VepCustomAnnotation] vep_custom_annotations = []  # !UnverifiedStruct
    Boolean vep_annotate_coding_only = false
    Boolean vep_everything = true
    # Frameshift and Wildtype are REQUIRED by pVACseq -- they write the
    # downstream frameshift peptide and the full wild-type protein sequence into
    # the CSQ. pvacseq run aborts immediately without them.
    Array[String] vep_plugins = ["Frameshift", "Wildtype"]
    # flag_pick, NOT pick. pVACseq applies its own transcript prioritization
    # (--transcript-prioritization-strategy, --maximum-transcript-support-level)
    # and needs to see every transcript to do it. --pick would collapse each
    # variant to a single consequence before pVACseq ever sees it.
    String vep_pick = "flag_pick"

    # ---------- Post-calling fixups ----------
    # DeepSomatic writes FORMAT/VAF; pVACtools reads FORMAT/AF for its
    # --tdna-vaf / --normal-vaf filters. Only affects the tumor column.
    Boolean rename_vaf_to_af = true

    # Generate DNA AD/AF/DP for BOTH tumor and normal with bam-readcount, the
    # same way detect_variants.wdl does for the short-read arm.
    #
    # This is on by default because DeepSomatic does not reliably emit normal
    # genotype fields, which leaves pVACseq's --normal-cov / --normal-vaf with
    # nothing to filter on. See the ensureNormalSample / dnaReadcounts calls
    # below for the mechanics.
    Boolean add_dna_readcounts = true
    Int? readcount_minimum_base_quality
    Int? readcount_minimum_mapping_quality
  }

  File align_reference = select_first([reference_mmi, reference])

  # ---------- 1. Alignment ----------

  call pb.alignHifiBams as alignTumor {
    input:
    bams=tumor_bams,
    fastqs=tumor_fastqs,
    read_group_ids=tumor_read_group_ids,
    reference=align_reference,
    reference_fai=reference_fai,
    sample_name=tumor_sample_name,
    skip_align=skip_align,
    strip_kinetics=strip_kinetics,
    pbmm2_cores=pbmm2_cores,
    samtools_cores=samtools_cores
  }

  call pb.alignHifiBams as alignNormal {
    input:
    bams=normal_bams,
    fastqs=normal_fastqs,
    read_group_ids=normal_read_group_ids,
    reference=align_reference,
    reference_fai=reference_fai,
    sample_name=normal_sample_name,
    skip_align=skip_align,
    strip_kinetics=strip_kinetics,
    pbmm2_cores=pbmm2_cores,
    samtools_cores=samtools_cores
  }

  # ---------- 2. Coverage QC ----------
  # Read these before choosing pVACseq's --tdna-cov / --normal-cov. HiFi WGS
  # typically runs 20-30x, where the exome-tuned short-read defaults over-filter.

  if (run_mosdepth) {
    call md.mosdepth as mosdepthTumor {
      input:
      bam=alignTumor.final_bam,
      bam_bai=alignTumor.final_bam_bai,
      output_basename=tumor_sample_name
    }
    call md.mosdepth as mosdepthNormal {
      input:
      bam=alignNormal.final_bam,
      bam_bai=alignNormal.final_bam_bai,
      output_basename=normal_sample_name
    }
  }

  # ---------- 3. Somatic small variants ----------

  call sc.splitContigs {
    input:
    reference_fai=reference_fai,
    chunk_size=chunk_size
  }

  scatter (ctg in splitContigs.contigs) {
    call ds.deepSomatic {
      input:
      tumor_bam=alignTumor.final_bam,
      tumor_bam_bai=alignTumor.final_bam_bai,
      normal_bam=alignNormal.final_bam,
      normal_bam_bai=alignNormal.final_bam_bai,
      reference=reference,
      reference_fai=reference_fai,
      regions_bed=ctg,
      tumor_sample_name=tumor_sample_name,
      normal_sample_name=normal_sample_name,
      output_basename=tumor_sample_name + "." + basename(ctg, ".bed"),
      cores=deepsomatic_cores
    }
  }

  call ds.gatherDeepSomatic {
    input:
    vcfs=deepSomatic.vcf,
    output_basename=tumor_sample_name + ".deepsomatic"
  }

  # ---------- 4. Germline small variants (phasing backbone + proximal variants) ----------

  call c3.clair3 as clair3Tumor {
    input:
    bam=alignTumor.final_bam,
    bam_bai=alignTumor.final_bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    sample_name=tumor_sample_name,
    clair_model=clair3_model,
    cores=clair3_cores
  }

  call c3.clair3 as clair3Normal {
    input:
    bam=alignNormal.final_bam,
    bam_bai=alignNormal.final_bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    sample_name=normal_sample_name,
    clair_model=clair3_model,
    cores=clair3_cores
  }

  # ---------- 5. Phasing ----------
  # Both vcfs are phased against the tumor bam in one invocation, so they come
  # out sharing PS ids -- that is what makes step 7 a simple concat.

  call hp.hiphaseSomatic as phaseTumor {
    input:
    bam=alignTumor.final_bam,
    bam_bai=alignTumor.final_bam_bai,
    germline_vcf=clair3Tumor.vcf,
    germline_vcf_tbi=clair3Tumor.vcf_tbi,
    somatic_vcf=gatherDeepSomatic.vcf,
    somatic_vcf_tbi=gatherDeepSomatic.vcf_tbi,
    reference=reference,
    reference_fai=reference_fai,
    output_basename=tumor_sample_name,
    cores=hiphase_cores
  }

  call hp.hiphaseGermline as phaseNormal {
    input:
    bam=alignNormal.final_bam,
    bam_bai=alignNormal.final_bam_bai,
    vcf=clair3Normal.vcf,
    vcf_tbi=clair3Normal.vcf_tbi,
    reference=reference,
    reference_fai=reference_fai,
    output_basename=normal_sample_name,
    cores=hiphase_cores
  }

  # ---------- 6. Normalize + annotate ----------

  if (rename_vaf_to_af) {
    call rn.renameVcfAnnotations as fixAf {
      input:
      input_vcf=phaseTumor.phased_somatic_vcf,
      renames=["FORMAT/VAF\tAF"]
    }
  }

  File somatic_for_norm = select_first([fixAf.renamed_vcf, phaseTumor.phased_somatic_vcf])

  call bn.bcftoolsNorm as normSomatic {
    input:
    input_vcf=somatic_for_norm,
    reference=reference,
    reference_fai=reference_fai
  }

  # Reuses tools/vep.wdl unmodified. The important part is the inputs:
  # Frameshift + Wildtype plugins, flag_pick (not pick), and an ENSEMBL cache --
  # do not point vep_cache_dir_zip at a --refseq cache, pVACtools expects
  # ENST/ENSG ids.
  call v.vep as annotateSomatic {
    input:
    vcf=normSomatic.normalized_vcf,
    cache_dir_zip=vep_cache_dir_zip,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    plugins=vep_plugins,
    ensembl_assembly=vep_ensembl_assembly,
    ensembl_version=vep_ensembl_version,
    ensembl_species=vep_ensembl_species,
    synonyms_file=vep_synonyms_file,
    custom_annotations=vep_custom_annotations,
    coding_only=vep_annotate_coding_only,
    everything=vep_everything,
    pick=vep_pick
  }

  # ---------- 6b. DNA readcounts for tumor AND normal ----------
  #
  # bam-readcount and vcf-readcount-annotator both look their sample up by name
  # and raise on a missing one, so the normal column has to exist before we can
  # write into it. ensureVcfSample creates an all-missing column when
  # DeepSomatic did not emit one, and passes through untouched when it did.
  #
  # vcfReadcount is subworkflows/vcf_readcount.wdl, reused verbatim from the
  # short-read arm -- it is technology agnostic, and it returns a bgzipped,
  # indexed vcf so no separate compress/index step is needed on this branch.

  if (add_dna_readcounts) {
    call evs.ensureVcfSample as ensureNormalSample {
      input:
      vcf=annotateSomatic.annotated_vcf,
      sample_name=normal_sample_name,
      output_basename=tumor_sample_name + ".annotated.with_normal"
    }

    call vr.vcfReadcount as dnaReadcounts {
      input:
      vcf=ensureNormalSample.vcf_with_sample,
      reference=reference,
      reference_fai=reference_fai,
      reference_dict=reference_dict,
      tumor_sample_name=tumor_sample_name,
      tumor_bam=phaseTumor.phased_bam,
      tumor_bam_bai=phaseTumor.phased_bam_bai,
      normal_sample_name=normal_sample_name,
      normal_bam=phaseNormal.phased_bam,
      normal_bam_bai=phaseNormal.phased_bam_bai,
      minimum_base_quality=readcount_minimum_base_quality,
      minimum_mapping_quality=readcount_minimum_mapping_quality
    }
  }

  # Fallback when add_dna_readcounts is false: tools/vep.wdl emits a plain
  # uncompressed .vcf with no index, and every downstream consumer wants
  # bgzip + tabix.
  if (!add_dna_readcounts) {
    call bg.bgzip as bgzipAnnotated {
      input: file=annotateSomatic.annotated_vcf
    }

    call iv.indexVcf as indexAnnotated {
      input: vcf=bgzipAnnotated.bgzipped_file
    }
  }

  # ---------- 7. Proximal variants vcf for pVACseq -p ----------

  call cp.combinePhasedProximalVcf as proximalVariants {
    input:
    germline_vcf=phaseTumor.phased_germline_vcf,
    germline_vcf_tbi=phaseTumor.phased_germline_vcf_tbi,
    somatic_vcf=phaseTumor.phased_somatic_vcf,
    somatic_vcf_tbi=phaseTumor.phased_somatic_vcf_tbi,
    tumor_sample_name=tumor_sample_name,
    output_basename=tumor_sample_name + ".phased_proximal_variants"
  }

  output {
    # ---- Alignments ----
    String tumor_input_mode = alignTumor.input_mode
    String normal_input_mode = alignNormal.input_mode
    File tumor_bam = alignTumor.final_bam
    File tumor_bam_bai = alignTumor.final_bam_bai
    File normal_bam = alignNormal.final_bam
    File normal_bam_bai = alignNormal.final_bam_bai
    File tumor_bam_phased = phaseTumor.phased_bam
    File tumor_bam_phased_bai = phaseTumor.phased_bam_bai
    File normal_bam_phased = phaseNormal.phased_bam
    File normal_bam_phased_bai = phaseNormal.phased_bam_bai

    # ---- Small variants ----
    File somatic_vcf = gatherDeepSomatic.vcf
    File somatic_vcf_tbi = gatherDeepSomatic.vcf_tbi
    File somatic_vcf_phased = phaseTumor.phased_somatic_vcf
    File somatic_vcf_phased_tbi = phaseTumor.phased_somatic_vcf_tbi
    File tumor_germline_vcf = clair3Tumor.vcf
    File normal_germline_vcf = clair3Normal.vcf
    File normal_germline_vcf_phased = phaseNormal.phased_vcf

    # ---- The files immuno_longread.wdl actually consumes ----
    File final_annotated_vcf = select_first([dnaReadcounts.readcounted_vcf, indexAnnotated.indexed_vcf])
    File final_annotated_vcf_tbi = select_first([dnaReadcounts.readcounted_vcf_tbi, indexAnnotated.indexed_vcf_tbi])
    File vep_summary = annotateSomatic.vep_summary

    # ---- DNA readcounts ----
    # "true" here means DeepSomatic did NOT emit a normal column and one was
    # created. Worth reading on run 1.
    String? normal_sample_column_was_added = ensureNormalSample.sample_was_added
    File? tumor_snv_bam_readcount_tsv = dnaReadcounts.tumor_snv_bam_readcount_tsv
    File? tumor_indel_bam_readcount_tsv = dnaReadcounts.tumor_indel_bam_readcount_tsv
    File? normal_snv_bam_readcount_tsv = dnaReadcounts.normal_snv_bam_readcount_tsv
    File? normal_indel_bam_readcount_tsv = dnaReadcounts.normal_indel_bam_readcount_tsv
    File phased_proximal_variants_vcf = proximalVariants.phased_vcf
    File phased_proximal_variants_vcf_tbi = proximalVariants.phased_vcf_tbi

    # ---- QC ----
    File? mosdepth_tumor_summary = mosdepthTumor.summary
    File? mosdepth_normal_summary = mosdepthNormal.summary
    File hiphase_tumor_stats = phaseTumor.stats
    File hiphase_tumor_summary = phaseTumor.summary
    File hiphase_tumor_blocks = phaseTumor.blocks
    File hiphase_normal_stats = phaseNormal.stats
  }
}
