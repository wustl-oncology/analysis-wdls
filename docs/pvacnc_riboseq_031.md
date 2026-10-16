# pVACnc 0.3.1: noncanonical ORF coverage

Both ORFanage and TransDecoder use the same implementation. Preparation saves
each mapped ORF's CDS blocks after subtracting the union of reference CDS on
the same chromosome and strand. Coordinates are 1-based inclusive. ORFs are
keyed by caller, with-e/no-e origin and ORF ID. Missing/invalid models are
recorded in `model_status` without changing the generated VCF or mapping.

After the existing coordinate checker, the supplemental checker reports:

- Level 1: existing SNV ALT support.
- Level 2: existing variant-coordinate coverage when level 1 is absent.
- Level 3: any aligned-base coverage in the same ORF's noncanonical CDS blocks
  when levels 1 and 2 are absent.

No P-site, frame, read-strand or region-base-quality check is performed.
MAPQ and alignment-flag settings match the current coordinate checker.
CIGAR deletions/skips are not aligned bases. An alignment spanning several
blocks is counted once per ORF. Coverage is evidence at a region, not proof
that the predicted ORF or mutant peptide is translated. Levels are independent
of pVACseq filtering and pVACview tiers. Existing tables/counts are preserved.
All mapping candidates are assessed, retaining the existing filtered-membership
columns; the legacy filtered-report join is unchanged.

New outputs per caller:

- Preparation: `results/<sample>.noncanonical_orf_blocks.tsv`.
- Ribo: `*.riboseq_orf_support_by_bam.tsv`.
- Ribo: `*.riboseq_orf_candidate_summary.tsv`.

## Build on the cluster

```bash
cd /storage1/fs1/mgriffit/Active/griffithlab/adhoc/allen/MS_TA_work/Noncanonical_ORF_Variant_Projection_wdl/pvacnc_container/pvacnc
bsub -G compute-oncology -g /jinglun/wdl -q oncology-interactive -Is \
  -a 'docker_build(jinglunli/pvacnc)' -- --tag jinglunli/pvacnc:0.3.1 .
```

After the build succeeds:

```bash
bsub -G compute-oncology -g /jinglun/wdl -q general-interactive -Is \
  -n 1 -M 4GB -R 'rusage[mem=4GB]' \
  -a 'docker(jinglunli/pvacnc:0.3.1)' -- /bin/sh -c '
    PATH=/usr/local/bin:/usr/bin:/bin
    export PATH
    pvacnc-noncanonical-orf-blocks --help >/dev/null &&
    pvacnc-riboseq-orf-support --help >/dev/null &&
    python3 -B /opt/pvacnc/test_riboseq_031.py &&
    echo OK
  '
```

## Reuse existing results

Existing pVACseq results do not need recomputing. Inside the new image, use
`pvacnc-noncanonical-orf-blocks --caller orfanage` (or `transdecoder`) with
`--mapping`, `--reference-gtf`, `--with-e`, `--no-e`, and `--output`.
For ORFanage the two model inputs are GTFs; for TransDecoder they are genome
GFF3 files. Then use `pvacnc-riboseq-orf-support` with the existing
`--candidate-summary`, new `--noncanonical-blocks`, repeated
`--bam`/`--bam-index`/`--label` arguments, and a new `--output-prefix`.
Use the same BAMs and MAPQ/flag settings as the existing coordinate check.
This requires neither ORF calling nor pVACseq to run again.

For a new full workflow, the modified WDL supplies the block file automatically
when Ribo support is enabled. For a standalone Ribo WDL run, provide the new
required `noncanonical_orf_blocks` input. Mouse module runs continue to use
the standalone `definitions/pvacnc.wdl` entry point.

## Validation status

Six synthetic tests passed locally, covering subtraction, both strands,
both callers, skipped models, origin isolation, CIGAR overlap, filters and
level precedence. Container build, full WDL engine validation and real BAM
integration remain for the cluster. No existing results were overwritten.
Use a new output directory and tag 0.3.1; do not replace tag 0.3.0.
