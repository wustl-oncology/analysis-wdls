# pVACnc

## Inputs and selection

Run `definitions/pvacnc.wdl` standalone, or enable the human immuno branch
with `enable_pvacnc: true`.

Required upstream inputs are the tumor RNA BAM, StringTie `-e` GTF,
expression-annotated VEP VCF, MHC alleles and matched genome, annotation and
peptide references. A phased VCF is optional.

Both callers:

1. Retain `FILTER=PASS` ALT alleles from the source VCF.
2. Exclude alleles whose reference span overlaps any annotated reference CDS.
3. Project variants into predicted ORFs from both transcript branches.
4. Retain distinct no-`-e` protein-changing events; prefer successful `-e`
   events for identical WT/MT protein sequences of the same allele and strand.
5. Generate a derived VCF and run pVACseq and pVACview annotation.

Routing is independent between callers. An ORFanage hit does not suppress
a TransDecoder result.

Canonical pVACseq prediction reports are not required. The integrated human
workflow still obtains the expression-annotated VCF from the canonical
pVACseq subworkflow. Moving that annotation upstream is needed for independent
scheduling.

## Containers

| Tool | Default image |
|---|---|
| ORFanage | `quay.io/biocontainers/orfanage:1.2.0--heaafb18_2` |
| TransDecoder | `quay.io/biocontainers/transdecoder:6.0.0--pl5321hdfd78af_0` |
| Preparation and pVACview annotation | `jinglunli/pvacnc:0.3.5` |
| Ribo-seq support | `jinglunli/pvacnc:0.3.5` |

Images can be overridden through workflow inputs. Preparation and Ribo-seq
image versions are independent of the WDL interface version.

The preparation image is maintained outside this repository. It includes
Python, `bgzip`, `tabix` and three entry points:

- `pvacnc-prepare`
- `pvacnc-transdecoder-prepare`
- `pvacnc-annotate-pvacview`

The wrappers use image-local executables rather than a host Conda environment.

## ORF prediction and translation

ORFanage uses `BEST`. TransDecoder uses strand-specific prediction,
a configurable minimum ORF length (default 30 aa), partial and complete
ORFs, and multiple ORFs per transcript. No BLAST/Pfam retention evidence
is supplied.

Both callers translate frameshift alleles through the full spliced parent
transcript to the first new stop codon or transcript end. ORFanage checks
that reconstructed wild-type translation agrees with its CDS translation.

Candidate-specific model, projection and VCF-construction failures are
recorded in caller-specific skip tables. Conditional skipping of downstream
prediction for zero candidates remains unimplemented.

## Annotation and expression

With-`-e` ORFs matched to a source VEP transcript retain its reference
transcript and gene metadata, including MANE Select, Canonical and TSL.
They retain the source WDL/Kallisto `TX` and `GX` values.

No-`-e` rescue ORFs retain sample-specific identifiers and StringTie expression.
Unique coordinate-level VEP gene symbols are joined with `/`. Reference
transcript fields remain unset; full VEP contexts are retained in the mapping.

Synthetic translated ORFs use `BIOTYPE=protein_coding` for pVACseq
compatibility. The source reference biotype does not replace that field.
ORF provenance remains in the mapping; predicted-ORF consequences are not
ordinary reference-CDS VEP annotations.

In pVACview, no-`-e` MANE Select, Canonical and TSL fields remain empty.
`Transcript Pass` displays `de novo` and is exempt only from the
reference-transcript-quality tier check. Other tier criteria are unchanged.

## Outputs

Integrated outputs are separated into `pvacnc/orfanage/` and
`pvacnc/transdecoder/`.

Standalone model outputs use matching branch folders:

```text
with_e_orfanage_outputs/<sample>.with_e.orfanage.gtf
no_e_orfanage_outputs/<sample>.no_e.orfanage.gtf
with_e_transdecoder_outputs/<sample>.with_e.transdecoder.genome.gff3
no_e_transdecoder_outputs/<sample>.no_e.transdecoder.genome.gff3
```

TransDecoder also exports cDNA, CDS, peptide, transcript-GFF3 and summary
files. Preparation files use `<sample>.orfanage.*` or
`<sample>.transdecoder.*`, including VCF/index, mapping, routing, effects,
noncanonical blocks, skip reports and audits.

Prediction reports are in caller-specific MHC-I, MHC-II and combined
directories. Standard pVACseq/pVACview filenames are retained within those
directories. Prediction logs use `orfanage_prediction_logs/` and
`transdecoder_prediction_logs/`.

## Optional Ribo-seq support

Enable `enable_pvacnc_riboseq_support` in human immuno or
`enable_riboseq_support` in the standalone workflow. Supply equal-length
BAM, BAM-index and sample-label arrays.

Each caller publishes:

- `*.riboseq_support_by_bam.tsv`
- `*.riboseq_allele_summary.tsv`
- `*.riboseq_candidate_summary.tsv`
- `*.riboseq_run_summary.tsv`
- `*.riboseq_orf_support_by_bam.tsv`
- `*.riboseq_orf_candidate_summary.tsv`

Evidence levels distinguish SNV ALT support (level 1), variant-coordinate
coverage (level 2), and aligned-base coverage in the same ORF's noncanonical
blocks (level 3). Level 3 requires matching BAM/ORF strand; there is no
P-site or frame check. Coordinate coverage is not proof of ORF translation.

Insertions use the VCF anchor. Deletions use any aligned-base overlap across
`POS` through `POS + len(REF) - 1`; full-span coverage is not required.
REF/ALT assessment applies only to SNVs.

Default MAPQ and base-quality thresholds are zero; secondary, supplementary,
QC-fail and duplicate reads are retained. Stricter filters are available as
WDL inputs. Ribo evidence does not change pVACseq filtering or pVACview tiers.

Build and test instructions:
[0.3.1](pvacnc_riboseq_031.md),
[0.3.2](pvacnc_riboseq_032.md),
[standalone 0.3.3 correction](pvacnc_riboseq_033.md).
The standalone 0.3.3 correction does not change WDL image defaults.

## Non-human runs

`definitions/nonhuman_immuno.wdl` does not yet call pVACnc. Use the standalone
workflow with matched non-human references and species-supported predictors.

The EO771 configuration uses GRCm39/Ensembl 113 and `H-2-Db`, `H-2-Kb`,
`H2-IAb`. Its completed upstream run did not publish a phased VCF.

Do not use `prediction_algorithms: [all]` for the tested pVACtools 7.0.1
mouse configuration: PRIME failed on H-2 alleles. The EO771 configuration uses:

- Class I: MHCflurryEL, MHCflurry, NetMHC, NetMHCcons, NetMHCpan,
  NetMHCpanEL, PickPocket, SMM, SMMPMBEC.
- Class II: NetMHCIIpan, NetMHCIIpanEL, NNalign, SMMalign.

Check contig consistency, species-specific transcript identifiers, VEP fields,
expression/read counts, optional phasing and empty-result behavior before
integrating the non-human parent workflow.
# Official same-noncanonical-block support (2026-10-16)

WDL interface `0.3.5` uses `jinglunli/pvacnc:0.3.5` for preparation,
pVACview annotation and Ribo support. The pVACtools image is unchanged. Level 3
requires coverage in the residual noncanonical block containing the full REF
span, not any block in the ORF. Existing reference_annotation is forwarded
internally to both callers for CDS+stop_codon masking. Coordinate/ALT and read
quality settings are unchanged. See [configuration, contract and tests](pvacnc_riboseq_035.md).
