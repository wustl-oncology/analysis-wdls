# Workflows for Genomic Analysis

A collection of WDL workflows for analysis of genomic sequencing data. These enable primary analysis of many common experiments, including somatic and germline variant calling, RNA sequencing, epigenomic assays, and ommunogenomics approaches for cancer vaccine development.

## Optional pVACnc branch

`definitions/immuno.wdl` includes an opt-in `pVACnc` branch for
reference-noncanonical variant peptides. Set `enable_pvacnc: true` to run it.
Workflow-version changes are recorded in `PVACNC_WDL_CHANGELOG.md`.
The current pVACnc WDL interface version is `0.3.4`; its preparation image
remains `jinglunli/pvacnc:0.3.1`. No new Docker build is required for this
output organization change.
The branch uses the existing WDL RNA BAM, StringTie `-e` GTF,
expression-annotated VCF, HLA calls, phasing VCF and pVACseq parameters. It
adds one shared StringTie no-`-e` assembly and runs two independent ORF
callers on both transcript-model branches:

* ORFanage 1.2.0 in `BEST` mode;
* TransDecoder 6.0.0 with a 30-aa minimum ORF length, strand-specific mode,
  partial and complete ORFs, multiple ORFs per transcript, and no BLAST/Pfam
  retention evidence.

Within each caller, coordinate-level `-e` ORF coverage takes priority and
no-`-e` is used only as rescue. Routing is caller-specific: an ORFanage hit
does not suppress a TransDecoder result, or vice versa. Each caller then has
its own variant projection, derived VCF, pVACseq and pVACview path. The final
output is separated as `pvacnc/orfanage/` and `pvacnc/transdecoder/`.

The branch requires three pinned images before execution:

* `pvacnc_orfanage_docker`: an ORFanage image whose executable accepts the
  command encoded in `definitions/tools/orfanage.wdl`.
* `pvacnc_transdecoder_docker`: the TransDecoder Biocontainers image used by
  `definitions/tools/transdecoder.wdl`.
* `pvacnc_preparation_docker`: the separately maintained pVACnc preparation image, which
  packages the ORFanage and TransDecoder projection, routing, derived-VCF and
  pVACview annotation code.

### Optional Ribo-seq support module

Ribo region support now targets `jinglunli/pvacnc:0.3.2`: only BAM
alignments whose strand matches the predicted noncanonical ORF are counted.
Build this image before running the updated Ribo task. The standalone
`riboseq_docker` and human `pvacnc_riboseq_docker` inputs select the Ribo image
separately from preparation. See [0.3.2 build and rerun instructions](docs/pvacnc_riboseq_032.md).
Existing coordinate and ALT checks retain their rules; old 0.3.1 region
counts include both alignment strands and must be recalculated for this rule.

The current source targets `jinglunli/pvacnc:0.3.1`; build and verify the image
before submitting new workflows. See [0.3.1 instructions](docs/pvacnc_riboseq_031.md).
Both preparation lanes now export strand-aware noncanonical CDS blocks per
mapped ORF. After pVACseq, the optional Ribo module retains the original four
tables and adds `*.riboseq_orf_support_by_bam.tsv` and
`*.riboseq_orf_candidate_summary.tsv`. Level 3 means aligned-base coverage in
the same ORF's noncanonical blocks when neither SNV ALT support (level 1) nor
variant-coordinate coverage (level 2) is present. There is no P-site or frame
check. This evidence level does not change pVACseq filtering or pVACview tiers.

Image and WDL version `0.3.0` add an optional Ribo-seq evidence module after
the ORFanage and TransDecoder pVACseq calls. Enable it with
`enable_pvacnc_riboseq_support: true` in the complete immuno workflow or
`enable_riboseq_support: true` in the standalone pVACnc workflow, and provide
equal-length BAM, BAM-index and sample-label arrays. The module runs
independently for the two ORF callers and publishes four TSVs under each
caller's pVACnc output bundle:

* `*.riboseq_support_by_bam.tsv` — one row per DNA allele and Ribo-seq BAM;
* `*.riboseq_allele_summary.tsv` — one row per unique DNA allele;
* `*.riboseq_candidate_summary.tsv` — one row per pVACnc allele–ORF mapping,
  annotated with combined `filtered.tsv` membership;
* `*.riboseq_run_summary.tsv` — input, threshold and support counts.

Coordinate coverage, REF-read support and ALT-read support are distinct. Only
`has_alt_support=yes` represents mutant-allele Ribo-seq evidence for SNVs only.
Insertions and deletions receive coordinate-level assessment without REF/ALT
classification. Insertions use the VCF anchor coordinate; deletions use `POS`
through `POS + len(REF) - 1`, with any aligned-base overlap counted as support.
Full deletion-span coverage is not required. The default
MAPQ/base-quality thresholds are zero and secondary, supplementary, QC-fail
and duplicate alignments are retained to reproduce the original EO771 check.
All filters can be enabled through explicit WDL inputs. The implementation is
packaged as `pvacnc-riboseq-support` in `jinglunli/pvacnc:0.3.0`; the WDL
remains in this repository at `definitions/tools/pvacnc_riboseq_support.wdl`.

The ORFanage default is pinned to
`quay.io/biocontainers/orfanage:1.2.0--heaafb18_2`. TransDecoder is pinned to
`quay.io/biocontainers/transdecoder:6.0.0--pl5321hdfd78af_0`. The preparation
image defaults to `jinglunli/pvacnc:0.3.1`. Images and the TransDecoder
minimum ORF length can be overridden per YAML input.

Version `0.2.7` retains the `0.2.6` candidate-skip handling, adds the missing
ORFanage `annotation_status` mapping field to match TransDecoder, and allows
pVACview annotation to infer de novo handling from `orf_origin=no_e_rescue`
for older mappings. Version `0.2.6` introduced candidate-specific ORFanage and
TransDecoder ORF/model
loading, projection and VCF-construction failures without terminating the
lane. Every skipped allele--ORF pair is reported in a lane-specific TSV with a
machine-readable reason and details.
Conditional WDL skipping of downstream prediction for a zero candidate count
is not yet implemented and must be validated before empty-result handling is
considered complete.

### pVACnc implementation status and required containers

The pVACnc preparation image is maintained outside this WDL repository. It
provides Python 3, `bgzip`, `tabix`, and the projection/derived-VCF scripts.
Its public task entry points are installed as
`/usr/local/bin/pvacnc-prepare` for ORFanage,
`/usr/local/bin/pvacnc-transdecoder-prepare` for TransDecoder and
`/usr/local/bin/pvacnc-annotate-pvacview` for both callers. The WDL invokes
these commands by name, consistent with other containerized tools in this
repository.

Each preparation task materializes only the canonical
`Combined.all_epitopes.tsv` as a regular task-local file. The image-local
wrappers select that bounded input using a shell filename pattern and do not
recursively search task directories.

The planned preparation image is `jinglunli/pvacnc:0.3.1`. Its entry-point
wrapper resets its internal `PATH` so LSF cannot substitute a host-side Conda
Python for the Python installed in the image. The ORFanage lane distinguishes
reference and assembled-transcript annotation: with-e candidates retain the
exact ENST-matched VEP transcript metadata, whereas no-e candidates use the
sample-specific StringTie transcript and gene IDs. For no-e candidates, all
unique coordinate-level VEP gene symbols are reported in the existing gene
name field with `/` as the delimiter, while MANE Select, Canonical, TSL,
APPRIS, PICK, CCDS and reference-protein identifiers are left unset. Complete
source-VEP contexts remain available in the mapping TSV for provenance.

TransDecoder follows the same reference-provenance policy as ORFanage. A
with-e-priority ORF whose parent StringTie transcript exactly matches a source
VEP ENST is reported under that ENST and retains its VEP gene, MANE Select,
Canonical and TSL fields. Only no-e-rescue ORFs use the predicted TransDecoder
ORF Feature and leave reference-transcript fields unset. The mapping TSV also
retains the parent StringTie transcript, its expression, ORF origin and
complete coordinate-level VEP contexts. With-e candidates keep the original
WDL/Kallisto `TX` and `GX` values from the expression-annotated source VCF in
place. Synthetic TransDecoder CSQ entries retain the pVACseq-compatible
`BIOTYPE=protein_coding`; source reference biotypes do not replace this
workflow-interface field.
both caller lanes; StringTie expression is injected only for no-e candidates.

Frameshift translation is also harmonized between the callers. For ORFanage,
the preparation image reconstructs the complete spliced transcript from the
ORFanage exon model and reference genome, maps the CDS start and variant into
cDNA coordinates, and translates the mutant sequence from the ORF start to the
first new stop codon or transcript end. The builder verifies that the
reconstructed wild-type translation agrees with the original ORFanage CDS
translation before emitting a `FrameshiftSequence`. TransDecoder applies the
same full-parent-cDNA translation boundary.

After pVACseq prediction, a pVACnc-specific annotation task updates the
existing pVACseq/pVACview transcript fields for de novo Features. `MANE
Select`, `Canonical`, and `TSL` remain empty and the pVACview-computed
`Transcript Pass` is shown as `de novo`; no additional report column is
introduced. In tier calculation, only the reference-transcript-quality check
is exempt. Binding, expression, VAF, clonality, anchor, reference-match,
problematic-position and other tier criteria remain unchanged. In the
ORFanage and TransDecoder lanes this applies only to no-e rescue Features.
With-e features retain normal reference-transcript filtering. The annotation
task reads the explicit mapping contract rather than relying on transcript ID
patterns.

#### Canonical pVACseq independence

pVACnc does not consume canonical `Combined.all_epitopes.tsv` output. Both ORF
callers start from all ALT alleles in the expression-annotated VCF, retain
`FILTER=PASS`, exclude alleles whose reference span overlaps any annotated
reference CDS, and then retain protein-changing effects in their predicted
ORFs. The integrated immuno workflow currently obtains the shared
expression-annotated VCF from the canonical pVACseq subworkflow; separating
that upstream annotation task remains necessary before pVACnc can be scheduled
entirely independently of the canonical subworkflow.

Future pVACnc code, image, input-contract, and runtime changes must be
recorded in this README at the same time as the corresponding implementation
change.

#### Standalone module-boundary integration test

Starting with WDL interface version `0.3.2`, standalone ORF model files use
the same caller/branch folder structure:

```text
with_e_orfanage_outputs/<sample>.with_e.orfanage.gtf
no_e_orfanage_outputs/<sample>.no_e.orfanage.gtf
with_e_transdecoder_outputs/<sample>.with_e.transdecoder.genome.gff3
no_e_transdecoder_outputs/<sample>.no_e.transdecoder.genome.gff3
```

TransDecoder folders also contain cDNA, CDS, peptide, transcript-GFF3 and
summary files. ORFanage emits one GTF per branch. The integrated human
workflow includes each caller's models in its existing `pvacnc/orfanage/`
or `pvacnc/transdecoder/` bundle. Previous published results retain their
original layout. This change replaces the standalone scalar outputs
`with_e_orfanage_gtf` and `no_e_orfanage_gtf` with
`with_e_orfanage_outputs` and `no_e_orfanage_outputs` arrays.

Version `0.3.3` completes caller-prefix symmetry across preparation outputs:
ORFanage uses `<sample>.orfanage.*`, and TransDecoder uses
`<sample>.transdecoder.*`. This includes input VCF/index, mapping, routing,
effects, noncanonical blocks, skip reports and preparation audit files.
The WDL applies these publication names after each image command completes;
file contents and identifiers are unchanged. Both annotation logs now have
the same caller prefixes, and prediction input YAML logs publish under
`orfanage_prediction_logs/` and `transdecoder_prediction_logs/`.
Standard pVACseq report and pVACview asset names remain inside the symmetric
`orfanage_mhc_i/`, `orfanage_mhc_ii/`, `orfanage_combined/` and corresponding
`transdecoder_*` directories so pVACview's expected file names remain valid.
Existing published runs retain their prior names. This WDL-only release
continues to use preparation image `0.3.1`.

The pVACnc workflow can be tested without rerunning the complete immunoNX
workflow by supplying the same files that its upstream calls would normally
produce: the final tumour RNA BAM, reference-only (`-e`) StringTie GTF,
expression-annotated VEP VCF, optional phased VCF, MHC alleles, and matched
reference files. The workflow then starts at
the no-`-e` StringTie call and runs both the ORFanage and TransDecoder lanes
through pVACseq and pVACview annotation.

The NTR014 module test is configured in:

```text
/storage1/fs1/mgriffit/Active/natera/human/short_read_analysis/immuno/immuno_runs/NTR014_pvacnc_test/yamls/NTR014_pvacnc_module.yaml
```

Submit it on compute1 with:

```bash
bash /storage1/fs1/mgriffit/Active/natera/human/short_read_analysis/immuno/immuno_runs/NTR014_pvacnc_test/run_pvacnc_module_compute1.sh
```

The `0.2.3` module outputs are written to:

```text
/storage1/fs1/mgriffit/Active/natera/human/short_read_analysis/immuno/immuno_runs/NTR014_pvacnc_test/pvacnc_module_outputs_023
```

Its retained Cromwell execution directory is
`/scratch1/fs1/mgriffit/allen/NTR014_pvacnc_module_023`. These paths are
separate from the earlier tests, whose outputs are not overwritten.

#### Non-human validation status

`definitions/nonhuman_immuno.wdl` does not yet call pVACnc. Before adding that
integration, the unchanged standalone `definitions/pvacnc.wdl` is being
validated against the completed EO771 non-human immunoNX outputs. The test
configuration and submission wrapper are isolated from the historical mouse
projection analyses at:

```text
/storage1/fs1/mgriffit/Active/griffithlab/adhoc/allen/MS_TA_work/Noncanonical_ORF_Variant_Projection_mouse/current_wdl_runs/EO771_unstimulated_pvacnc_025
```

This validation uses the published EO771 tumour RNA BAM, reference-only
StringTie GTF and expression-annotated VEP VCF together with the exact
GRCm39/Ensembl 113 genome, annotation and peptide reference. It uses
`H-2-Db`, `H-2-Kb` and `H2-IAb` and intentionally omits the optional phased
VCF because the completed non-human run did not publish one. The current WDL
remains authoritative: StringTie without `-e` is run from the BAM, ORFanage
uses `BEST`, both ORFanage and TransDecoder use `-e` priority with no-`-e`
rescue, all annotated reference-CDS overlaps are excluded, and preparation
uses `jinglunli/pvacnc:0.2.7`. Historical mouse aggregated-maximal GTFs,
ORFanage-only scripts, canonical pVACseq-output dependency and pVACtools 6.0.2
wrapper are not reused.

Potential non-human failure modes that must be audited before integration and
pull request submission are:

* mouse `ENSMUST`, `ENSMUSG` and `ENSMUSP` identifiers must retain correct
  reference provenance without any implicit human `ENST`/`ENSG` assumption;
* VCF, BAM, FASTA and GTF contig names must agree exactly;
* the mouse VEP CSQ schema must provide every field required by derived-VCF
  construction, while absent human-centric MANE or TSL values remain safe;
* H-2 class-I and class-II alleles must be routed to supported prediction
  algorithms and pVACtools must identify the run as mouse;
* missing phased-proximal-variant input must remain a valid optional case;
* expression and read-count annotations must survive into derived VCFs and
  final reports under mouse transcript identifiers;
* ORFanage and TransDecoder must produce valid results against the exact
  GRCm39/Ensembl 113 reference pair; and
* a lane with zero candidates must finish cleanly rather than invoking
  pVACseq on an unusable empty derived VCF.

Only after the standalone EO771 run and these audits succeed should pVACnc be
added to `definitions/nonhuman_immuno.wdl` as an opt-in peer of pVACseq.

The first EO771 standalone attempt exposed a pVACtools 7.0.1 compatibility
failure when `prediction_algorithms: [all]` selected PRIME for mouse H-2
alleles. The rerun configuration therefore uses only predictors represented
in the successful canonical EO771 reports. Class I uses MHCflurryEL,
MHCflurry, NetMHC, NetMHCcons, NetMHCpan, NetMHCpanEL, PickPocket, SMM and
SMMPMBEC. Class II uses NetMHCIIpan, NetMHCIIpanEL, NNalign and SMMalign.

These pipelines were developed by the [Washington University Division of Oncology](https://oncology.wustl.edu/), in collaboration with the [McDonnell Genome Institute](https://genome.wustl.edu/) and other members of the [Washington University School of Medicine](https://medicine.wustl.edu/).  Much of the code and structure is adapted from a sister repository of [workflows in the CWL language](https://github.com/genome/analysis-workflows/).


## Getting Started

### Workflows
Download our repository with `git clone https://github.com/wustl-oncology/analysis-wdls.git`.

Terra hosts an excellent guide to [getting started with WDL](https://support.terra.bio/hc/en-us/articles/360037117492-Overview-Getting-started-with-WDL).

### Workflow Execution
These are primarily designed to be run with the cromwell workflow system, either on a self-hosted server through the Google Cloud Platform ([details here](https://github.com/wustl-oncology/cloud-workflows)), through [Terra](https://terra.bio/)), or on local compute clusters. They have been extensively used in these contexts, though the portability of WDL means that running on other platforms should be straightforward as well.

### Docker
In order to provide a portable environment, each tool in our workflow has a designated Docker container. [Download Docker here](https://www.docker.com/products/docker-desktop).

### Data
Input "bundles" containing reference genomes and annotations are available for use/download along with documentation on how they were created.

### Help
We answer questions through github issues on this repository, and also have compiled a list of [common errors](https://github.com/wustl-oncology/analysis-wdls/blob/main/docs/common_errors.md) that may be useful.

## Contributions

A big thanks to all of the developers, bioinformaticians, and scientists who built this resource. For a complete list of software contributions, i.e. commits, to this repository, please see the GitHub Contributors both to [this repository](https://github.com/wustl-oncology/analysis-wdls/graphs/contributors) as well as to the [analysis-workflows](https://github.com/genome/analysis-workflows/graphs/contributors) repo.

Genomic data evolves rapidly and pull requests and improvements that add new features are welcome!
### EO771 0.2.6 failure-boundary recovery

The EO771 standalone 0.2.6 execution completed StringTie, both ORFanage calls,
both TransDecoder calls, pVACnc ORFanage preparation, and ORFanage pVACseq
prediction before failing in `orfanageAnnotatePvacview` because the 0.2.6
mapping lacked `annotation_status`. A case-local recovery chain reuses those
completed artifacts, annotates the existing ORFanage prediction with 0.2.7,
and continues the TransDecoder lane from pVACnc preparation. It does not rerun
StringTie, ORFanage, TransDecoder, or ORFanage pVACseq. Recovery outputs are
written separately under `pvacnc_module_outputs_026_resume_027`; the original
0.2.6 execution remains unchanged.
### Standalone Ribo-seq 0.3.3 correction

Source and 12-set test commands are documented in
`docs/pvacnc_riboseq_033.md`. Selected-variant support is restricted to exact
selected ORFs; all evidence levels require matching BAM alignment strand;
Level 3 excludes reference annotated CDS on either strand. This standalone
test does not change existing WDL image defaults or overwrite prior results.
Image build and production validation remain pending.
