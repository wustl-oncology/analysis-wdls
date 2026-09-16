# Workflows for Genomic Analysis

A collection of WDL workflows for analysis of genomic sequencing data. These enable primary analysis of many common experiments, including somatic and germline variant calling, RNA sequencing, epigenomic assays, and ommunogenomics approaches for cancer vaccine development.

## Optional pVACnc branch

`definitions/immuno.wdl` includes an opt-in `pVACnc` branch for
reference-noncanonical variant peptides. Set `enable_pvacnc: true` to run it.
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

The ORFanage default is pinned to
`quay.io/biocontainers/orfanage:1.2.0--heaafb18_2`. TransDecoder is pinned to
`quay.io/biocontainers/transdecoder:6.0.0--pl5321hdfd78af_0`. The preparation
image defaults to `jinglunli/pvacnc:0.2.2`. Images and the TransDecoder
minimum ORF length can be overridden per YAML input.

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

The planned preparation image is `jinglunli/pvacnc:0.2.2`. Its entry-point
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

#### TODO: remove the canonical pVACseq-output dependency

The current implementation still passes `pvacseq.combined` into pVACnc and
uses the canonical `Combined.all_epitopes.tsv` to label alleles as
canonical-pVACseq candidates or non-candidates. This is a temporary,
incorrect dependency. The intended implementation must classify eligibility
directly from VEP CSQ consequences in the expression-annotated VCF, then use
that VCF without requiring canonical pVACseq prediction outputs. This also
requires separating VCF expression/read-count annotation from the monolithic
canonical pVACseq prediction subworkflow, so pVACnc can run independently.

Future pVACnc code, image, input-contract, and runtime changes must be
recorded in this README at the same time as the corresponding implementation
change.

#### Standalone module-boundary integration test

The pVACnc workflow can be tested without rerunning the complete immunoNX
workflow by supplying the same files that its upstream calls would normally
produce: the final tumour RNA BAM, reference-only (`-e`) StringTie GTF,
expression-annotated VEP VCF, canonical pVACseq combined report, phased VCF,
consensus HLA alleles, and matched reference files. The workflow then starts at
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

The corrected `0.2.2` module outputs are written to:

```text
/storage1/fs1/mgriffit/Active/natera/human/short_read_analysis/immuno/immuno_runs/NTR014_pvacnc_test/pvacnc_module_outputs_022
```

Its retained Cromwell execution directory is
`/scratch1/fs1/mgriffit/allen/NTR014_pvacnc_module_022`. These paths are
separate from the earlier `0.2.0` and `0.2.1` tests, whose outputs are not
overwritten.

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
