# Workflows for Genomic Analysis

A collection of WDL workflows for analysis of genomic sequencing data. These enable primary analysis of many common experiments, including somatic and germline variant calling, RNA sequencing, epigenomic assays, and ommunogenomics approaches for cancer vaccine development.

## pVACnc

Optional noncanonical variant-peptide analysis in `definitions/immuno.wdl`.
Enable with `enable_pvacnc: true`, or run `definitions/pvacnc.wdl` standalone.

- WDL interface: `0.3.4`.
- Preparation image: `jinglunli/pvacnc:0.3.1`.
- Ribo-seq image: `jinglunli/pvacnc:0.3.2` (configured separately).
- ORF callers: ORFanage 1.2.0 (`BEST`) and TransDecoder 6.0.0
  (30-aa minimum by default).

Both callers use StringTie `-e` models first, with no-`-e` assembly as
coordinate-level rescue. Each has separate variant projection, VCF,
pVACseq and pVACview outputs under `pvacnc/orfanage/` and
`pvacnc/transdecoder/`.

Inputs include the tumor RNA BAM, StringTie `-e` GTF, expression-annotated
VEP VCF, MHC alleles and matched references; phased VCF is optional.
Preparation retains `FILTER=PASS` alleles outside all annotated reference
CDSs. Canonical pVACseq prediction reports are not required.

Optional Ribo-seq support produces caller-specific evidence tables without
changing candidate filtering or pVACview tiers.

### Protein-event rescue: image 0.3.4 source prepared

EO771 WDL034 30-aa example: `11:54757309 G>A` is synonymous in the
selected with-`-e` ORF `ENSMUST00000117710.p1`, but missense **G34S**
(`GGC` to `AGC`, CDS position 100) in no-`-e` ORF
`EO771_unstimulated.no_e.4747.1.p3`. Current coordinate-only priority
suppresses the no-`-e` result and omits this candidate before pVACseq.
Nearby `11:54757279 T>C` is retained as C24R in that same no-`-e` ORF.

- Docker source now admits supported protein-changing no-`-e` events even with
  with-`-e` overlap. Exact final WT/MT sequence equivalents for the same allele
  and strand prefer with-`-e`; both origins remain traceable in mapping/skip logs.
- Build/publish/smoke-test image `jinglunli/pvacnc:0.3.4`, then explicitly select
  it for NEW runs and audit candidate differences. WDL defaults are unchanged.
- Nearby C24R/G34S variants are still applied individually, not phased or combined.
- The same pending image also masks reference `CDS + stop_codon` for Level-3
  block generation/reclipping; variant eligibility is unchanged. Regenerate
  blocks and region-support counts rather than reuse old results.
- Source handoff: `MS_TA_work/Noncanonical_ORF_Variant_Projection_wdl/pvacnc_container/pvacnc/ROUTING_034_HANDOFF.md`
  under the shared allen project. Local regressions pass; production validation pending.

See [configuration and outputs](docs/pvacnc.md) and the
[WDL changelog](PVACNC_WDL_CHANGELOG.md). Record interface, image and runtime
changes in both documents.

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
