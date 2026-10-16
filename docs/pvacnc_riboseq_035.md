# pVACnc 0.3.5

Preparation, pVACview annotation and Ribo-seq tasks default to
`jinglunli/pvacnc:0.3.5`. ORFanage, TransDecoder and pVACtools use their
existing images. Image build and runtime validation are required before use.

## Candidate routing

Both callers retain supported protein-changing events from no-e transcripts,
even when a with-e ORF overlaps the allele. A successful with-e representation
takes priority only when allele, strand and complete WT/MT protein sequences
are identical. Mapping and skip tables retain equivalent-ORF provenance.

## Ribo-seq evidence

Level 3 requires same-strand aligned-base coverage in the residual noncanonical
block containing the variant's full REF span (anchor for insertions). Evidence
from other blocks in the ORF is not used. Reference CDS and stop-codon features
on either strand are excluded from these blocks.

Coordinate/ALT checks, read-quality defaults and pVACseq filters are unchanged.
Regenerate block and support tables when rerunning older results; do not reuse
counts calculated under earlier evidence rules.

## Image overrides

Standalone inputs are `pvacnc.preparation_docker` and `pvacnc.riboseq_docker`.
Integrated inputs are `immuno.pvacnc_preparation_docker` and
`immuno.pvacnc_riboseq_docker`. Existing YAML overrides take precedence over
the WDL defaults; update them explicitly when starting a new run.
