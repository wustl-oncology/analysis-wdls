# pVACnc WDL changelog

This file tracks pVACnc workflow-interface changes. Container implementation
changes are tracked in the preparation-image `CHANGELOG.md`.

| Version | WDL changes |
|---|---|
| `0.3.5` | All pVACnc preparation, annotation and Ribo image defaults use 0.3.5. Both callers retain distinct no-e protein events. Ribo support uses the variant-containing noncanonical block and reference CDS/stop-codon masking. No new top-level inputs; pVACtools image unchanged. Build and runtime validation pending. |
| `0.3.4` | Separate Ribo image inputs (`riboseq_docker`, `pvacnc_riboseq_docker`), defaulting to 0.3.2. Region support requires matching BAM/ORF strand. Preparation and coordinate/ALT checks unchanged. Runtime validation pending. |
| `0.3.3` | Caller-prefixed preparation files and annotation logs; separate prediction-log arrays. Standard report names retained within caller directories. Preparation image remains 0.3.1; biological settings unchanged. |
| `0.3.2` | ORFanage model outputs now use caller/branch folders matching TransDecoder. Replaces scalar GTF outputs with `with_e_orfanage_outputs` and `no_e_orfanage_outputs` arrays. Preparation image remains 0.3.1; existing results unchanged. |
| `0.3.1` | Adds required `noncanonical_orf_blocks` input and two ORF-coverage Ribo tables per caller. Level 3 uses aligned-base coverage without P-site/frame checks. Image defaults: 0.3.1. Cluster validation pending. |
| `0.3.0` | Optional caller-specific Ribo support after pVACseq, with BAM/index/label arrays and four TSVs. SNVs receive REF/ALT assessment; indels receive coordinate assessment. Preparation/post-processing image defaults: 0.3.0. |
| `0.2.7` | Matching caller mapping contracts for pVACview, with support for older ORFanage mappings. |

## Versioning rule

Record interface, image-default, output-contract and runtime changes here
with the implementation. Record image-internal changes in the image changelog.
