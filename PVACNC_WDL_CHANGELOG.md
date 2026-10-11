# pVACnc WDL changelog

This file tracks pVACnc workflow-interface changes. Container implementation
changes are tracked in the preparation-image `CHANGELOG.md`.

| Version | WDL changes |
|---|---|
| `0.3.4` | Ribo task defaults to image 0.3.2, which requires region reads to match the ORF alignment strand. Adds `riboseq_docker` in the standalone module and `pvacnc_riboseq_docker` in human immuno so both caller Ribo tasks use this version independently of preparation. Preparation defaults and all coordinate/ALT/prediction rules are unchanged. Image build/runtime validation pending. |
| `0.3.3` | Completes the partial 0.3.2 naming fix: both preparation lanes apply `<sample>.<caller>.*` prefixes to generated files after preparation, with matching VCF/index/mapping/routing/block/skip-report names and caller-prefixed annotation logs. Prediction input logs now publish in matching caller-specific arrays/folders rather than colliding scalar filenames. Human parent consumes both log arrays. Standard pVACseq/pVACview filenames remain inside caller-separated report directories. Docker stays 0.3.1; file contents and biological settings are unchanged. |
| `0.3.2` | WDL-only output organization change: standalone ORFanage models now publish as `with_e_orfanage_outputs/<sample>.with_e.orfanage.gtf` and `no_e_orfanage_outputs/<sample>.no_e.orfanage.gtf`, matching the TransDecoder array/folder structure. Replaces the scalar `with_e_orfanage_gtf` and `no_e_orfanage_gtf` workflow outputs with caller-specific output arrays; the human immuno parent consumes the new arrays. Preparation image remains 0.3.1; ORF prediction settings, candidate selection and filtering are unchanged. Existing results are not reorganized. |
| `0.3.1` | Both preparation tasks expose `noncanonical_orf_blocks`; pVACnc passes caller-specific blocks into the Ribo task and publishes two additional ORF-coverage TSVs per caller. The standalone Ribo task requires that block file. New workflow defaults use image 0.3.1. Level 3 is same-ORF noncanonical aligned-base coverage with no P-site/frame requirement. Existing variant-coordinate/SNV checks and prediction/filtering settings are unchanged. Build and cluster validation pending. |
| `0.2.7` | Aligned the ORFanage and TransDecoder mapping contracts for pVACview annotation and retained backward-compatible annotation of older ORFanage mappings. |
| `0.3.0` | Added `definitions/tools/pvacnc_riboseq_support.wdl`. The standalone pVACnc workflow and optional human immuno pVACnc branch accept equal-length Ribo-seq BAM, BAM-index and label arrays and run caller-specific support tasks after pVACseq. Four Ribo-seq TSVs are published under both the ORFanage and TransDecoder bundles. SNVs receive coordinate and REF/ALT assessment; insertions use anchor-coordinate coverage only; deletions use any aligned-base overlap across the VCF REF interval and do not receive REF/ALT classification. Existing pVACnc execution is unchanged unless the new Ribo-seq enable flag is true. All pVACnc preparation/post-processing task defaults now use `jinglunli/pvacnc:0.3.0`. |

## Versioning rule

Every future pVACnc interface, default-image, output-contract or runtime change
must add a row here at the same time as the implementation. Image-internal
changes must also be recorded in the image changelog before publication.
