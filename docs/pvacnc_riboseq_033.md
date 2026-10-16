# Ribo-seq 0.3.3: selected ORF and strict noncanonical coordinates

All levels require BAM alignment strand matching the predicted ORF. Level 3
blocks exclude reference annotated CDS on either strand. No frame/P-site test.
MAPQ/base-quality defaults and alignment inclusion policy are unchanged.

The 12 standalone jobs cover three groups (original filtered, pVACview Pass,
TransDecoder-100aa filtered) × four datasets (EO771 six BAMs, RAPPL sorted,
seven combined; KP2 RAPPL sorted). Each job runs both callers. Only exact
selected allele–ORF matches contribute to selected variant support. Filtered
events use transcript/HGVSp/REF/ALT and uniquely resolved POS. Pass sets use
the displayed representative Index from the aggregated Pass row, recovered
in all_epitopes and matched to its ORF. This does not claim support for every
transcript attached to a Pass allele. Join errors are explicit, not silently
interpreted as no support.

Old candidate summaries supply mapping provenance only. All read counts are
recomputed. Old noncanonical blocks are further clipped against the union of
reference CDS on both strands; this is equivalent to subtracting that union
from the original ORF CDS. Existing blocks/models and pVACseq results are not
overwritten. New model block tables and selected-event audit TSVs are saved.
The full-mapping denominator is retained for Unique variants; only matched
selected ORFs contribute to levels 1–3. Region bars remain cumulative.

## Build, test, submit

```bash
cd /storage1/fs1/mgriffit/Active/griffithlab/adhoc/allen/MS_TA_work/Noncanonical_ORF_Variant_Projection_wdl/pvacnc_container/pvacnc
bsub -G compute-oncology -g /jinglun/wdl -q oncology-interactive -Is \
  -a 'docker_build(jinglunli/pvacnc)' -- --tag jinglunli/pvacnc:0.3.3 .

bsub -G compute-oncology -g /jinglun/wdl -q general-interactive -Is \
  -n 1 -M 4GB -R 'rusage[mem=4GB]' -a 'docker(jinglunli/pvacnc:0.3.3)' -- \
  /usr/local/bin/python3 -B /opt/pvacnc/test_riboseq_033.py

cd /storage1/fs1/mgriffit/Active/griffithlab/adhoc/allen/MS_TA_work/Noncanonical_ORF_Variant_Projection_mouse
bash submit_all_riboseq_033.sh
```

The submitter prints the exact new output root and plotting command. After all
12 jobs succeed, run `bash plot_all_riboseq_033.sh OUTPUT_ROOT` with a Python
environment containing matplotlib (override PYTHON_BIN if needed). Logs are
under OUTPUT_ROOT/logs; support tables under support/GROUP/DATASET/CALLER;
PNG/PDF/SVG and membership/count tables under figures/GROUP/DATASET.

## Verification

Ten mock tests passed, including opposite-strand exclusion, both callers,
exact ORF selection, and indel coverage. Small real filtered-report joins
matched 17/122 and 17/70 EO771 ORFs, 4/15 and 4/8 KP2 ORFs (original/100aa).
No production BAM scans performed. Image publication/digest, Pass Index joins,
production counts and full WDL runtime remain UNCONFIRMED.

Existing WDL image defaults are not changed by this standalone 12-job test.
The new jobs explicitly use 0.3.3 and strictly reclip old blocks. Future WDL
adoption must update both block production and Ribo counting; changing only
the Ribo image does not regenerate old blocks.
