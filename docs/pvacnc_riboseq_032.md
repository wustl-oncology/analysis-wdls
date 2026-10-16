# Ribo region support: image 0.3.2

Only reads whose BAM alignment strand matches the noncanonical ORF strand
contribute to region coverage. SAM flag 16 means minus strand; an unset flag
means plus strand. This implements the user-selected alignment-strand rule,
without library-orientation conversion, P-site offsets or frame inference.
Both ORFanage and TransDecoder use the same checker. Coordinate and ALT
tables are unchanged. Plot region bars remain cumulative: region evidence
union coordinate evidence union ALT evidence.

WDL interface 0.3.4 uses image 0.3.2 for Ribo tasks. The standalone workflow
input is `riboseq_docker`; human immuno uses `pvacnc_riboseq_docker`.
Preparation remains on its existing image version.

Build on the cluster:

```bash
cd /storage1/fs1/mgriffit/Active/griffithlab/adhoc/allen/MS_TA_work/Noncanonical_ORF_Variant_Projection_wdl/pvacnc_container/pvacnc
bsub -G compute-oncology -g /jinglun/wdl -q oncology-interactive -Is \
  -a 'docker_build(jinglunli/pvacnc)' -- --tag jinglunli/pvacnc:0.3.2 .
```

Smoke test:

```bash
bsub -G compute-oncology -g /jinglun/wdl -q general-interactive -Is \
  -n 1 -M 4GB -R 'rusage[mem=4GB]' \
  -a 'docker(jinglunli/pvacnc:0.3.2)' -- \
  /usr/local/bin/python3 -B /opt/pvacnc/test_riboseq_031.py
```

Recheck the four 100-aa filtered datasets in new timestamped folders:

```bash
cd /storage1/fs1/mgriffit/Active/griffithlab/adhoc/allen/MS_TA_work/Noncanonical_ORF_Variant_Projection_mouse
bash submit_riboseq_100aa_032_same_strand.sh
```

Eight jobs reuse the existing 100-aa candidate summaries and noncanonical
block tables from `riboseq_100aa_031_filtered_20261005_113458`. Only region
coverage is recomputed: EO771 six BAMs, EO771 RAPPL, EO771 seven combined,
KP2 RAPPL; two callers each. No upstream workflow is rerun. Logs live inside
each new timestamped root under `logs/`. The submitter prints the exact path
manifest and the plotting command to run after jobs finish. Source results
and earlier figures are preserved.

Local synthetic tests pass for matching plus/minus strands, opposite-strand
exclusion, CIGAR skips/deletions, flag filters and evidence-tier precedence.
Image publication, digest and production results remain UNCONFIRMED until
the user completes the cluster build and rerun.

## Complete 12-set rerun

Use `submit_all_riboseq_032_same_strand.sh` in the mouse project to cover
all three collections: original filtered, original pVACview Pass, and
100-aa TransDecoder filtered. Each contains EO771 six BAMs, EO771 RAPPL,
EO771 seven combined and KP2 RAPPL. It submits 14 caller-specific region
checks: six for original evidence and eight for 100-aa evidence. Original
filtered/Pass selections reuse the same recalculated support; their seven-BAM
sets are unions of six-BAM and RAPPL support. No support count is summed as
a variant count. Use the printed `plot_all_riboseq_032_same_strand.sh` command
after completion to create all 12 folders, each with PNG/PDF/SVG and two CSVs.
Everything is under a new `riboseq_all_032_same_strand_TIMESTAMP` root;
earlier results remain preserved. Shell syntax checks passed; jobs have not
been submitted by the assistant.
