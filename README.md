# Next Major Steps

Our main concern so far has been getting CWLs converted to
WDL. Following this will be efforts on optimization of the workflows,
and cleanup of the repository.

In the future we may rework the structure of this repository to
a format that Dockstore supports and leverage that tool.


# Differences from CWL

## Directory types must be a zip file, or Array[File]

There is not yet a supported Directory type in WDL. Instances of this
like `Directory vep_cache_dir` which involve nested directory structure are
replaced with `File vep_cache_dir_zip`. Instances of this like
`Directory hla_call_files` which are just a flat collection of files are
replaced with `Array[File] hla_call_files`.


## Input files must prefix arguments with the name of the workflow

Input files must prefix each argument with the name of the workflow
they're going to run, because a WDL file can contain multiple
workflows or pass inputs over a layer if they aren't propagated
through in the definition. e.g. to call workflow `somaticExome` with
input `foo`, yaml key must be `somaticExome.foo`

If WDLs are being used leveraging the
[`cloud-workflows/scripts/cloudize-workflow.py` helper
script](https://github.com/griffithlab/cloud-workflows/tree/main/scripts),
the generated input file will have this handled already.


# Conversions
## Pipelines

- [x] alignment\_exome
- [x] alignment\_exome\_nonhuman
- [ ] alignment\_umi\_duplex
- [ ] alignment\_umi\_molecular
- [x] alignment\_wgs
- [x] alignment\_wgs\_nonhuman
- [ ] aml\_trio\_cle
- [ ] aml\_trio\_cle\_gathered
- [x] bisulfite
- [ ] chipseq
- [ ] chipseq\_alignment\_nonhuman
- [x] detect\_variants
- [x] detect\_variants\_nonhuman
- [x] detect\_variants\_wgs
- [ ] downsample\_and\_recall
- [ ] gathered\_downsample\_and\_recall
- [x] germline\_exome
- [x] germline\_exome\_gvcf
- [x] germline\_exome\_hla\_typing
- [x] germline\_wgs
- [x] germline\_wgs\_gvcf
- [x] immuno
- [x] rnaseq
- [x] rnaseq\_star\_fusion
- [x] rnaseq\_star\_fusion\_with\_xenosplit
- [x] somatic\_exome
- [x] somatic\_exome\_cle
- [ ] somatic\_exome\_cle\_gathered  # This doesn't make sense in cloud
- [ ] somatic\_exome\_gathered       # This doesn't make sense in cloud
- [x] somatic\_exome\_nonhuman
- [x] somatic\_wgs
- [x] tumor\_only\_detect\_variants
- [x] tumor\_only\_exome
- [x] tumor\_only\_wgs


## Subworkflows

- [ ] align
- [ ] align\_sort\_markdup
- [x] bam\_readcount
- [x] bam\_to\_trimmed\_fastq
- [x] bam\_to\_trimmed\_fastq\_and\_biscuit\_alignments
- [x] bam\_to\_trimmed\_fastq\_and\_hisat\_alignments
- [x] bgzip\_and\_index
- [x] bisulfite\_qc
- [ ] cellranger\_mkfastq\_and\_count
- [x] cnvkit\_single\_sample
- [ ] cram\_to\_bam\_and\_index
- [ ] cram\_to\_cnvkit
- [x] docm\_cle
- [x] docm\_germline
- [ ] duplex\_alignment
- [x] filter\_vcf
- [x] filter\_vcf\_nonhuman
- [x] fp\_filter
- [x] gatk\_haplotypecaller\_iterator
- [x] germline\_detect\_variants
- [x] germline\_filter\_vcf
- [x] hs\_metrics
- [ ] joint\_genotype
- [x] merge\_svs
- [ ] molecular\_alignment
- [ ] molecular\_qc
- [x] mutect
- [x] phase\_vcf
- [x] pindel
- [x] pindel\_cat
- [ ] pindel\_region
- [x] pvacseq
- [x] qc\_exome
- [x] qc\_exome\_no\_verify\_bam
- [x] qc\_wgs
- [x] qc\_wgs\_nonhuman
- [ ] sequence\_align\_and\_tag\_adapter
- [x] sequence\_to\_bqsr
- [x] sequence\_to\_bqsr\_nonhuman
- [ ] single\_cell\_rnaseq
- [x] single\_sample\_sv\_callers
- [x] strelka\_and\_post\_processing
- [x] strelka\_process\_vcf
- [x] sv\_depth\_caller\_filter
- [x] sv\_paired\_read\_caller\_filter
- [ ] umi\_alignment
- [x] varscan
- [x] varscan\_germline
- [x] varscan\_pre\_and\_post\_processing
- [ ] vcf\_eval\_cle\_gold
- [ ] vcf\_eval\_concordance
- [x] vcf\_readcount\_annotator


## Tools

- [x] add\_strelka\_gt
- [x] add\_string\_at\_line
- [x] add\_string\_at\_line\_bgzipped
- [x] add\_vep\_fields\_to\_table
- [ ] align\_and\_tag
- [x] annotsv
- [x] annotsv\_filter
- [x] apply\_bqsr
- [x] bam\_readcount
- [x] bam\_to\_bigwig
- [x] bam\_to\_cram
- [x] bam\_to\_fastq
- [ ] bam\_to\_sam
- [x] bcftools\_merge
- [x] bedgraph\_to\_bigwig
- [ ] bedtools\_intersect
- [x] bgzip
- [x] biscuit\_align
- [x] biscuit\_markdup
- [x] biscuit\_pileup
- [x] bisulfite\_qc\_conversion
- [x] bisulfite\_qc\_coverage\_stats
- [x] bisulfite\_qc\_cpg\_retention\_distribution
- [x] bisulfite\_qc\_mapping\_summary
- [x] bisulfite\_vcf2bed
- [x] bqsr
- [ ] call\_duplex\_consensus
- [ ] call\_molecular\_consensus
- [x] cat\_all
- [x] cat\_out
- [ ] cellmatch\_lineage
- [ ] cellranger\_atac\_count
- [ ] cellranger\_count
- [ ] cellranger\_feature\_barcoding
- [ ] cellranger\_mkfastq
- [ ] cellranger\_vdj
- [ ] cle\_aml\_trio\_report\_alignment\_stat
- [ ] cle\_aml\_trio\_report\_coverage\_stat
- [ ] cle\_aml\_trio\_report\_full\_variants
- [ ] clip\_overlap
- [x] cnvkit\_batch
- [x] cnvkit\_vcf\_export
- [x] cnvnator
- [x] collect\_alignment\_summary\_metrics
- [x] collect\_gc\_bias\_metrics
- [x] collect\_hs\_metrics
- [x] collect\_insert\_size\_metrics
- [x] collect\_wgs\_metrics
- [ ] combine\_gvcfs
- [x] combine\_variants
- [ ] combine\_variants\_concordance
- [x] combine\_variants\_wgs
- [x] concordance
- [ ] cram\_to\_bam
- [x] docm\_add\_variants
- [x] docm\_gatk\_haplotype\_caller
- [ ] downsample
- [x] duphold
- [ ] duplex\_seq\_metrics
- [ ] eval\_cle\_gold
- [ ] eval\_vaf\_report
- [x] extract\_hla\_alleles
- [ ] extract\_umis
- [ ] fastq\_to\_bam
- [ ] filter\_consensus
- [x] filter\_known\_variants
- [x] filter\_sv\_vcf\_blocklist\_bedpe
- [x] filter\_sv\_vcf\_depth
- [x] filter\_sv\_vcf\_read\_support
- [x] filter\_sv\_vcf\_size
- [x] filter\_vcf\_cle
- [x] filter\_vcf\_coding\_variant
- [x] filter\_vcf\_custom\_allele\_freq
- [x] filter\_vcf\_depth
- [x] filter\_vcf\_docm
- [x] filter\_vcf\_mapq0
- [x] filter\_vcf\_somatic\_llr
- [ ] fix\_vcf\_header
- [x] fp\_filter
- [ ] gather\_to\_sub\_directory
- [ ] gatherer
- [ ] gatk\_genotypegvcfs
- [x] gatk\_haplotype\_caller
- [x] generate\_qc\_metrics
- [x] germline\_combine\_variants
- [ ] grolar
- [ ] group\_reads
- [x] hisat2\_align
- [x] hla\_consensus
- [ ] homer\_tag\_directory
- [x] index\_bam
- [x] index\_cram
- [x] index\_vcf
- [x] intersect\_known\_variants
- [x] interval\_list\_expand
- [x] intervals\_to\_bed
- [x] kallisto
- [ ] kmer\_size\_from\_index
- [x] manta\_somatic
- [x] mark\_duplicates\_and\_sort
- [ ] mark\_illumina\_adapters
- [x] merge\_bams
- [ ] merge\_bams\_samtools
- [x] merge\_vcf
- [x] mutect
- [x] name\_sort
- [x] normalize\_variants
- [x] optitype\_dna
- [x] picard\_merge\_vcfs
- [x] pindel
- [ ] pindel2vcf
- [x] pindel\_somatic\_filter
- [ ] pizzly
- [ ] pvacbind
- [ ] pvacfuse
- [x] pvacseq
- [x] pvacseq\_combine\_variants
- [ ] pvacvector
- [ ] read\_backed\_phasing
- [ ] realign
- [x] remove\_end\_tags
- [ ] rename
- [x] replace\_vcf\_sample\_name
- [x] samtools\_flagstat
- [x] samtools\_sort
- [x] select\_variants
- [x] sequence\_align\_and\_tag
- [ ] sequence\_to\_bam
- [x] set\_filter\_status
- [x] single\_sample\_docm\_filter
- [x] smoove
- [ ] somatic\_concordance\_graph
- [ ] sompy
- [x] sort\_vcf
- [x] split\_interval\_list
- [x] split\_interval\_list\_to\_bed
- [x] staged\_rename
- [x] star\_align\_fusion
- [x] star\_fusion\_detect
- [x] strandedness\_check
- [x] strelka
- [x] stringtie
- [x] survivor
- [x] transcript\_to\_gene
- [x] trim\_fastq
- [ ] umi\_align
- [x] variants\_to\_table
- [x] varscan\_germline
- [x] varscan\_process\_somatic
- [x] varscan\_somatic
- [x] vcf\_expression\_annotator
- [x] vcf\_readcount\_annotator
- [x] vcf\_sanitize
- [x] vep
- [x] verify\_bam\_id
- [x] vt\_decompose
- [x] xenosplit
