#!/bin/bash

SRC_DIR=$(dirname "$0")

show_help () {
    cat <<EOF
usage: sh $0 --ARGUMENT <value>

arguments:
  -h, --help      prints this block and immediately exists

  --phlat-dir     the path to where phlat is located, DEFAULT \$SRC_DIR/phlat-release
  --data-dir      the path to the fastq data, DEFAULT example
  --samtools      the path to the samtools executable, DEFAULT /usr/local/bin/samtools
  --bam           the path to where the bam file is located, NO DEFAULT 
  --index-dir     the path to b2folder, DEFAULT \$PHLAT_DIR/b2folder
  --rs-dir        the path where the results will be, DEFAULT \$DATA_DIR/results
  --b2url         the path to where bowtie2 is located, DEFAULT /usr/bin/bowtie2 
  --fastq1        the name of seq num 1, DEFAULT 'example_1.fastq.gz' 
  --fastq2        the name of seq num 2, DEFAULT 'example_2.fastq.gz'
  --ref-fasta     the path to where the reference fasta file is located, NO DEFAULT


EOF
}

# die and opts based on this example
# http://mywiki.wooledge.org/BashFAQ/035
# --long-opt* example here
# https://stackoverflow.com/a/7069755
function die {
    printf '%s\n' "$1" >&2 && show_help && exit 1
}

# check arguments
while test $# -gt 0; do
    case $1 in
        -h|--help)
            show_help
            exit
            ;;
         --phlat-dir*)
            if [ ! "$2" ]; then
                PHLAT_DIR=""
            else
                PHLAT_DIR=$2
                shift
            fi
            ;;
         --data-dir*)
	    if [ ! "$2" ]; then
		DARA_DIR=""
	    else
		DATA_DIR=$2
		shift
	    fi
	    ;;
         --samtools*)
	    if [ ! "$2" ]; then
		SAMTOOLS=""
	    else
		SAMTOOLS=$2
		shift
	    fi
	    ;;
         --bam*)
            if [ ! "$2" ]; then 
                die 'ERROR: "--bam" requires a non-empty argument.'
            else
                BAM=$2
                shift
            fi
            ;;
        --index-dir*)
	    if [ ! "$2" ]; then
		INDEX_DIR=""
	    else
		INDEX_DIR=$2
		shift
	    fi
       	    ;;
	--rs-dir*)
	    if [ ! "$2" ]; then
		RS_DIR=""
	    else
		RS_DIR=$2
		shift
	    fi
       	    ;;
	--b2url*)
	    if [ ! "$2" ]; then
		B2URL=""
	    else
		B2URL=$2
		shift
	    fi
       	    ;;
        --fastq1*)
	    if [ ! "$2" ]; then
		FASTQ1=""
	    else
		FASTQ1=$2
		shift
	    fi
       	    ;;
       --fastq2*)
	    if [ ! "$2" ]; then
		FASTQ2=""
	    else
		FASTQ2=$2
		shift
	    fi
       	    ;;
        --ref-fasta*)
	    if [ ! "$2" ]; then
		REF_FASTA="ERROR: --ref-fasta requires non-empty argument"
	    else
		REF_FASTA=$2
		shift
	    fi
       	    ;; 

        *)
            break
            ;;
    esac
    shift
done

# double check all vars are set up
[ -z $PHLAT_DIR    ] && PHLAT_DIR="$SRC_DIR/phlat-release"
[ -z $DATA_DIR     ] && DATA_DIR="example"
[ -z $SAMTOOLS     ] && SAMTOOLS="/usr/local/bin/samtools"
[ -z $BAM          ] && die "Missing argument --bam"
[ -z $INDEX_DIR    ] && INDEX_DIR="$PHLAT_DIR/b2folder"
[ -z $RS_DIR       ] && RS_DIR="$DATA_DIR/results" 
[ -z $B2URL        ] && B2URL="/usr/bin/bowtie2"
[ -z $FASTQ1       ] && FASTQ1="example_1.fastq.gz"
[ -z $FASTQ2       ] && FASTQ2="example_2.fastq.gz"
[ -z $REF_FASTA    ] && die "Missing argument --ref-fasta"


mkdir -p $DATA_DIR/results
tmpdir=$DATA_DIR/tmp
mkdir -p $tmpdir

# extract hla regions and unmapped reads
echo "extracting hla region and unmapped reads ..."
$SAMTOOLS view -h -T $REF_FASTA $BAM chr6:29836259-33148325 >$tmpdir/reads.sam

$SAMTOOLS view -H -T $REF_FASTA $BAM | grep "^@SQ" | cut -f 2 | cut -f 2- -d : | grep HLA | while read chr;do 
# echo "checking $chr:1-9999999"
$SAMTOOLS view -T $REF_FASTA $BAM "$chr:1-9999999" >>$tmpdir/reads.sam
done

$SAMTOOLS view -f 4 -T $REF_FASTA $BAM >>$tmpdir/reads.sam
$SAMTOOLS view -Sb -o $tmpdir/reads.bam $tmpdir/reads.sam 

echo "running pircard..."
/usr/bin/java -Xmx6g -jar /usr/picard/picard.jar SamToFastq VALIDATION_STRINGENCY=LENIENT F=$DATA_DIR/hlaPlusUnmapped_1.fastq.gz F2=$DATA_DIR/hlaPlusUnmapped_2.fastq.gz I=$tmpdir/reads.bam R=$REF_FASTA FU=$DATA_DIR/unpaired.fastq.gz

#workaround to get everything passed in appropriately
echo "running PHLAT ..."
python2 -O ${PHLAT_DIR}/dist/PHLAT.py -1 ${DATA_DIR}/${FASTQ1} -2 ${DATA_DIR}/${FASTQ2} -index $INDEX_DIR -b2url $B2URL -orientation "--fr" -tag $DATA_DIR -e $PHLAT_DIR -o $RS_DIR -tmp 0 -p 4 >$DATA_DIR/run_phlat.sh
