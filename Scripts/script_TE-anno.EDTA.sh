#!/bin/bash

# Tsu-0 EDTA annotation

set -xvu


# EXAMPLE # ./Scripts/script_TE-anno.EDTA.sh \
#           ./Data/Tsu-0.reference.fa \
#           ./Data/TE-anno \
#			5 \
#           16


# Parse parameteres
REF=$1
OUTDIR=$2
DIV=${3:-5}
CORES=${4:-${NSLOTS:-12}}
LIB="./Data/Col-CC.cur-lib.fa"

# Variables to be set
EDTApath= 


# Create necessary directories and change dir
mkdir -p $OUTDIR
cd $OUTDIR

# Virtual env activation:
# containing EDTA v2.2; gffread; bedtools
## source activate /ebio/abt6_projects/Ath_HiFi_met1/conda/envs/EDTA_v2.2

# CDS.fa generation
gffread -x ./Data/TAIR10.gene-anno.CDS.fa -g ./Data/TAIR10.fa ./Data/TAIR10.gene-anno.gff

# inplace change for specification of divergence
echo "Divergence is set at ${DIV}..."
sed -i "s/genome.out.new -maxdiv [0-9]*\.[0-9]*/genome.out.new -maxdiv ${DIV}/" ${EDTApath} #/ebio/abt6_projects/Ath_HiFi_met1/conda/software/EDTA/EDTA.pl


echo "Running EDTA"
perl /ebio/abt6_projects/Ath_HiFi_met1/conda/software/EDTA/EDTA.pl --threads $CORES \
--genome $REF \
--sensitive 1 \
--anno 1 \
--overwrite 0 \
--evaluate 1 \
--debug 1 \
--cds ./Data/TAIR10.gene-anno.CDS.fa \
--curatedlib $LIB




# cleaning TE annotation from repeats overlap
INPUT=$OUTDIR/Tsu-0.reference.fa.mod.EDTA.TEanno.gff3
outnorepeats=$WORKDIR/`basename $INPUT ".gff3"`.no-repeats.gff3
FREPEATS=./Data/Tsu-0.repeats.merged.gff

## TEs have to overlap at least 50% to be removed
bedtools intersect -header -v -f 0.5 -a $INPUT -b $FREPEATS > $outnorepeats




# merging TE annotations by name
## This step has been done to ensure VANDAL merging and might not be good for other uses
INPUT=$outnorepeats
names=`grep "Name="  ${INPUT} |grep -v  "target_site_duplication" |grep -v "long_terminal_repeat"  | awk '{if ($3!="repeat_region"){print $0}}'| grep -v "ptg" | sed -E 's/Parent=repeat_region_[0-9]*;//' |  sed 's/;/\t/g'  | cut -f10 | sort | uniq`
output=${INPUT%.gff*}.merged-overlapping-elements.gff

grep "^#" $INPUT > $output.tmp

for i in $names ; do
  echo $i | sed 's/Name=//' 
  cat $INPUT   |grep -v  "target_site_duplication" |grep -v "long_terminal_repeat" | awk '{if ($3!="repeat_region"){print $0}}'| grep -v "^#" | grep -v "ptg" |
  grep ${i} | 
  bedtools merge -s -i stdin -c 1,2,3,4,5,6,7,8,9 -o distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct -delim "|" |
  awk -F'\t' -v OFS='\t' '{sub(/\|.*/, "",$2)} 1' | awk -F'\t' -v OFS='\t' '{sub(/.*\|/, "",$3)} 1' |#tr "|" "\t" |
  awk '{OFS="\t"}{ if ($3 < $2){print $1, $5, $6, $3, $2, $9, $10 , $11 , $12} else if ($3 > $2){print $1, $5, $6, $2, $3, $9, $10 , $11 , $12}}' >> $output.tmp
done  

awk '{OFS="\t"}{ if ($4==0){print $1 , $2, $3 , 1 , $5 , $6 ,$7 , $8 ,$9} else {print $0}}'  $output.tmp |
 bedtools intersect -v -a stdin -b $FREPEATS |
 sortBed   | uniq | sed -E 's/Parent=repeat_region_[0-9]*;//' > $output && rm $output.tmp



# Virtual env deactivation
## conda deactivate



############## TE annotation has been manually refined for those TEs that have been found to be mobile ##############
# Tsu-0.TE-anno.manual.gff has been generated
# Manually annotated TEs can be retrieved with `grep MANUAL Tsu-0.TE-anno.manual.gff`