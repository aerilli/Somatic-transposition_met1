#!/bin/bash

# amovilli

## for detecting excised DNA TEs


set -xvu
# set -o pipefail

### EXAMPLE ###
## for SP in `cat ./list.SINGLE-PLANTS.txt` ; do
### ./Scripts/script_calling_somatic-insertions.CIGAR.sh \
###  -s  ${SP} \
###  -b  ./Data/ALN/${SP}-aln2-WT-Chr_winnowmap.sortedRG.bam \
###  -o  ./Data/Somatic_events/insertions_CIGAR \
###  -t  ./Data/Tsu-0.TE-anno.manual.gff \
###  -f  ./Data/Tsu-0.repeats.merged.gff \
###  -R  ./Reads/Tsu-0.reference.fa \
###   2>&1 | tee -a ./Data/logs/$(date '+%Y-%m-%d_%H-%M')_calling_TE_insertions_${SP}.log
## done



# setting variables


while getopts "s:b:o:t:f:R:h" opt; do
 case $opt in
s)
  singleplant=$OPTARG >&2
  ;;
b)
  BAMIN=$OPTARG >&2
  ;;
o)
  outdir=$OPTARG >&2
  ;;
t)
  TEanno=$OPTARG >&2
  ;;
f)
  FREPEATS=$OPTARG >&2
  ;;
R)
  REF=$OPTARG >&2
  ;;

h)
  echo "Help"
	exit 1
 esac
done


# Creating dirs and setting other variables 
mkdir -p $outdir
cd $outdir

BAMOUT=`basename ${BAMIN%.*}`
# FREPEATSNS=${FREPEATS%.gff*}.Nslop50.nocentr.gff

# Virtual env activation:
# containing bioawk; bedtools v2.30.0; samtools v1.16.1
## source activate /ebio/abt6_projects/met1_somatic_transpositions/conda/so_somatic-events



############## INS ##############


# Getting reads with CIGAR strings with >500bp deletions
samtools view -@ 32 -q 50 -h $BAMIN | awk 'BEGIN{OFS="\t"} { if ($6 ~ /[5-9][0-9][0-9]+I |[1-9][0-9][0-9][0-9]+I/) {print $0} else if ($1 ~ /^@/){print $0} }' |
 samtools view -@ 32 -bS > $BAMOUT.insertion.500.bam

samtools index $BAMOUT.insertion.500.bam

# reformatting using '|' ; substituting with tab ; removing parts of string that are not needed
## then retrieving coordinates of insertion within the read: A & B
samtools view -@32  $BAMOUT.insertion.500.bam | sed -E 's/[0-9]{3,5}I/\|\|\|/g' |  
 awk '{if ($6!~/.*\|\|\|.*\|\|\|.*/){print $0}}'  | tee ${BAMOUT}.insertion.500.mod.formatted.sam |
 sed 's/|||/\t/g' | cut -f6 | sed 's/[MIS=X]/\+/g;s/[HDPN]/\*0+/g;s/$/0/g' | bc > distance-I.A.tmp
 
cat ${BAMOUT}.insertion.500.mod.formatted.sam | sed 's/|||/\t/' | cut -f7 |
 sed 's/[MIS=X]/\+/g;s/[HDPN]/\*0+/g;s/$/0/g'  | bc > distance-I.B.tmp

cat ${BAMOUT}.insertion.500.mod.formatted.sam | sed 's/|||/\t/' | cut -f6 |
 sed 's/[M=XD]/\+/g;s/[HSPIN]/\*0+/g;s/$/0/g'  | bc > distance4coord-I.A.tmp

cat ${BAMOUT}.insertion.500.mod.formatted.sam | sed 's/|||/\t/' | cut -f7 |
 sed 's/[M=XD]/\+/g;s/[HSPIN]/\*0+/g;s/$/0/g'  | bc > distance4coord-I.B.tmp



# With coordinates, getting bed and fasta of insertions

##bed 
paste distance-I.A.tmp distance-I.B.tmp distance4coord-I.A.tmp distance4coord-I.B.tmp ${BAMOUT}.insertion.500.mod.formatted.sam  |
 awk '{OFS="\t"}{if ($6==0 || $6==2048){print $5 "|" $7 ":" $8 + $3 "-" $8 + $3 "__", $1,length($14) - $2 } else if ($6==16 || $6==2064){print $5 "|" $7 ":" $8 + $3 "-" $8 + $3 "__", $1,length($14)-$2}}' \
  > ${BAMOUT}.CIGAR.insertions.bed

##fasta
paste distance-I.A.tmp distance-I.B.tmp distance4coord-I.A.tmp distance4coord-I.B.tmp ${BAMOUT}.insertion.500.mod.formatted.sam |
 awk '{OFS="\t"}{if ($6==0 || $6==2048){print ">" $5 "|" $7 ":" $8 + $3 "-" $8 + $3 "__" "\n" $14} else if ($6==16 || $6==2064){print ">" $5 "|" $7 ":" $8 + $3 "-" $8 + $3 "__" "\n" $14}}' \
  >  ${BAMOUT}.CIGAR.reads.insertions.fa
  
samtools faidx ${BAMOUT}.CIGAR.reads.insertions.fa
bedtools getfasta -fi ${BAMOUT}.CIGAR.reads.insertions.fa  -bed ${BAMOUT}.CIGAR.insertions.bed |
 sed 's/__.*//' > ${BAMOUT}.CIGAR.insertions.fa


# Mapping insertion fasta to the reference genome
minimap2 -t 32 -ax map-pb $REF \
 ${BAMOUT}.CIGAR.insertions.fa |
 samtools view -h -b -@32 | #-F 256 |
 samtools sort -O SAM - |
 paftools.js sam2paf -L - \
 >  $BAMOUT.reads.withINS.paf


# awk rearrange for preparing bed file for intersecting with EDTA intact aka TE annotation
awk 'BEGIN{OFS="\t"; print "#" "tname","tstart","tend","qname","qstart","qend","qlen","strand","tlen","matches","lalnblock","Q","tp","mm","gn","go", "cg"}
    {print $6,$8,$9,  $1, $3,$4, $2, $5, $7,$10,$11,$12,$13,$14,$15,$16,$17}' $BAMOUT.reads.withINS.paf | sortBed \
    > $BAMOUT.reads.withINS.TE-format.bed # header is lost anyways # bed file for input bedtools intersect with EDTA intacts.



#intersecting with EDTA intact
grep -v -e ChrM -e ptg  $TEanno | # Have to remove ChrM and ptg entries first
bedtools intersect -v -a - -b $FREPEATS | # removing TEs that overlap with FRepeats (rDNA, slop Nstretches, organellar, telomeres) from TE intact; centromeres included though
bedtools intersect -a $BAMOUT.reads.withINS.TE-format.bed  -b stdin -wa -wb -loj |
awk '{if($12>50){print $0}}'> $BAMOUT.reads.withINS.TE-format.re-overlap.bed   # loj==left outer join


sed 's/|/\t/;s/:/\t/;s/-/\t/' $BAMOUT.reads.withINS.TE-format.re-overlap.bed |
 awk 'BEGIN{ print "#ChrIns\tStartIns\tEndIns\tReadID\tChrOrigin\tStartOrigin\tEndOrigin\tTEannoOrigin"} ;  {print $5 "\t" $6 "\t" $7 "\t"  $4 "\t" $1 "\t" $2 "\t" $3 "\t" $21"|"$24"|"$25-1"|"$29}' |
 sortBed -header |
 bedtools merge -i stdin -c 4,5,6,7,8 -o distinct,distinct,distinct,distinct,collapse, -delim "__" -header \
  > $BAMOUT.CG.insertions.final.bed




## conda deactivate

