#!/bin/bash

source activate /ebio/abt6_projects/Ath_HiFi_met1/conda/envs/gffutils

set -xv

cd /ebio/abt6_projects/met1_somatic_transpositions/data/3_somatic-events/r_mutations/insertion-bias/genic-feature_bias


somins=SomaticInsertions_CURATED-list.no-comm.mod.24-06-14.tsv


genes=liftoff-Genes.gff
genesintr=liftoff-Genes.added-introns.AGAT.gff
repeats=../Repeats.gff
TEs=TEs.gff
essentials=../essential-genes_bias/Lloyd_et_al.2015_known-predicted_lethal-genes.list

cat ../TEs.gff | sort -k1,1 -k4,4n > $TEs 

cat /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/6_geneAnno/gene-anno/Tsu-0/2407_liftoff/all-iso/Tsu-0.TAIR10.all-iso.liftoff.gff3 |
 sed 's/,AT.G.*Protein//' > $genes
[ -s $genesintr ] && [ -f $genesintr ] || agat_sp_add_introns.pl --gff $genes --out $genesintr


awk '{if ($3!="exon" && $3!="gene" && $3 != "mRNA" && $3 != "RNA" && $3 != "protein" ){print $0}}' $genesintr > ${genesintr%.gff}.no-exons.gff

# within genes
mkfifo a 
grep CONFIRMED $somins | cut -f1-13 | sort -k1,1 -k2,2n > a &
awk '{if ($3!="exon" && $3!="gene"  && $3 != "RNA" && $3 != "protein" ){print $0}}' $genes |
 bedtools intersect -a <(cat a) -b stdin  -wo > tmp
grep -f $essentials tmp  | awk '{OFS="\t"}{print $11, "essential","gene", $16}'  > genes_nointr.genic-feat.overlap.txt 
grep -v -f $essentials tmp  | awk '{OFS="\t"}{print $11, "nonessential","gene", $16}'  >> genes_nointr.genic-feat.overlap.txt 
rm tmp & rm a

# within genes with introns
grep CONFIRMED $somins | cut -f1-13 | sort -k1,1 -k2,2n | bedtools intersect -a stdin -b ${genesintr%.gff}.no-exons.gff -wo > tmp
grep -f $essentials tmp  | awk '{OFS="\t"}{print $11, "essential","gene", $16}'  > genes.genic-feat.overlap.txt 
grep -v -f $essentials tmp  | awk '{OFS="\t"}{print $11, "nonessential","gene", $16}'  >> genes.genic-feat.overlap.txt 
rm tmp

# next to genes
grep CONFIRMED $somins | cut -f1-13 | sort -k1,1 -k2,2n | bedtools closest -a stdin -b ${genesintr%.gff}.only-gene.gff -io -iu -D b > tmp
grep -f $essentials tmp  | awk '{if ($23<0){$23 *= -1 ; print $0}else{print $0}}' | awk '{OFS="\t"}{if ($23<100){print $11, "essential", $16, "-100"} else if ($23>100 && $23<500){print $11, "essential", $16, "-500"} else if ($23>500 && $23<1000){print $11, "essential", $16, "-1k"} else {print $11, "essential", $16, ">-1k"}}'  > genes.genic-feat.downstream.txt
grep -v -f $essentials tmp  | awk '{if ($23<0){$23 *= -1 ; print $0}else{print $0}}' | awk '{OFS="\t"}{if ($23<100){print $11, "nonessential", $16, "-100"} else if ($23>100 && $23<500){print $11, "nonessential", $16, "-500"} else if ($23>500 && $23<1000){print $11, "nonessential", $16, "-1k"} else {print $11, "nonessential", $16, ">-1k"}}'  >> genes.genic-feat.downstream.txt
rm tmp

grep CONFIRMED $somins | cut -f1-13 | sort -k1,1 -k2,2n | bedtools closest -a stdin -b ${genesintr%.gff}.only-gene.gff -io -id -D b > tmp
grep -f $essentials tmp  | awk '{if ($23<0){$23 *= -1 ; print $0}else{print $0}}' | awk '{OFS="\t"}{if ($23<100){print $11, "essential", $16, "+100"} else if ($23>100 && $23<500){print $11, "essential", $16, "+500"} else if ($23>500 && $23<1000){print $11, "essential", $16, "+1k"} else {print $11, "essential", $16, ">+1k"}}'  > genes.genic-feat.upstream.txt
grep -v -f $essentials tmp  | awk '{if ($23<0){$23 *= -1 ; print $0}else{print $0}}' | awk '{OFS="\t"}{if ($23<100){print $11, "nonessential", $16, "+100"} else if ($23>100 && $23<500){print $11, "nonessential", $16, "+500"}  else if ($23>500 && $23<1000){print $11, "nonessential", $16, "+1k"} else {print $11, "nonessential", $16, ">+1k"}}'  >> genes.genic-feat.upstream.txt
rm tmp


# within TEs
grep CONFIRMED $somins | cut -f1-13 | sort -k1,1 -k2,2n | bedtools intersect -a stdin -b ${TEs} -wo | sort -k12,12 -u |  awk '{OFS="\t"}{print $11,$16,"TE","within"}'  > TEs.genic-feat.overlap.txt

# next to TEs
grep CONFIRMED $somins | cut -f1-13 | sort -k1,1 -k2,2n | bedtools closest -a stdin -b ${TEs} -io -iu -D b | sort -k12,12 -u > tmp
cat tmp |awk '{if ($23<0){$23 *= -1 ; print $0}else{print $0}}' | awk '{OFS="\t"}{if ($23<100){print $11, $16,"TE", "-100"} 
 else if ($23>100 && $23<500){print $11,  $16, "TE", "-500"} 
 else if ($23>500 && $23<1000){print $11, $16, "TE", "-1k"} 
 else {print $11, $16, "TE", ">-1k"}}'  > TEs.genic-feat.downstream.txt
rm tmp

grep CONFIRMED $somins | cut -f1-13 | sort -k1,1 -k2,2n | bedtools closest -a stdin -b ${TEs} -io -id -D b | sort -k12,12 -u > tmp
cat tmp |awk '{if ($23<0){$23 *= -1 ; print $0}else{print $0}}' | awk '{OFS="\t"}{if ($23<100){print $11, $16,"TE", "+100"} 
 else if ($23>100 && $23<500){print $11,  $16, "TE", "+500"} 
 else if ($23>500 && $23<1000){print $11, $16, "TE", "+1k"} 
 else {print $11, $16, "TE", ">+1k"}}'  > TEs.genic-feat.upstream.txt
rm tmp



# within REPEATS
grep CONFIRMED $somins | cut -f1-13 | sort -k1,1 -k2,2n | bedtools intersect -a stdin -b $repeats -wo | sort -k1,1 -k12,12 -u |
  awk '{OFS="\t"}{print $11,$16,"Repeat","within"}'  > repeats.genic-feat.overlap.txt

# next to REPEATS
grep CONFIRMED $somins | cut -f1-13 | sort -k1,1 -k2,2n | bedtools closest -a stdin -b $repeats -io  -D b | sort -k1,1 -k12,12 -u > tmp
cat tmp |awk '{if ($23<0){$23 *= -1 ; print $0}else{print $0}}' | awk '{OFS="\t"}{if ($23<100){print $11, $16,"Repeat", "100"} 
 else if ($23>100 && $23<500){print $11,  $16, "Repeat", "500"} 
 else if ($23>500 && $23<1000){print $11, $16, "Repeat", "1k"} 
 else {print $11, $16, "Repeat", ">1k"}}'  > repeats.genic-feat.down-upstream.txt
rm tmp


# completely intergenic
grep CONFIRMED $somins | cut -f1-13 | sort -k1,1 -k2,2n | bedtools intersect -v -a stdin -b $TEs ${genesintr%.gff}.no-exons.gff $repeats -wo | sort -k1,1 -k12,12 -u |
 awk '{OFS="\t"}{print $11,"intergenic","intergenic","within"}'  > intergenic.genic-feat.overlap.txt

set +xv
conda deactivate