#!/bin/bash

source activate /ebio/abt6_projects/Ath_HiFi_met1/conda/envs/gffutils

set -xv

cd /ebio/abt6_projects/met1_somatic_transpositions/data/3_somatic-events/r_mutations/insertion-bias/genic-feature_bias/V21_hypervariable_BND_insertion-bias


somins="./merged-runs_V21_hypermutable-region/manual_classification_combined.txt"


genes=liftoff-Genes.gff
genesintr=liftoff-Genes.added-introns.AGAT.gff
repeats="../../Repeats.gff"
TEs=TEs.gff
essentials="../../essential-genes_bias/Lloyd_et_al.2015_known-predicted_lethal-genes.list"

cat ../../TEs.gff | sort -k1,1 -k4,4n > $TEs 

cat /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/6_geneAnno/gene-anno/Tsu-0/2407_liftoff/all-iso/Tsu-0.TAIR10.all-iso.liftoff.gff3 |
 sed 's/,AT.G.*Protein//' | sortBed -i stdin > $genes
[ -s $genesintr ] && [ -f $genesintr ] || agat_sp_add_introns.pl --gff $genes --out $genesintr


awk '{if ($3!="exon" && $3!="gene" && $3 != "mRNA" && $3 != "RNA" && $3 != "protein" ){print $0}}' $genesintr | sortBed -i - > ${genesintr%.gff}.no-exons.gff

# within genes
mkfifo a 
awk '{OFS="\t"}{if ($3!="NA"){print $2,$3,$3,$4,$5,$6,$7}}' $somins | tail -n +2 | grep -v MappingOut | sort -k1,1 -k2,2n > a &
awk '{if ($3!="exon" && $3!="gene"  && $3 != "RNA" && $3 != "protein" ){print $0}}' $genes |
 bedtools intersect -a <(cat a) -b stdin  -wo > tmp
grep -f $essentials tmp  | awk '{OFS="\t"}{print "V21-hypervariable_BND", "essential","gene", $10}'  > genes_nointr.genic-feat.overlap.txt 
grep -v -f $essentials tmp  | awk '{OFS="\t"}{print "V21-hypervariable_BND", "nonessential","gene", $10}'  >> genes_nointr.genic-feat.overlap.txt 
rm tmp & rm a

# exit 1

# within genes with introns
awk '{OFS="\t"}{if ($3!="NA"){print $2,$3,$3,$4,$5,$6,$7}}' $somins | tail -n +2 | grep -v MappingOut | sort -k1,1 -k2,2n | bedtools intersect -a stdin -b ${genesintr%.gff}.no-exons.gff -wo > tmp
grep -f $essentials tmp  | awk '{OFS="\t"}{print "V21-hypervariable_BND", "essential","gene", $10}'  > genes.genic-feat.overlap.txt 
grep -v -f $essentials tmp  | awk '{OFS="\t"}{print "V21-hypervariable_BND", "nonessential","gene", $10}'  >> genes.genic-feat.overlap.txt 
rm tmp

awk '{if ($3=="gene"){print $0}}' $genes > ${gene%.gff}.only-gene.gff

# next to genes
awk '{OFS="\t"}{if ($3!="NA"){print $2,$3,$3,$4,$5,$6,$7}}' $somins | tail -n +2 | grep -v MappingOut | sort -k1,1 -k2,2n | bedtools closest -a stdin -b ${gene%.gff}.only-gene.gff -io -iu -D b > tmp
grep -f $essentials tmp  | awk '{if ($17<0){$17 *= -1 ; print $0}else{print $0}}' | awk '{OFS="\t"}{if ($17<100){print "V21-hypervariable_BND", "essential", $10, "-100"} else if ($17>100 && $17<500){print "V21-hypervariable_BND", "essential", $10, "-500"} else if ($17>500 && $17<1000){print "V21-hypervariable_BND", "essential", $10, "-1k"} else {print "V21-hypervariable_BND", "essential", $10, ">-1k"}}'  > genes.genic-feat.downstream.txt
grep -v -f $essentials tmp  | awk '{if ($17<0){$17 *= -1 ; print $0}else{print $0}}' | awk '{OFS="\t"}{if ($17<100){print "V21-hypervariable_BND", "nonessential", $10, "-100"} else if ($17>100 && $17<500){print "V21-hypervariable_BND", "nonessential", $10, "-500"} else if ($17>500 && $17<1000){print "V21-hypervariable_BND", "nonessential", $10, "-1k"} else {print "V21-hypervariable_BND", "nonessential", $10, ">-1k"}}'  >> genes.genic-feat.downstream.txt
rm tmp

awk '{OFS="\t"}{if ($3!="NA"){print $2,$3,$3,$4,$5,$6,$7}}' $somins | tail -n +2 | grep -v MappingOut | sort -k1,1 -k2,2n | bedtools closest -a stdin -b ${gene%.gff}.only-gene.gff -io -id -D b > tmp
grep -f $essentials tmp  | awk '{if ($17<0){$17 *= -1 ; print $0}else{print $0}}' | awk '{OFS="\t"}{if ($17<100){print "V21-hypervariable_BND", "essential", $10, "+100"} else if ($17>100 && $17<500){print "V21-hypervariable_BND", "essential", $10, "+500"} else if ($17>500 && $17<1000){print "V21-hypervariable_BND", "essential", $10, "+1k"} else {print "V21-hypervariable_BND", "essential", $10, ">+1k"}}'  > genes.genic-feat.upstream.txt
grep -v -f $essentials tmp  | awk '{if ($17<0){$17 *= -1 ; print $0}else{print $0}}' | awk '{OFS="\t"}{if ($17<100){print "V21-hypervariable_BND", "nonessential", $10, "+100"} else if ($17>100 && $17<500){print "V21-hypervariable_BND", "nonessential", $10, "+500"}  else if ($17>500 && $17<1000){print "V21-hypervariable_BND", "nonessential", $10, "+1k"} else {print "V21-hypervariable_BND", "nonessential", $10, ">+1k"}}'  >> genes.genic-feat.upstream.txt
rm tmp

# exit 1

# within TEs
awk '{OFS="\t"}{if ($3!="NA"){print $2,$3,$3,$4,$5,$6,$7}}' $somins | tail -n +2 | grep -v MappingOut | sort -k1,1 -k2,2n | bedtools intersect -a stdin -b ${TEs} -wo | sort |  awk '{OFS="\t"}{print "V21-hypervariable_BND",$10,"TE","within"}'  > TEs.genic-feat.overlap.txt

# next to TEs
awk '{OFS="\t"}{if ($3!="NA"){print $2,$3,$3,$4,$5,$6,$7}}' $somins | tail -n +2 | grep -v MappingOut | sort -k1,1 -k2,2n | bedtools closest -a stdin -b ${TEs} -io -iu -D b | sort > tmp
cat tmp |awk '{if ($17<0){$17 *= -1 ; print $0}else{print $0}}' | awk '{OFS="\t"}{if ($17<100){print "V21-hypervariable_BND", $10,"TE", "-100"} 
 else if ($17>100 && $17<500){print "V21-hypervariable_BND",  $10, "TE", "-500"} 
 else if ($17>500 && $17<1000){print "V21-hypervariable_BND", $10, "TE", "-1k"} 
 else {print "V21-hypervariable_BND", $10, "TE", ">-1k"}}'  > TEs.genic-feat.downstream.txt
rm tmp

awk '{OFS="\t"}{if ($3!="NA"){print $2,$3,$3,$4,$5,$6,$7}}' $somins | tail -n +2 | grep -v MappingOut | sort -k1,1 -k2,2n | bedtools closest -a stdin -b ${TEs} -io -id -D b | sort > tmp
cat tmp |awk '{if ($17<0){$17 *= -1 ; print $0}else{print $0}}' | awk '{OFS="\t"}{if ($17<100){print "V21-hypervariable_BND", $10,"TE", "+100"} 
 else if ($17>100 && $17<500){print "V21-hypervariable_BND",  $10, "TE", "+500"} 
 else if ($17>500 && $17<1000){print "V21-hypervariable_BND", $10, "TE", "+1k"} 
 else {print "V21-hypervariable_BND", $10, "TE", ">+1k"}}'  > TEs.genic-feat.upstream.txt
rm tmp



# within REPEATS
awk '{OFS="\t"}{if ($3!="NA"){print $2,$3,$3,$4,$5,$6,$7}}' $somins | tail -n +2 | grep -v MappingOut | sort -k1,1 -k2,2n | bedtools intersect -a stdin -b $repeats -wo | sort -k1,1 |
  awk '{OFS="\t"}{print "V21-hypervariable_BND",$10,"Repeat","within"}'  > repeats.genic-feat.overlap.txt

# next to REPEATS
awk '{OFS="\t"}{if ($3!="NA"){print $2,$3,$3,$4,$5,$6,$7}}' $somins | tail -n +2 | grep -v MappingOut | sort -k1,1 -k2,2n | bedtools closest -a stdin -b $repeats -io  -D b | sort -k1,1 > tmp
cat tmp |awk '{if ($17<0){$17 *= -1 ; print $0}else{print $0}}' | awk '{OFS="\t"}{if ($17<100){print "V21-hypervariable_BND", $10,"Repeat", "100"} 
 else if ($17>100 && $17<500){print "V21-hypervariable_BND",  $10, "Repeat", "500"} 
 else if ($17>500 && $17<1000){print "V21-hypervariable_BND", $10, "Repeat", "1k"} 
 else {print "V21-hypervariable_BND", $10, "Repeat", ">1k"}}'  > repeats.genic-feat.down-upstream.txt
rm tmp


# completely intergenic
awk '{OFS="\t"}{if ($3!="NA"){print $2,$3,$3,$4,$5,$6,$7}}' $somins | tail -n +2 | grep -v MappingOut | sort -k1,1 -k2,2n | bedtools intersect -v -a stdin -b $TEs ${genesintr%.gff}.no-exons.gff $repeats -wo | sort -k1,1 |
 awk '{OFS="\t"}{print "V21-hypervariable_BND","intergenic","intergenic","within"}'  > intergenic.genic-feat.overlap.txt

set +xv
conda deactivate