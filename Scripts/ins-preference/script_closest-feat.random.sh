#!/bin/bash

source activate /ebio/abt6_projects/Ath_HiFi_met1/conda/envs/gffutils

set -xv

cd /ebio/abt6_projects/met1_somatic_transpositions/data/3_somatic-events/r_mutations/insertion-bias/genic-feature_bias


bedtools random -l 0 -n 1000 -g /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/2_hifiasm_ragtag/Tsu-0/Tsu-0.Chr_scaffolds.v2.genome > random-insertions.1k.bed

ran=random-insertions.1k.bed


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
cat $ran | sort -k1,1 -k2,2n > a &
awk '{if ($3!="exon" && $3!="gene" && $3 != "mRNA" && $3 != "RNA" && $3 != "protein" ){print $0}}' $genes |
 bedtools intersect -a <(cat a) -b stdin  -wo > tmpr
grep -f $essentials tmpr  | awk '{OFS="\t"}{print "random", "essential","gene", $9}'  > genes_nointr.genic-feat.overlap.RANDOM.txt 
grep -v -f $essentials tmpr  | awk '{OFS="\t"}{print "random", "nonessential","gene", $9}'  >> genes_nointr.genic-feat.overlap.RANDOM.txt 
rm tmpr & rm a

# within genes with introns
cat $ran | sort -k1,1 -k2,2n | bedtools intersect -a stdin -b ${genesintr%.gff}.no-exons.gff -wo > tmpr
grep -f $essentials tmpr  | awk '{OFS="\t"}{print "random", "essential","gene", $9}'  > genes.genic-feat.overlap.RANDOM.txt 
grep -v -f $essentials tmpr  | awk '{OFS="\t"}{print "random", "nonessential","gene", $9}'  >> genes.genic-feat.overlap.RANDOM.txt 
rm tmpr

# next to genes
cat $ran | sort -k1,1 -k2,2n | bedtools closest -a stdin -b ${genesintr%.gff}.only-gene.gff -io -iu -D b > tmpr
grep -f $essentials tmpr  | awk '{if ($16<0){$16 *= -1 ; print $0}else{print $0}}' | awk '{OFS="\t"}{if ($16<100){print "random", "essential", $9, "-100"} else if ($16>100 && $16<500){print "random", "essential", $9, "-500"} else if ($16>500 && $16<1000){print "random", "essential", $9, "-1k"} else {print "random", "essential", $9, ">-1k"}}'  > genes.genic-feat.downstream.RANDOM.txt
grep -v -f $essentials tmpr  | awk '{if ($16<0){$16 *= -1 ; print $0}else{print $0}}' | awk '{OFS="\t"}{if ($16<100){print "random", "nonessential", $9, "-100"} else if ($16>100 && $16<500){print "random", "nonessential", $9, "-500"} else if ($16>500 && $16<1000){print "random", "nonessential", $9, "-1k"} else {print "random", "nonessential", $9, ">-1k"}}'  >> genes.genic-feat.downstream.RANDOM.txt
rm tmpr

cat $ran | sort -k1,1 -k2,2n | bedtools closest -a stdin -b ${genesintr%.gff}.only-gene.gff -io -id -D b > tmpr
grep -f $essentials tmpr  | awk '{if ($16<0){$16 *= -1 ; print $0}else{print $0}}' | awk '{OFS="\t"}{if ($16<100){print "random", "essential", $9, "+100"} else if ($16>100 && $16<500){print "random", "essential", $9, "+500"} else if ($16>500 && $16<1000){print "random", "essential", $9, "+1k"} else {print "random", "essential", $9, ">+1k"}}'  > genes.genic-feat.upstream.RANDOM.txt
grep -v -f $essentials tmpr  | awk '{if ($16<0){$16 *= -1 ; print $0}else{print $0}}' | awk '{OFS="\t"}{if ($16<100){print "random", "nonessential", $9, "+100"} else if ($16>100 && $16<500){print "random", "nonessential", $9, "+500"}  else if ($16>500 && $16<1000){print "random", "nonessential", $9, "+1k"} else {print "random", "nonessential", $9, ">+1k"}}'  >> genes.genic-feat.upstream.RANDOM.txt
rm tmpr


# within TEs
cat $ran | sort -k1,1 -k2,2n | bedtools intersect -a stdin -b ${TEs} -wo  |  awk '{OFS="\t"}{print "random",$9,"TE","within"}'  > TEs.genic-feat.overlap.RANDOM.txt

# next to TEs
cat $ran | sort -k1,1 -k2,2n | bedtools closest -a stdin -b ${TEs} -io -iu -D b  > tmpr
cat tmpr |awk '{if ($16<0){$16 *= -1 ; print $0}else{print $0}}' | awk '{OFS="\t"}{if ($16<100){print "random", $9,"TE", "-100"} 
 else if ($16>100 && $16<500){print "random",  $9, "TE", "-500"} 
 else if ($16>500 && $16<1000){print "random", $9, "TE", "-1k"} 
 else {print "random", $9, "TE", ">-1k"}}'  > TEs.genic-feat.downstream.RANDOM.txt
rm tmpr

cat $ran | sort -k1,1 -k2,2n | bedtools closest -a stdin -b ${TEs} -io -id -D b  > tmpr
cat tmpr |awk '{if ($16<0){$16 *= -1 ; print $0}else{print $0}}' | awk '{OFS="\t"}{if ($16<100){print "random", $9,"TE", "+100"} 
 else if ($16>100 && $16<500){print "random",  $9, "TE", "+500"} 
 else if ($16>500 && $16<1000){print "random", $9, "TE", "+1k"} 
 else {print "random", $9, "TE", ">+1k"}}'  > TEs.genic-feat.upstream.RANDOM.txt
rm tmpr



# within REPEATS
cat $ran | sort -k1,1 -k2,2n | bedtools intersect -a stdin -b $repeats -wo  |
  awk '{OFS="\t"}{print "random",$9,"Repeat","within"}'  > repeats.genic-feat.overlap.RANDOM.txt

# next to REPEATS
cat $ran | sort -k1,1 -k2,2n | bedtools closest -a stdin -b $repeats -io  -D b  > tmpr
cat tmpr |awk '{if ($16<0){$16 *= -1 ; print $0}else{print $0}}' | awk '{OFS="\t"}{if ($16<100){print "random", $9,"Repeat", "100"} 
 else if ($16>100 && $16<500){print "random",  $9, "Repeat", "500"} 
 else if ($16>500 && $16<1000){print "random", $9, "Repeat", "1k"} 
 else {print "random", $9, "Repeat", ">1k"}}'  > repeats.genic-feat.down-upstream.RANDOM.txt
rm tmpr


# completely intergenic
cat $ran | sort -k1,1 -k2,2n | bedtools intersect -v -a stdin -b $TEs ${genesintr%.gff}.no-exons.gff $repeats -wo  |
 awk '{OFS="\t"}{print "random","intergenic","intergenic","within"}'  > intergenic.genic-feat.overlap.RANDOM.txt

set +xv
conda deactivate