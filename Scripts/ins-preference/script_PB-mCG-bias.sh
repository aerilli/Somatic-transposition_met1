#!/bin/bash


cd /ebio/abt6_projects/met1_somatic_transpositions/data/3_somatic-events/r_mutations/insertion-bias/mCG_bias

set -xv

IND=`echo "met1_0"{1..9}`" met1_10"

for i in $IND ; do
    echo ""
    echo $i

    B=`grep "$i" /ebio/abt6_projects/met1_somatic_transpositions/data/3_somatic-events/r_mutations/insertion-bias/mCG_bias/PB_mCG/demux.txt | cut -f1`
    
    # random generation of genomic positions
    bedtools random -g /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/2_hifiasm_ragtag/Tsu-0/Tsu-0.Chr_scaffolds.v2.genome -n 10000 -l 0 |
     bedtools slop -i stdin -g /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/2_hifiasm_ragtag/Tsu-0/Tsu-0.Chr_scaffolds.v2.genome -b 25 |
     bedtools intersect -a stdin -b ./PB_mCG/tsumet.merged.ccs.hifi.demux.${B}--${B}.jasmine.pbmm2.combined.bedGraph -loj | 
     awk -v I="${i}" '{OFS="\t"}{if ($10=="."){print $1,$2,$3,"random","random","0",I} else {print $1,$2,$3,"random","random",$10,I}}' \
        > PBmCG.10k_random-positions.bslop25bp.${i}.bed

    ## excluding repeat regions
    bedtools random -g /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/2_hifiasm_ragtag/Tsu-0/Tsu-0.Chr_scaffolds.v2.genome -n 1000000 -l 0 |
     bedtools slop -i stdin -g /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/2_hifiasm_ragtag/Tsu-0/Tsu-0.Chr_scaffolds.v2.genome -b 25 |
     bedtools intersect -v -a stdin  -b /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/3_annotation_REPEATS/Chr_scaffolds/Tsu-0/Tsu-0.Repeats_merged.gff |
     shuf -n 10000 | sortBed -i stdin |
     bedtools intersect -a stdin -b ./PB_mCG/tsumet.merged.ccs.hifi.demux.${B}--${B}.jasmine.pbmm2.combined.bedGraph -loj |
     awk -v I="${i}" '{OFS="\t"}{if ($10=="."){print $1,$2,$3,"random_NORPT","random_NORPT","0",I} else {print $1,$2,$3,"random_NORPT","random_NORPT",$10,I}}' \
      > PBmCG.10k_random-positions.NO-REPEATS.bslop25bp.${i}.bed



    grep $i ./SomaticInsertions_CURATED-list.no-comm.mod.tsv | grep CONFIRMED | cut -f1,2,3,11,13 |
        bedtools slop -i stdin -g /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/2_hifiasm_ragtag/Tsu-0/Tsu-0.Chr_scaffolds.v2.genome -b 25 |
        bedtools intersect -a stdin -b ./PB_mCG/tsumet.merged.ccs.hifi.demux.${B}--${B}.jasmine.pbmm2.combined.bedGraph -loj |
        awk '{OFS="\t"}{if ($9=="."){print $1,$2,$3,$4,"TE","0",$5} else {print $1,$2,$3,$4,"TE",$9,$5}}' \
        > PBmCG.insertion-sites.bslop25bp.${i}.bed


    # exit 1

done
