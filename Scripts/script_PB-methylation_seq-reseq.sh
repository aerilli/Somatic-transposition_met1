#!/bin/bash


# amovilli on 24-04-25

set -xv

OUTDIR=$1 #/ebio/scratch/amovilli/mC
movie_seq=$2 # /ebio/seq/lra/runs/runs/64079/r64079_20221220_102429/1_A01/m64079_221220_112036.subreads.bam
movie_reseq=$3 #/ebio/seq/lra/runs/runs/64079/r64079_20240212_104439/1_A01/m64079_240212_113350.subreads.bam
OUT=$4 #tsumet
REF=$5 #/ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/2_hifiasm_ragtag/Tsu-0/Tsu-0.Chr_scaffolds.v2.fa
CORES=$6
listreads=/ebio/abt6_projects/met1_somatic_transpositions/data/3_somatic-events/r_mutations/SomaticInsertions_CURATED-list.no-comm.mod.24-06-14.READ-ID-LIST.txt  # list explaining somatic insertions



ccsout_seq=$OUT.seq.ccs.bam
ccsout_reseq=$OUT.reseq.ccs.bam
ccsout_merged=$OUT.merged.ccs.bam
ccsout_merged_hifi=$OUT.merged.ccs.hifi.bam
ccsout_merged_hifi_demux=$OUT.merged.ccs.hifi.demux.bam
ccsout_merged_hifi_4jasmine=$OUT.merged.ccs.hifi.demux.bc*--bc*.bam


source activate /ebio/abt6_projects/Ath_HiFi_met1/conda/envs/pbcoretools_2023

mkdir -p $OUTDIR

cd $OUTDIR

ccs --all  --hifi-kinetics -j $CORES --report-file $OUTDIR/${ccsout_seq%.bam}.ccs_report.txt \
 $movie_seq \
 $OUTDIR/$ccsout_seq

ccs --all  --hifi-kinetics -j $CORES --report-file $OUTDIR/${ccsout_reseq%.bam}.ccs_report.txt \
 $movie_reseq \
 $OUTDIR/$ccsout_reseq


pbmerge -o $OUTDIR/$ccsout_merged  $OUTDIR/$ccsout_seq $OUTDIR/$ccsout_reseq



echo "Indexing BAM ..."
pbindex  $OUTDIR/$ccsout_merged


extracthifi $OUTDIR/$ccsout_merged \
 $OUTDIR/$ccsout_merged_hifi


lima $OUTDIR/$ccsout_merged_hifi /ebio/abt6_projects/Ath_HiFi_met1/data/assembly/barcodes/barcodes.fasta $OUTDIR/$ccsout_merged_hifi_demux \
--same --ccs --min-score 70 --min-scoring-regions 2 --min-ref-span 0.8 --peek-guess --split-bam-named -j $CORES

pbmm2 index $REF $REF.mmi

> $OUTDIR/reads-explaining-somatic-ins.bam
for i in  $ccsout_merged_hifi_4jasmine ; do

	# only reads not explaining insertions
	samtools view -@ $CORES -h -b -q 50 -N $listreads $OUTDIR/$i > $OUTDIR/${i%.bam}.filtered.bam
	# only reads not explaining insertions
	samtools view -@ $CORES -h -q 50 $OUTDIR/$i | grep -v -f $listreads  | samtools view -@ $CORES -h -b  >> $OUTDIR/reads-explaining-somatic-ins.bam

	jasmine -j $CORES  $OUTDIR/${i%.bam}.filtered.bam \
	 $OUTDIR/${i%.bam}.jasmine.bam

	pbmm2 align -j $CORES --preset CCS --sort \
		$REF.mmi \
		$OUTDIR/${i%.bam}.jasmine.bam \
		$OUTDIR/${i%.bam}.jasmine.pbmm2.bam

	/ebio/abt6_projects/Ath_HiFi_met1/conda/software/pb-CpG-tools/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores \
	 --bam $OUTDIR/${i%.bam}.jasmine.pbmm2.bam \
	 --output-prefix $OUTDIR/${i%.bam}.jasmine.pbmm2 \
	 --model /ebio/abt6_projects/Ath_HiFi_met1/conda/software/pb-CpG-tools/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite \
	 --threads $CORES

done


i=reads-explaining-somatic-ins.merged.ccs.hifi.demux.ALL.bam

jasmine -j $CORES  $OUTDIR/reads-explaining-somatic-ins.bam \
	 $OUTDIR/${i%.bam}.jasmine.bam

pbmm2 align -j $CORES --preset CCS --sort \
	$REF.mmi \
	$OUTDIR/${i%.bam}.jasmine.bam \
	$OUTDIR/${i%.bam}.jasmine.pbmm2.bam

/ebio/abt6_projects/Ath_HiFi_met1/conda/software/pb-CpG-tools/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores \
	--bam $OUTDIR/${i%.bam}.jasmine.pbmm2.bam \
	--output-prefix $OUTDIR/${i%.bam}.jasmine.pbmm2 \
	--model /ebio/abt6_projects/Ath_HiFi_met1/conda/software/pb-CpG-tools/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite \
	--threads $CORES




conda deactivate