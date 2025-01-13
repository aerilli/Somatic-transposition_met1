#!/bin/bash

# liftoff gene annotation of Tsu-0 reference genome using TAIR10 gene annotation

# Variables
REF=/Data/Tsu-0.reference.fa
generef=./Data/TAIR10.gene-anno.gff
generefSingleISO=./Data/TAIR10.gene-anno.single-iso.gff
outliftover=/Data/Tsu-0.gene-anno.gff

grep -Ev '\.[02-9]+' $generef > $generefSingleISO

# Virtual env activation:
# containing liftoff v1.6.3 (https://github.com/agshumate/Liftoff) 
## source activate /ebio/abt6_projects/Ath_HiFi_met1/conda/envs/Annotators

# Using primary isoforms for liftoff
liftoff -g $generefSingleISO \
 -o $outliftover \
 -f ./Data/list.liftoff.txt \
 -p 64  -copies   $REF \
 ./Data/TAIR10.fa

# Virtual env deactivation
## conda deactivate