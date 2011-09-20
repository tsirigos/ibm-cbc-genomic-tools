#!/bin/tcsh -f

##
## USAGE: example01.tcsh
##

if ($#argv != 0) then
  grep '^##' $0
  exit
endif

gunzip exons.bed.gz
gunzip rnaseq.reads.bed.gz 

cat rnaseq.reads.bed | genomic_overlaps density -v exons.bed >! rnaseq.density.txt

echo
echo '********************************'
echo "RNA-seq densities in exons: "
echo '********************************'
head rnaseq.density.txt
echo ...

