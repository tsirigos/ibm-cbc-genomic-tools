#!/bin/tcsh -f

##
## USAGE: example05.tcsh
##

if ($#argv != 0) then
  grep '^##' $0
  exit
endif

set win_size = 500
set win_dist = 25
set min_reads = 20
set pval = 1e-05

gunzip chipseq.bed.gz
gunzip control.bed.gz 
gunzip genome.bed.gz

genomic_scans peaks -v -cmp -w $win_size -d $win_dist -min $min_reads -pval $pval -g genome.bed chipseq.bed control.bed \
  | genomic_regions bed -t "peaks" -c '0,150,0' >! peaks.bed

echo
echo '*********************************'
echo "ChIP-seq peaks in BED format: "
echo '*********************************'
head peaks.bed
echo ...

