#!/bin/tcsh -f

##
## USAGE: example04.tcsh
##

if ($#argv != 0) then
  grep '^##' $0
  exit
endif

set win_size = 500
set win_dist = 25
set min_reads = 20

gunzip chipseq.bed.gz
gunzip genome.bed.gz

cat chipseq.bed | genomic_scans counts -v -min $min_reads -w $win_size -d $win_dist -g genome.bed | genomic_regions center \
  | genomic_regions wig -t "densities" -s $win_dist -c '0,0,150' >! densities.wig

echo 
echo '****************************************************'
echo "Genome-wide ChIP-seq densities in WIGGLE format: "
echo '****************************************************'
head densities.wig
echo ...

