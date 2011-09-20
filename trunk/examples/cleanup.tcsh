#!/bin/tcsh

##
## USAGE: cleanup.tcsh
##

if ($#argv != 0) then
  grep '^##' $0
  exit
endif

rm -f \
  rnaseq.density.txt \
  TSS.10kb.bed \
  profile.txt \
  TSS.10kb.bed \
  heatmap.txt \
  densities.wig \
  peaks.bed \
  tss.val \
  tss.val+go \
  peaks.enriched.go.in.tss

gzip *.bed

