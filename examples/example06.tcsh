#!/bin/tcsh -f

##
## USAGE: example06.tcsh 
##

if ($#argv != 0) then
  grep '^##' $0
  exit
endif

cat peaks.bed | genomic_overlaps -density -v -i TSS.10kb.bed | grep -v '^track' >! tss.val

cat tss.val | tr ':' '\t' | cut -f1,3 | sort | uniq | vectors -merge -n 6 | vectors -max -n 6 | join -a1 -t '	' - gene.go >! tss.val+go

permutation_test -v -h -S n -a -p 10000 -q 0.05 tss.val+go >! peaks.enriched.go.in.tss

echo
echo "******************************************"
echo "Enriched Gene Ontology terms in peaks: "
echo "******************************************"
head peaks.enriched.go.in.tss
echo ...

