#!/bin/tcsh -f

##
## USAGE: example02.tcsh
##

if ($#argv != 0) then
  grep '^##' $0
  exit
endif

gunzip genes.bed.gz
gunzip chipseq.bed.gz

cat genes.bed | genomic_regions pos -op 5p | genomic_regions shiftp -5p -10000 -3p +10000 >! TSS.10kb.bed
cat chipseq.bed | genomic_overlaps offset -v -i -op 5p -a TSS.10kb.bed | vectors -hist -n 6 -b 100 >! profile.txt

echo
echo '****************************************'
echo "ChIP-seq read profile in TSS regions: "
echo '****************************************'
head profile.txt
echo ...

