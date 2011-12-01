#!/bin/tcsh -f

##
## USAGE: example03.tcsh
##

if ($#argv != 0) then
  grep '^##' $0
  exit
endif

gunzip genes.bed.gz
gunzip chipseq.bed.gz

cat genes.bed | genomic_regions pos -op 5p | genomic_regions shiftp -5p -10000 -3p +10000 | genomic_regions fix >! TSS.10kb.bed
cat chipseq.bed | genomic_overlaps offset -v -i -op 5p -a TSS.10kb.bed | sort | vectors -merge | vectors -bins -n 6 -b 200 >! heatmap.txt

echo
echo '******************************************************'
echo "Heatmap of ChIP-seq densities for gene TSS regions: "
echo '******************************************************'
head heatmap.txt
echo ...

