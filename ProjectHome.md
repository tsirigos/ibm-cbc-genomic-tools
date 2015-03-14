<p>
<h1><a>Contents</a></h1>
<ul>
<li><a href='#Introduction.md'>Introduction</a></li>
<li><a href='#Latest_news.md'>Latest news</a></li>
<li><a href='#License.md'>License</a></li>
<li><a href='#Installation.md'>Installation</a></li>
<li><a href='#Documentation.md'>Documentation</a></li>
<li><a href='#Examples.md'>Examples</a></li>
<li><a href='#Citation.md'>Citation</a></li>
<li><a href='#References.md'>References</a></li>
<li><a href='#Contact.md'>Contact</a></li>


<p>
<h1>Introduction</h1>
<p>GenomicTools is a flexible computational platform for the analysis and manipulation of high-throughput sequencing data such as RNA-seq and ChIP-seq. GenomicTools implements a variety of mathematical operations between sets of genomic regions thereby enabling the prototyping of computational pipelines that can address a wide spectrum of tasks from preprocessing and quality control to meta-analyses. More specifically, the user can easily create average read profiles across transcriptional start sites or enhancer sites, quickly prototype customized peak discovery methods for ChIP-seq experiments, perform genome-wide statistical tests such as enrichment analyses, design controls via appropriate randomization schemes, among other applications. In addition to enabling rapid prototyping, the GenomicTools platform is designed to analyze large-datasets in a single-pass fashion in order to minimize memory and intermediate file requirements. Finally, the GenomicTools platform supports the widely used BED format to facilitate visualization as well as integration with existing platforms and pipelines such as <a href='http://galaxy.psu.edu/'>Galaxy</a> or <a href='http://www.bioconductor.org/'>BioConductor</a>.<br>
<br>
<br>
<br>
<br>
<p>
<h1>Latest news</h1>


<p>-- April 2, 2014: version 2.8.0 is now available:<br>
<ul>
<li> new tool: gtools_hic for Hi-C sequencing data analysis (alignment and filtering)</li>
<li> new option: genomic_overlaps annotate --query-op {center,overlap}</li>
<li> undocumented feature: in genomic_apps heatmap reordered heatmap rows top-to-bottom instead of bottom-to-top</li>
<li> new option: genomic_scans counts/peaks, genomic_apps profile/heatmap, genomic_overlaps count/coverage/density: --max-label-value to replace --labels-as-values/-n/-val </li>
</ul>
</p>

<p>-- January 15, 2013: version 2.7.2 is now available:<br>
<ul>
<li> minor bug fix: genomic_regions search was crashing when input sequences were fed from stdin </li>
<li> new option: genomic_regions wig -scale </li>
<li> new option: genomic_overlaps annotate --proximal-dist/--print-header </li>
<li> output format: added header and row names to the output dat file of genomic_apps heatmap </li>
</ul>
</p>

<p>-- October 9, 2012: version 2.7.0 is now available:<br>
<ul>
<li> minor bug fix: a small number of operations was rejecting BAM files with header </li>
<li> new operation: genomic_overlaps annotate </li>
<li> new operation (beta): genomic_overlaps bin </li>
<li> new option: genomic_regions split --by-chrstrand </li>
<li> new option: genomic_regions link --label-func </li>
<li> new option: genomic_apps peakdiff -fold </li>
</ul>
</p>

<p>-- September 12, 2012: version 2.6.0 is now available. It includes various improvements in the peakdiff pipeline and its output files. </p>

<p>-- July 1, 2012: version 2.5.0 is now available. In summary, new options were implemented for 'genomic_overlaps offset', 'genomic_apps profile' and 'genomic_apps heatmap'. Additionally, a new operation, 'genomic_regions gsort' can be used to globally sort genomic region sets faster than UNIX sort. Finally, 'genomic_apps peakdiff' makes use of the UnsortedGenomicRegionSetScanner class, therefore it does not require sorted genomic regions as inputs. </p>

<p>-- June 8, 2012: version 2.4.0 contains improvements and minor bug fixes related to the GenomicRegionSetScanner class, affecting genomic_scans (operations counts/peaks) and genomic_apps (operation peakdiff). A new class, UnsortedGenomicRegionSetScanner, can now be used to scan unsorted sets of genomic regions. Currently, this class can be invoked only in 'genomic_scans counts', but in the future we will enable it for 'genomic_scans peaks' as well. </p>

<p>-- June 5, 2012: version 2.3.2 is released which allows non-standard intervals (start<=0) and ignores invalid intervals (start>stop or stop<=0) when computing overlaps using the algorithm for unsorted intervals (the algorithm for sorted intervals is able to handle these cases without changes).</p>

<p>-- April 5, 2012: version 2.3.0(beta) now accepts gzipped and BAM files as input using the <a href='http://www.cs.unc.edu/Research/compgeom/gzstream/'>gzstream library</a> and code from <a href='http://samtools.sourceforge.net/'>SAMtools</a>. Note: this functionality is only available for files; if input is fed through the standard input, then simply use 'gunzip' or 'samtools view' respectively. </p>

<p>-- March 27, 2012: version 2.2.0 introduces an application for differential peak discovery between two samples. This is an adaptation of the method described in <a href='http://www.nature.com/nm/journal/vaop/ncurrent/full/nm.2651.html'>Ntziachristos et al.</a>

<img src='http://ibm-cbc-genomic-tools.googlecode.com/files/peakdiff.snapshot.v2.jpg' />

<p>-- March 19, 2012: version 2.1.0 introduces <i>genomic_apps</i>, a series of operations whose goal is to create plots as output (in tiff format). Currently implemented operations can be used to create tif images of heatmap representations and read profiles for sequencing data. The basic implementation idea is to first compute the necessary data using the GenomicTools API and then call an R script. Users are also given the option to modify the script so as to obtain customized behavior. R version 2.14 (or later) is required. To get a quick summary of the options and a few examples, simply type 'genomic_apps heatmap' or 'genomic_apps profile'. Source code and/or linux precompiled binaries can be downloaded from the 'Featured Downloads' list. The following snapshot was created using <i>genomic_apps</i> (sequencing data was obtained from <a href='http://www.pnas.org/content/108/36/14908.short'>Wang et al.</a>).<br>
<br>
<img src='http://ibm-cbc-genomic-tools.googlecode.com/files/snapshot.v1.jpg' />

<p>-- November 23, 2011: published as  <a href='http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btr646?ijkey=3dsKvzqvhmEq9WQ&keytype=ref'>Applications Note in Bioinformatics</a>.<br>
<br>
<p>-- October 24, 2011: uploaded logo created by <a href='http://www.behance.net/AntoninoSerafino'>Antonino Serafino</a>.<br>
<br>
<p>-- September 27, 2011: version v2.0.0 is released; source, documentation and pre-compiled executables are now available.<br>
<br>
<p>-- September 23, 2011: uploaded user's manual for v2.0.0.<br>
<br>
<p>-- September 12, 2011: version v2.0.0 is under development; source code is now on the Google Code source repository <a href='http://code.google.com/p/ibm-cbc-genomic-tools/source/checkout'>http://code.google.com/p/ibm-cbc-genomic-tools/source/checkout</a>.<br>
<br>
<p>-- August 19, 2011: version v1.3.0(beta) was released (allows GFF/SAM formats as input).<br>
<br>
<p>-- July 15, 2011: version v1.2.2 was released (genomic_overlaps bug fix and improved performance).<br>
<br>
<p>-- April 19, 2011: version v1.0 of GenomicTools was released.<br>
<br>
<p>
<h1>License</h1>
<p>The GenomicTools source code is released under the GNU General Public License.<br>
<br>
<br>
<br>
<br>
<br>
<p>
<h1>Installation</h1>
<p>Precompiled binaries for linux/x86_64 are now available in the svn repository (bin/linux-x86_64). To compile GenomicTools from source, first download the source code from the <a href='https://code.google.com/p/ibm-cbc-genomic-tools/source/checkout'>svn repository</a>, and then follow these steps (see also requirements below):<br>
<pre><code>cd &lt;genomic-tools-directory&gt;<br>
make<br>
sudo cp bin/* /usr/local/bin<br>
</code></pre>
<b>Requirements</b>:<br>
<p>-- GenomicTools require installation of the <a href='http://www.gnu.org/software/gsl/'>GNU Scientific Library (GSL)</a> before compilation of the source code.</p>
<p>-- genomic_apps requires R 2.14 or later and the following R libraries: gplots and BioConductor's preprocessCore.</p>
<b>Notes</b>:<br>
<p>-- we suggest that GSL be installed by compiling the source code, otherwise it needs to be installed with its headers (gsl-devel rather than plain gsl; also on some systems you will need gslcblas-devel)</p>
<p>-- on Mac OS X there maybe an issue with the header files from GSL: if they are not installed in /usr/include and are instead in the /opt/local/include, the compiler may not be able to locate them (as mentioned above, our suggestion is to install GSL from source instead of using a repository such as MacPorts)</p>




<p>
<h1>Documentation</h1>
<p>The GenomicTools platform was implemented in C/C++. We provide the following documentation and tutorials:<br>
<ul>
<li><a href='http://ibm-cbc-genomic-tools.googlecode.com/files/genomic-tools-v2.0.2-user-manual.pdf'>User's manual</a> for using the command-line version of the tools</li>
<li>Documentation of the C++ classes of the GenomicTools API is available with the source code distribution (subdirectory "doc") and online <a href='http://www.cs.nyu.edu/~tsirigos/genomic-tools-doc/'>here</a>.</li>
</ul>




<p>
<h1>Examples</h1>
The examples are now part of the repository (subdirectory "examples"). A detailed tutorial is available in the User's Manual in the chapter titled "Case study: a simple ChIPseq pipeline".<br>
<br>
<br>
<br>
<br>
<br>
<p>
<h1>Citation</h1>
<a href='http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btr646?ijkey=3dsKvzqvhmEq9WQ&keytype=ref'>GenomicTools: a computational platform for developing high-throughput analytics in genomics</a>, Tsirigos A, Haiminen N, Bilal E, Utro F, Bioinformatics 2011.<br>
<br>
<br>
<br>
<br>
<p>
<h1>References</h1>
GenomicTools was presented as a poster at the CSHL Biology of the Genomes meeting (May 2011). Part of the GenomicTools suite is the genomic_overlaps tool for computing overlaps between genomic regions. When compared to BEDTools and Bioconductor it yields multi-fold improvements both in time and memory performance (download the flyer <a href='http://ibm-cbc-genomic-tools.googlecode.com/files/GenomicTools.CSHL.flyer.pdf'>here</a>). The algorithm operates on sorted inputs, scans the files sequentially and computes all overlaps essentially using a mergesort algorithm adapted to handle intervals instead of numbers.<br>
<br>
GenomicTools has been used for extensive genome-wide computations in the following publications:<br>
<ul>
<li>Ntziachristos et al.<br>
<a href='http://www.nature.com/nm/journal/vaop/ncurrent/full/nm.2651.html'>
Genetic inactivation of the polycomb repressive complex 2 in T cell acute lymphoblastic leukemia</a> - <i>Nature Medicine</i>, January 2012.</li>
<li>Laurent L, Wong E, Li G, Huynh T, Tsirigos A, Ong CT, Low HM, Kin Sung KW, Rigoutsos I, Loring J, Wei CL.<br>
<a href='http://www.ncbi.nlm.nih.gov/pubmed/20133333'>
Dynamic Changes in the Human Methylome During Differentiation</a> - <i>Genome Research</i>, February 2010.</li>
<li>Tsirigos A, Rigoutsos I.<br>
<a href='http://www.ncbi.nlm.nih.gov/pubmed/20019790'>
Alu and B1 repeats have been selectively retained in upstream and intronic regions of genes of specific functional classes</a> - <i>PLoS Computational Biology</i>, December 2009.</li>
<li>Weiss A, Charbonnier E, Ellertsdóttir E, Tsirigos A, Wolf C, Schuh R, Pyrowolakis G, Affolter M.<br>
<a href='http://www.ncbi.nlm.nih.gov/pubmed/20010841'>
A conserved activation element in BMP signaling during Drosophila development</a> - <i>Nature Structural & Molecular Biology</i>, December 2009.</li>
<li>Ochoa-Espinosa A, Yu D, Tsirigos A, Struffi P, Small S. <a href='http://www.ncbi.nlm.nih.gov/pubmed/19237583'>Anterior-posterior positional information in the absence of a strong Bicoid gradient</a> - <i>PNAS</i>, February 2009.</li>
<li>Tsirigos A, Rigoutsos I. <a href='http://www.ncbi.nlm.nih.gov/pubmed/18450818'>Human and mouse introns are linked to the same processes and functions through each genome's most frequent non-conserved motifs</a> - <i>Nucleic Acids Research</i>, May 2008.</li>
</ul>




<p>
<h1>Contact</h1>
<p>Contact Aris Tsirigos at tsirigos@gmail.com. For research projects that use GenomicTools <a href='http://www.cs.nyu.edu/~tsirigos'>visit my website</a>.<br>
<br>
<br>
<br>
<p>Copyright (c) 2011 IBM Corporation. All rights reserved.<br>
<br>
<br>
<br>
<p>Disclaimer: the GenomicTools software is distributed as is with no guarantee whatsoever. It is strongly advised to independently verify the results.