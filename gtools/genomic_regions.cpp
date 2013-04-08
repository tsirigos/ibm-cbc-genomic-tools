//
// Copyright (c) 2011 IBM Corporation. 
// All rights reserved. This program and the accompanying materials are made available under the terms of the GNU General Public License v2.0 
// which accompanies this distribution, and is available at http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <map>
#include <string>
#include <list>
#include <vector>
#include <iostream>
#include <algorithm>
#include "core.h"
#include "genomic_intervals.h"
using namespace std;



//---------------------------------------------------------------------------------//
// Global variables & constants                                                    //
//---------------------------------------------------------------------------------//

const string PROGRAM = "genomic_regions";
const long int BUFFER_SIZE = 10000;
gsl_rng *RANDOM_GENERATOR;




  
//---------------------------------------------------------------------------------//
// COMMAND-LINE OPTIONS                                                            //
//---------------------------------------------------------------------------------//

bool VERBOSE;
bool HELP;
bool PRINT_COMPACT;
char *CHROMOSOME_FASTA_FILE = NULL;
char *CHROMOSOME_MAP_DIR = NULL;
char *CHROMOSOME_MAP_NAME = NULL;
char *GENOME_REG_FILE = NULL;
char *REF_REG_FILE = NULL;
char *CHROMOSOME_NAMES = NULL;
bool REPLACE;
bool IGNORE;
bool UNIQ;
long int SHIFT_START;
long int SHIFT_STOP;
long int SHIFT_5PRIME;
long int SHIFT_3PRIME;
bool DENSITY_EXACT;
char *DIST_OP1;
char *DIST_OP2;
long int WIN_STEP;
long int WIN_SIZE;
bool WIN_MAX_LABEL_VALUE;
bool WIN_IGNORE_REVERSE_STRAND;
char WIN_PREPROCESS;
long int WIN_MIN_READS;
char *INDEX_FILE;
bool OVERLAP_EXACT;
bool OVERLAP_ONE;
bool OVERLAP_NEGATIVE;
bool OVERLAP_EXTRACT;
char *OVERLAP_OFFSET;
bool SELECT_FIRST;
bool SELECT_LAST;
bool SELECT_5P;
bool SELECT_3P;
char *STRAND_OP;
bool STRAND_SORTED;
char *POSITION_OP;
long int POSITION_SHIFT;
char *PATTERN;
bool CLUSTER_MERGE;
char *OFFSET_OP;
bool OFFSET_FRACTION;
char *TRACK_TITLE;
char *TRACK_COLOR;
double TRACK_SCALE;
long int TRACK_SPAN;
char *TRACK_POSITION;
char *TRACK_OPTIONS;
bool HEADER;
bool SUMMARY;
bool CONVERT_CHROMOSOME;
bool SORTED_BY_STRAND;
long int LINK_MAX_DIFFERENCE;
char *LINK_LABEL_FUNC;
int BIN_BITS;
bool IGNORE_STRAND;
long int UPSTREAM_MIN_DISTANCE;
long int UPSTREAM_MAX_DISTANCE;
bool SPLIT_BY_CHROM_AND_STRAND;





/*
//-----OBSOLETE/Usage----------
//
void Usage(const string prog)
{
  fprintf(stderr, "\n");
  fprintf(stderr, "USAGE: \n");
  fprintf(stderr, "  %s OPERATION [OPTIONS] <REGION-SET>\n", prog.c_str()); 
  fprintf(stderr, "\n");
  fprintf(stderr, "Line-based operations for reg, bed and fasta (when applicable) input formats: \n");
  fprintf(stderr, "  * align     align sequences to reference genome (line-based)\n");  
  fprintf(stderr, "  * bed       convert input regions to BED format (line-based)\n");  
  fprintf(stderr, "  * bounds    check interval start/stop positions against genome chromosome regions and fix all rogue intervals (line-based)\n");
  fprintf(stderr, "  * center    print center interval\n");
  fprintf(stderr, "  * connect   connect intervals from minimum start to maximum stop\n");
  fprintf(stderr, "  * diff      compute the difference between successive intervals\n");
  fprintf(stderr, "  * dist      compute distances between successive region intervals\n");
  fprintf(stderr, "  * divide    divide intervals in the middle\n");
  fprintf(stderr, "  * fix       remove rogue intervals (i.e. start<1 or start>stop).\n");
  fprintf(stderr, "  * int       compute the intersection of input intervals\n");
  fprintf(stderr, "  * n         compute total interval length (including possible overlaps)\n");
  fprintf(stderr, "  * pos       modify interval start/stop positions\n");
  fprintf(stderr, "  * reg       convert to REG format\n");
  fprintf(stderr, "  * rnd       randomize intervals across entire genome\n");
  fprintf(stderr, "  * select    selects a subset of intervals according to their relative start positions\n");
  fprintf(stderr, "  * shift     shift interval start/stop positions\n");
  fprintf(stderr, "  * shuffle   shuffle intervals within given reference region\n");
  fprintf(stderr, "  * sort      sort intervals\n");
  fprintf(stderr, "  * split     split regions into their intervals which are printed on separate lines\n");
  fprintf(stderr, "  * strand    modify interval strand information\n");
  fprintf(stderr, "  * union     compute the interval union\n");
  fprintf(stderr, "  * wig       convert input reg to UCSC wiggle format\n");  
  fprintf(stderr, "  * win       create new invervals by sliding windows\n");
  fprintf(stderr, "  * x         extract corresponding sequences from DNA files\n");
#ifdef FULL_VERSION
  fprintf(stderr, "  * bedgraph  convert input reg to UCSC bedgraph format\n");  
  fprintf(stderr, "  * len       compute total net interval length excluding Ns in the corresponding sequence (but including possible overlaps)\n");
  fprintf(stderr, "  * noN       remove N-regions\n");
  fprintf(stderr, "  * offset    convert to offset format\n");
  fprintf(stderr, "  * rc        correct interval start/stop positions on reverse strand\n");
#endif
  fprintf(stderr, "\n");

  fprintf(stderr, "File-based operations for reg, bed and fasta (when applicable) input formats: \n");
  fprintf(stderr, "  * gdist     compute distances of successive regions\n");
  fprintf(stderr, "  * inv       invert regions given the genome chromosomal boundaries\n");
  fprintf(stderr, "  * link      link consecutive regions to produce a non-overlapping set\n");
  fprintf(stderr, "  * test      test whether genomic regions are sorted and non-overlapping\n");
#ifdef FULL_VERSION
  fprintf(stderr, "  * annot     annotate sorted regions\n");  
  fprintf(stderr, "  * cluster   cluster regions based on overlaps\n");
  fprintf(stderr, "  * pairs     print all pairs of consecutive regions\n");
  fprintf(stderr, "  * partition partition overlapping regions into non-overlapping parts\n");
  fprintf(stderr, "  * rev       reverse the order of reverse-strand regions\n");
  fprintf(stderr, "  * scan      scan input intervals in sliding windows and print read count distribution\n");		// NOTE: move this to scan_reads
#endif
  fprintf(stderr, "\n");

#ifdef FULL_VERSION
  fprintf(stderr, "Line-based operations for fasta format only: \n");
  fprintf(stderr, "  * nogap     remove gaps from sequences\n");
  fprintf(stderr, "  * search    search short pattern\n");
  fprintf(stderr, "  * verify    verify sequences\n");
  fprintf(stderr, "\n");
#endif
}
*/




//-------InitCmdLine-----------
//
CmdLineWithOperations *InitCmdLine(int argc, char *argv[], int *next_arg)
{
  // initialize
  CmdLineWithOperations *cmd_line = new CmdLineWithOperations(); 
  cmd_line->SetProgramName(PROGRAM,VERSION);

  // set operations
  cmd_line->AddOperation("align", "[OPTIONS] <REGION-SET>", \
  "Prints alignments of input sequences given the reference genome sequence.", \
  "* Input formats: SAM\n\
  * Operands: region, sequence, reference sequence\n\
  * Region requirements: chromosome/strand-compatible, sorted, non-overlapping\n\
  * Region-set requirements: none"\
  );
  
  cmd_line->AddOperation("annotator", "[OPTIONS] <REGION-SET>", \
  "Creates upstream regions given a set of reference regions.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operands: interval\n\
  * Region requirements: single-interval regions\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("bed", "[OPTIONS] <REGION-SET>", \
  "Converts input regions into BED format.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operand: region\n\
  * Region requirements: chromosome/strand-compatible, sorted, non-overlapping\n\
  * Region-set requirements: none"\
  ); 

  cmd_line->AddOperation("bounds", "[OPTIONS] <REGION-SET>", \
  "Trims start/stop positions given the chromosome bounds and removes invalid intervals.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operand: interval\n\
  * Region requirements: none\n\
  * Region-set requirements: none"\
  );
  
  cmd_line->AddOperation("center", "[OPTIONS] <REGION-SET>", \
  "Prints center interval.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operand: interval\n\
  * Region requirements: none\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("chrom", "[OPTIONS] <REGION-SET>", \
  "Modifies chromosome names.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operand: interval\n\
  * Region requirements: none\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("connect", "[OPTIONS] <REGION-SET>", \
  "Connects intervals from minimum start to maximum stop position.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operand: region\n\
  * Region requirements: chromosome/strand-compatible\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("diff", "[OPTIONS] <REGION-SET>", \
  "Computes the difference between successive intervals.", \
  "* Input formats: REG, BED, SAM\n\
  * Operand: interval pair\n\
  * Region requirements: chromosome/strand-compatible, sorted, non-overlapping\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("dist", "[OPTIONS] <REGION-SET>", \
  "Computes distances between pairs of successive intervals.", \
  "* Input formats: REG, BED, SAM\n\
  * Operand: interval pair\n\
  * Region requirements: none\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("divide", "[OPTIONS] <REGION-SET>", \
  "Divides intervals in the middle.", \
  "* Input formats: REG, BED\n\
  * Operand: interval\n\
  * Region requirements: none\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("fix", "[OPTIONS] <REGION-SET>", \
  "Removes invalid intervals, i.e. start<1 or start>stop.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operand: interval\n\
  * Region requirements: none\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("int", "[OPTIONS] <REGION-SET>", 
  "Computes the intersection of input intervals.", \
  "* Input formats: REG\n\
  * Operand: region\n\
  * Region requirements: chromosome/strand-compatible\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("n", "[OPTIONS] <REGION-SET>", \
  "Computes sum of lengths of intervals in a region (including possible overlaps).", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operand: region\n\
  * Region requirements: none\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("pos", "[OPTIONS] <REGION-SET>", \
  "Modifies interval start/stop positions.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operand: interval\n\
  * Region requirements: none\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("reg", "[OPTIONS] <REGION-SET>", \
  "Converts to REG format.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operand: interval\n\
  * Region requirements: chromosome/strand-compatible if option -c is set\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("rnd", "[OPTIONS] <REGION-SET>", \
  "Randomizes region across entire genome (relative interval distances are preserved).", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operand: region\n\
  * Region requirements: chromosome/strand-compatible, sorted, non-overlapping\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("search", "[OPTIONS] <REGION-SET>", \
  "Searches sequences for a short pattern.", \
  "* Input formats: SEQ\n\
  * Operand: sequence\n\
  * Region requirements: single-interval\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("select", "[OPTIONS] <REGION-SET>", \
  "Selects a subset of intervals according to their position in the region.", \
  "* Input formats: REG, BED, SAM\n\
  * Operand: region\n\
  * Region requirements: chromosome/strand-compatible\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("shift", "[OPTIONS] <REGION-SET>", \
  "Shifts interval start/stop positions.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operand: interval\n\
  * Region requirements: none\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("shiftp", "[OPTIONS] <REGION-SET>", \
  "Shifts interval 5\'/3\' positions.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operand: interval\n\
  * Region requirements: none\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("shuffle", "[OPTIONS] <REGION-SET>", 
  "Shuffles intervals within given reference region-set.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operand: interval\n\
  * Region requirements: single-interval\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("sort", "[OPTIONS] <REGION-SET>", 
  "Sorts region intervals by start position.", \
  "* Input formats: REG\n\
  * Operand: region\n\
  * Region requirements: none\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("split", "[OPTIONS] <REGION-SET>", 
  "Splits regions into their intervals which are printed on separate lines.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operand: region\n\
  * Region requirements: none\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("strand", "[OPTIONS] <REGION-SET>", \
  "Modifies interval strand information.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operand: interval\n\
  * Region requirements: chromosome/strand-compatible\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("union", "[OPTIONS] <REGION-SET>", \
  "Computes union of region intervals.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operand: region\n\
  * Region requirements: chromosome/strand-compatible\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("wig", "[OPTIONS] <REGION-SET>", \
  "Converts to UCSC wiggle format.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operand: region\n\
  * Region requirements: single-interval\n\
  * Region-set requirements: sorted by chromosome and start position"\
  );
  
  cmd_line->AddOperation("win", "[OPTIONS] <REGION-SET>", \
  "Creates new invervals by sliding windows.", \
  "* Input formats: REG, BED\n\
  * Operand: region\n\
  * Region requirements: single-interval\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("x", "[OPTIONS] <REGION-SET>", \
  "Extracts region sequence from DNA.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operand: region\n\
  * Region requirements: chromosome/strand-compatible\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("gdist", "[OPTIONS] <REGION-SET>", \
  "Computes distances of successive regions.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operand: region-pair\n\
  * Region requirements: single-interval\n\
  * Region-set requirements: sorted by chromosome/strand/start"\
  );

  cmd_line->AddOperation("gsort", "[OPTIONS] <REGION-SET>", \
  "Global sort: sorts the entire region set.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operand: region-set\n\
  * Region requirements: none\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("inv", "[OPTIONS] <REGION-SET>", \
  "Inverts regions given the genome chromosomal boundaries.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operand: region-set\n\
  * Region requirements: single-interval\n\
  * Region-set requirements: sorted by chromosome/strand/start"\
  );

  cmd_line->AddOperation("link", "[OPTIONS] <REGION-SET>", \
  "Links consecutive regions to produce a non-overlapping set.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operand: region-set\n\
  * Region requirements: single-interval\n\
  * Region-set requirements: sorted by chromosome/(strand)/start"\
  );
  
  cmd_line->AddOperation("test", "[OPTIONS] <REGION-SET>", \
  "Tests whether genomic regions are sorted and non-overlapping.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operand: region\n\
  * Region requirements: chromosome/strand-compatible, sorted, non-overlapping\n\
  * Region-set requirements: sorted by chromosome/(strand)/start"\
  );

  int min_files = 0;
  
  // print summary of operations
  if (argc<2) { 
    cmd_line->OperationSummary("OPERATION [OPTIONS] <REGION-SET>", \
    "Performs operations on genomic intervals and genomic regions. The input regions can be supplied as a file name or from the standard input. For detailed description and list of options choose an operation and use the --help option."); 
    delete cmd_line; 
    exit(1); 
  }

  // current operation
  string op = argv[1];
  if (op[0]=='-') op = op.substr(1);		// ensure compatibility with previous version
  cmd_line->SetCurrentOperation(op);
  
  // common options
  cmd_line->AddOption("--help", &HELP, false, "help");
  cmd_line->AddOption("-h", &HELP, false, "help");
  cmd_line->AddOption("-v", &VERBOSE, false, "verbose mode");

  // Main options
  if (op=="align") {
    cmd_line->AddOption("-Q", &CHROMOSOME_FASTA_FILE, "", "FASTA file containing all chromosomes (overrides options -D and -q)");
    cmd_line->AddOption("-D", &CHROMOSOME_MAP_DIR, ".", "chromosome DNA file and map directory");
    cmd_line->AddOption("-q", &CHROMOSOME_MAP_NAME, "chromosome.map", "sequence map file");
  }
  else if (op=="annotator") {
	cmd_line->AddOption("-g", &GENOME_REG_FILE, "", "genome region-set file");
    cmd_line->AddOption("-i", &IGNORE_STRAND, false, "ignore strand while finding overlaps");
    cmd_line->AddOption("--upstream-max", &UPSTREAM_MAX_DISTANCE, 1000000, "maximum allowed upstream region size ");
    cmd_line->AddOption("--upstream-min", &UPSTREAM_MIN_DISTANCE, 10000, "minimum allowed upstream region size (subject to genomic bounds)");
  }
  else if (op=="bed") {
    cmd_line->AddOption("-t", &TRACK_TITLE, "", "title");
    cmd_line->AddOption("-c", &TRACK_COLOR, "", "color");
    cmd_line->AddOption("-p", &TRACK_POSITION, "", "browser position");
    cmd_line->AddOption("-chr", &CONVERT_CHROMOSOME, false, "convert chromosome names from ENSEMBL to UCSC");
  }
  else if (op=="bedgraph") {
    cmd_line->AddOption("-t", &TRACK_TITLE, "", "title");
    cmd_line->AddOption("-c", &TRACK_COLOR, "0,0,0", "color");
    cmd_line->AddOption("-p", &TRACK_POSITION, "", "browser position");
    cmd_line->AddOption("-chr", &CONVERT_CHROMOSOME, false, "convert chromosome names from ENSEMBL to UCSC");
  }
  else if (op=="bounds") {
    cmd_line->AddOption("-g", &GENOME_REG_FILE, "", "genome region-set file");
  }
  else if (op=="center") {

  }
  else if (op=="chrom") {
    cmd_line->AddOption("-c", &CHROMOSOME_NAMES, "", "chromosome name conversion table");
  }
  else if (op=="cluster") {
    cmd_line->AddOption("-m", &CLUSTER_MERGE, false, "merge regions in each cluster");
  }
  else if (op=="connect") {

  }
  else if (op=="diff") {

  }
  else if (op=="dist") {
    cmd_line->AddOption("-op1", &DIST_OP1, "1", "reference point of 1st interval in pair (1=start, 2=stop, 5p=5'-end, 3p=3'-end)");
    cmd_line->AddOption("-op2", &DIST_OP2, "1", "reference point of 2nd interval in pair (1=start, 2=stop, 5p=5'-end, 3p=3'-end)");
  }
  else if (op=="divide") {

  }
  else if (op=="fix") {

  }
  else if (op=="gdist") {
    cmd_line->AddOption("-op1", &DIST_OP1, "1", "reference point of 1st interval in pair (1=start, 2=stop, 5p=5'-end, 3p=3'-end)");
    cmd_line->AddOption("-op2", &DIST_OP2, "1", "reference point of 2nd interval in pair (1=start, 2=stop, 5p=5'-end, 3p=3'-end)");
  }
  else if (op=="gsort") {
    cmd_line->AddOption("-s", &SORTED_BY_STRAND, false, "sort by strand in addition to chromosome and start position");
    cmd_line->AddOption("-b", &BIN_BITS, 12, "bucket size (in bits) used in bucket sort");
  }
  else if (op=="int") {

  }
  else if (op=="inv") {
    cmd_line->AddOption("-g", &GENOME_REG_FILE, "", "genome region-set file");
  }
  else if (op=="labels") {

  }
  else if (op=="len") {
    cmd_line->AddOption("-D", &CHROMOSOME_MAP_DIR, ".", "directory where chromosome sequences are located");
    cmd_line->AddOption("-q", &CHROMOSOME_MAP_NAME, "chromosome.map", "sequence map file");
  }
  else if (op=="link") {
    cmd_line->AddOption("-s", &SORTED_BY_STRAND, false, "input regions are sorted by strand");
    cmd_line->AddOption("-d", &LINK_MAX_DIFFERENCE, 0, "maximum difference between successive regions");
    cmd_line->AddOption("--label-func", &LINK_LABEL_FUNC, "", "label function = {min,max,sum,%c}, where %c is used as delimiter");
  }
  else if (op=="n") {

  }
  else if (op=="nogap") {

  }
  else if (op=="noN") {
    cmd_line->AddOption("-D", &CHROMOSOME_MAP_DIR, ".", "directory where chromosome sequences are located");
    cmd_line->AddOption("-q", &CHROMOSOME_MAP_NAME, "chromosome.map", "sequence map file");
  }
  else if (op=="overlap") {
    cmd_line->AddOption("-I", &INDEX_FILE, "", "sorted index region/sequence file");
    cmd_line->AddOption("-x", &OVERLAP_EXACT, false, "compute exact overlap boundaries");
    cmd_line->AddOption("-1", &OVERLAP_ONE, false, "print output in one line omitting index labels");
    cmd_line->AddOption("-neg", &OVERLAP_NEGATIVE, false, "consider also negative strand positions for index reg+ files");
    cmd_line->AddOption("-q", &OVERLAP_EXTRACT, false, "extract overlapping sequences from region labels");
    cmd_line->AddOption("-off", &OVERLAP_OFFSET, "", "offset reference point (1=start, 2=stop, 5p=5'-end, 3p=3'-end)");
  }
  else if (op=="offset") {
    cmd_line->AddOption("-op", &OFFSET_OP, "1", "reference point (1=start, 2=stop, 5p=5'-end, 3p=3'-end)");
    cmd_line->AddOption("-a", &OFFSET_FRACTION, false, "print distances as a fraction of total size");
  }
  else if (op=="pairs") {

  }
  else if (op=="partition") {

  }
  else if (op=="pos") {
    cmd_line->AddOption("-op", &POSITION_OP, "1", "position operation (1=start, 2=stop, 5p=5'-end, 3p=3'-end, c=center)");
    cmd_line->AddOption("-c", &POSITION_SHIFT, 0, "position shift");
  }
  else if (op=="reg") {
    cmd_line->AddOption("-c", &PRINT_COMPACT, false, "print in compact starts/ends format");
  }
  else if (op=="rc") {
    cmd_line->AddOption("-g", &GENOME_REG_FILE, "", "genome region-set file");
  }
  else if (op=="rev") {

  }
  else if (op=="rnd") {
    cmd_line->AddOption("-g", &GENOME_REG_FILE, "", "genome region-set file");
  }
  else if (op=="search") {
    cmd_line->AddOption("-p", &PATTERN, "C|G", "pattern");
    cmd_line->AddOption("-H", &HEADER, false, "print header");
    cmd_line->AddOption("-s", &SUMMARY, false, "print summary statistics only");
  }
  else if (op=="select") {
    cmd_line->AddOption("-first", &SELECT_FIRST, false, "select the first interval");
    cmd_line->AddOption("-last", &SELECT_LAST, false, "select the last interval");
    cmd_line->AddOption("-5p", &SELECT_5P, false, "select the first from the 5\' end");
    cmd_line->AddOption("-3p", &SELECT_3P, false, "select the first from the 3\' end");
  }
  else if (op=="shift") {
    cmd_line->AddOption("-start", &SHIFT_START, 0, "shift distance for start position");
    cmd_line->AddOption("-stop", &SHIFT_STOP, 0, "shift distance for stop position");
  }
  else if (op=="shiftp") {
    cmd_line->AddOption("-5p", &SHIFT_5PRIME, 0, "shift distance for 5\' end");
    cmd_line->AddOption("-3p", &SHIFT_3PRIME, 0, "shift distance for 3\' end");
  }
  else if (op=="shuffle") {
    cmd_line->AddOption("-r", &REF_REG_FILE, "", "reference region file");
  }
  else if (op=="split") {
    cmd_line->AddOption("--by-chrstrand", &SPLIT_BY_CHROM_AND_STRAND, false, "sort and split by chromosome and strand (REG format only)");
  }
  else if (op=="strand") {
    cmd_line->AddOption("-op", &STRAND_OP, "+", "strand operation (+=positive, -=negative, r=reverse)");
    cmd_line->AddOption("-s", &STRAND_SORTED, false, "returns a sorted region set (only works if input is sorted)");
  }
  else if (op=="sort") {

  }
  else if (op=="test") {
    cmd_line->AddOption("-s", &SORTED_BY_STRAND, false, "input regions are sorted by strand");
  }
  else if (op=="union") {

  }
  else if (op=="verify") {
    cmd_line->AddOption("-D", &CHROMOSOME_MAP_DIR, ".", "chromosome DNA file and map directory");
    cmd_line->AddOption("-q", &CHROMOSOME_MAP_NAME, "chromosome.map", "sequence map file");
    cmd_line->AddOption("-i", &IGNORE, false, "ignore errors");
  }
  else if (op=="wig") {
    cmd_line->AddOption("-t", &TRACK_TITLE, "", "title");
    cmd_line->AddOption("-c", &TRACK_COLOR, "200,0,0", "color");
    cmd_line->AddOption("-scale", &TRACK_SCALE, 1.0, "multiply values by constant factor");
    cmd_line->AddOption("-s", &TRACK_SPAN, 1, "span");
    cmd_line->AddOption("-p", &TRACK_POSITION, "", "browser position");
    cmd_line->AddOption("-o", &TRACK_OPTIONS, "", "track type options");
    cmd_line->AddOption("-chr", &CONVERT_CHROMOSOME, false, "convert chromosome names from ENSEMBL to UCSC");
  }
  else if (op=="win") {
    cmd_line->AddOption("-s", &WIN_SIZE, 1, "window size");
    cmd_line->AddOption("-d", &WIN_STEP, 1, "window distance");
  }
  else if (op=="x") {
    cmd_line->AddOption("-Q", &CHROMOSOME_FASTA_FILE, "", "FASTA file containing all chromosomes (overrides options -D and -q)");
    cmd_line->AddOption("-D", &CHROMOSOME_MAP_DIR, ".", "chromosome DNA file and map directory");
    cmd_line->AddOption("-q", &CHROMOSOME_MAP_NAME, "chromosome.map", "sequence map file");
    cmd_line->AddOption("-r", &REPLACE, false, "replace 'N' characters with 'a'");
    cmd_line->AddOption("-i", &IGNORE, false, "ignore boundary errors");
  }
  else {
    cerr << "Unknown operation '" << op << "'!\n";
    delete cmd_line;
    exit(1);
  }

  // process command line
  *next_arg = cmd_line->Read(argv+1,argc-1) + 1;
  if (HELP||(argc-*next_arg<min_files)) { cmd_line->OperationUsage(); delete cmd_line; exit(1); }
  
  return cmd_line;
}







//---------------------------------------------------------------------------------//
// MAIN	                                                                           //
//---------------------------------------------------------------------------------//
int main(int argc, char* argv[]) 
{
  // process command-line arguments
  int next_arg;
  CmdLineWithOperations *cmd_line = InitCmdLine(argc,argv,&next_arg); 
  char *REG_FILE = next_arg==argc ? NULL : argv[next_arg];
  _MESSAGES_ = VERBOSE;

  // initilize random generator
  RANDOM_GENERATOR = InitRandomGenerator(getpid()+time(NULL));

  // setup bounds & chromosomes if necessary
  StringLIntMap *bounds = GENOME_REG_FILE==NULL?NULL:ReadBounds(GENOME_REG_FILE,VERBOSE);
  Chromosomes *C = NULL;
  if (CHROMOSOME_MAP_DIR!=NULL) {
    if (strlen(CHROMOSOME_FASTA_FILE)>0) C = new Chromosomes(CHROMOSOME_FASTA_FILE,VERBOSE);
    else if (strlen(CHROMOSOME_MAP_DIR)>0) C = new Chromosomes(CHROMOSOME_MAP_DIR,CHROMOSOME_MAP_NAME); 
  }


  // open region sets
  bool hide_header = ((cmd_line->current_cmd_operation=="bed")||(cmd_line->current_cmd_operation=="reg"))?true:false;
  bool load_in_memory = ((cmd_line->current_cmd_operation=="gsort")||(cmd_line->current_cmd_operation=="annotator"))?true:false;
  GenomicRegionSet RegSet(REG_FILE,BUFFER_SIZE,VERBOSE,load_in_memory,hide_header);

  //----------------------------------------------//
  //  Line-based (horizontal) operations          //
  //----------------------------------------------//
  if (cmd_line->current_cmd_operation=="align") RegSet.RunAlign(C);
  else if (cmd_line->current_cmd_operation=="bed") RegSet.RunConvertToBED(TRACK_TITLE,TRACK_COLOR,TRACK_POSITION,CONVERT_CHROMOSOME);
  else if (cmd_line->current_cmd_operation=="bounds") RegSet.RunBounds(bounds);
  else if (cmd_line->current_cmd_operation=="center") RegSet.RunCenter();
  else if (cmd_line->current_cmd_operation=="chrom") RegSet.RunModifyChromosomeNames(CHROMOSOME_NAMES);
  else if (cmd_line->current_cmd_operation=="connect") RegSet.RunConnect();
  else if (cmd_line->current_cmd_operation=="diff") RegSet.RunDiff();
  else if (cmd_line->current_cmd_operation=="dist") RegSet.RunCalcDistances(DIST_OP1,DIST_OP2);
  else if (cmd_line->current_cmd_operation=="divide") RegSet.RunDivide();
  else if (cmd_line->current_cmd_operation=="fix") RegSet.RunFix();
  else if (cmd_line->current_cmd_operation=="int") RegSet.RunIntersection();
  else if (cmd_line->current_cmd_operation=="n") RegSet.RunSize();
  else if (cmd_line->current_cmd_operation=="pos") RegSet.RunModifyPos(POSITION_OP,POSITION_SHIFT);
  else if (cmd_line->current_cmd_operation=="reg") RegSet.RunConvertToREG(PRINT_COMPACT);
  else if (cmd_line->current_cmd_operation=="rnd") RegSet.RunRandomize(RANDOM_GENERATOR,bounds);
  else if (cmd_line->current_cmd_operation=="search") RegSet.PrintSearch(PATTERN,HEADER,SUMMARY);
  else if (cmd_line->current_cmd_operation=="select") RegSet.RunSelect(SELECT_FIRST,SELECT_LAST,SELECT_5P,SELECT_3P);
  else if (cmd_line->current_cmd_operation=="shift") RegSet.RunShiftPos(SHIFT_START,SHIFT_STOP,false);
  else if (cmd_line->current_cmd_operation=="shiftp") RegSet.RunShiftPos(SHIFT_5PRIME,SHIFT_3PRIME,true);
  else if (cmd_line->current_cmd_operation=="shuffle") RegSet.RunShuffle(RANDOM_GENERATOR,REF_REG_FILE);
  else if (cmd_line->current_cmd_operation=="sort") RegSet.RunSort();
  else if (cmd_line->current_cmd_operation=="split") RegSet.RunSplit(SPLIT_BY_CHROM_AND_STRAND);
  else if (cmd_line->current_cmd_operation=="strand") RegSet.RunModifyStrand(STRAND_OP,STRAND_SORTED);
  else if (cmd_line->current_cmd_operation=="union") RegSet.RunUnion();
  else if (cmd_line->current_cmd_operation=="wig") RegSet.RunConvertToWIG(TRACK_TITLE,TRACK_COLOR,TRACK_POSITION,TRACK_OPTIONS,TRACK_SCALE,TRACK_SPAN,CONVERT_CHROMOSOME);
  else if (cmd_line->current_cmd_operation=="win") RegSet.RunSlidingWindows(WIN_STEP,WIN_SIZE);
  else if (cmd_line->current_cmd_operation=="x") RegSet.RunExtractSeq(C,REPLACE);
  
  //----------------------------------------------//
  //  File-based (vertical) operations            //
  //----------------------------------------------//
  else if (cmd_line->current_cmd_operation=="annotator") RegSet.RunGlobalCreateAnnotator(bounds,IGNORE_STRAND,UPSTREAM_MAX_DISTANCE,UPSTREAM_MIN_DISTANCE);
  else if (cmd_line->current_cmd_operation=="gdist") RegSet.RunGlobalCalcDistances(DIST_OP1,DIST_OP2);
  else if (cmd_line->current_cmd_operation=="gsort") RegSet.RunGlobalSort(SORTED_BY_STRAND,BIN_BITS);
  else if (cmd_line->current_cmd_operation=="inv") RegSet.RunGlobalInvert(bounds);
  else if (cmd_line->current_cmd_operation=="link") RegSet.RunGlobalLink(SORTED_BY_STRAND,LINK_MAX_DIFFERENCE,LINK_LABEL_FUNC);
  else if (cmd_line->current_cmd_operation=="test") RegSet.RunGlobalTest(SORTED_BY_STRAND);




#ifdef FULL_VERSION
  //----------------------------------------------//
  //  UNDER DEVELOPMENT/TESTING                   //
  //----------------------------------------------//
  else if (cmd_line->current_cmd_operation=="bedgraph") RegSet.PrintBEDGraphFormat(TRACK_TITLE,TRACK_COLOR,TRACK_POSITION,CONVERT_CHROMOSOME);
  else if (cmd_line->current_cmd_operation=="cluster") RegSet.RunGlobalCluster(CLUSTER_MERGE);
  else if (cmd_line->current_cmd_operation=="len") RegSet.PrintSeqLength(C);
  //else if (cmd_line->current_cmd_operation=="nogap") RegSet.PrintNoGaps();
  else if (cmd_line->current_cmd_operation=="noN") { RegSet.PrintRemoveN(C); }
  //else if (cmd_line->current_cmd_operation=="offset") RegSet.PrintOffsetFormat(OFFSET_OP,OFFSET_FRACTION);
  //else if (cmd_line->current_cmd_operation=="overlap") { GenomicRegionSet INDEX(INDEX_FILE,BUFFER_SIZE,VERBOSE); RegSet.PrintOverlaps(INDEX,OVERLAP_EXACT,OVERLAP_ONE,OVERLAP_EXTRACT,OVERLAP_NEGATIVE,OVERLAP_OFFSET); }
  else if (cmd_line->current_cmd_operation=="partition") RegSet.RunGlobalPartition();
  else if (cmd_line->current_cmd_operation=="rc") RegSet.PrintReversePos(bounds);
  else if (cmd_line->current_cmd_operation=="rev") RegSet.RunGlobalReverseOrder(); 
  else if (cmd_line->current_cmd_operation=="verify") RegSet.PrintVerifySeq(C,IGNORE);
#endif

  else { cerr << "Unknown method '" << cmd_line->current_cmd_operation << "'!\n"; delete cmd_line; exit(1); }

  // clean up
  delete cmd_line;
  if (bounds!=NULL) delete bounds;
  if (C!=NULL) delete C;

  
  return 0;
}



