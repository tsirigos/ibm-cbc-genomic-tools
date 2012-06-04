
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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <map>
#include <string>
#include <list>
#include <queue>
#include <vector>
#include <iostream>
#include <algorithm>
#include <sstream>
#include "core.h"
#include "genomic_intervals.h"
#include "genomic_apps.rscripts.h"


using namespace std;


//---------------------------------------------------------------------------------//
// Global variables                                                                //
//---------------------------------------------------------------------------------//

const string PROGRAM = "genomic_apps";
const long int BUFFER_SIZE = 10000;
gsl_rng *RANDOM_GENERATOR;






  
//---------------------------------------------------------------------------------//
// Command-line options                                                            //
//---------------------------------------------------------------------------------//

bool VERBOSE;
bool HELP;

char *RSCRIPT_INPUT_FILE_NAME;
bool REUSE;
char *OUT_PREFIX; 
bool IGNORE_STRAND;
char *SHIFT;
char *LEGEND;
char *COLORS;
char *TITLE, *XLABEL, *YLABEL;
char *IMAGE_TYPE; 
char *IMAGE_SIZE;
int IMAGE_RESOLUTION;
bool NORMALIZE;

char *GENOME_REG_FILE;
long int WIN_SIZE, WIN_DIST;
double PVALUE_CUTOFF;
float FDR;
float OUTLIER_PROB;
unsigned long int PSEUDOCOUNT;
char *LABELS;






//-------InitCmdLine-----------
//
CmdLineWithOperations *InitCmdLine(int argc, char *argv[], int *next_arg)
{
  // initialize
  CmdLineWithOperations *cmd_line = new CmdLineWithOperations(); 
  cmd_line->SetProgramName(PROGRAM,VERSION);

  // set operations
  cmd_line->AddOperation("heatmap", "[OPTIONS] COMMA-SEPARATED-SIGNAL-REG-FILES REFERENCE-REG-FILE", \
  "Create heatmap profile of signal region in reference region.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operand: interval\n\
  * Region requirements: single-interval\n\
  * Region-set requirements: none", \
  "1. genomic_apps heatmap -v -o heatmap1 -i -colors red Notch1.bed tss.bed" \
  );

  cmd_line->AddOperation("peakdiff", "[OPTIONS] SAMPLE1-FILES SAMPLE2-FILES [SAMPLE1-CONTROL-FILES SAMPLE2-CONTROL-FILES]", \
"Finds differences in peak levels between two samples (two replicates per sample are required).", \
"This operation finds differences in read density levels of genomic windows between two samples. The algorithm used here is an adaptation of the method presented in Ntziachristos et al. (Nature Medicine, Feb 2012) and follows these steps: \n\
  1. identifies enriched genomic windows (i.e. peaks) in the input sample files based on the binomial distribution against random background or supplied control files (see options -i, -w, -d and -pval)\n\
  2. normalizes the genomic windows densities inside the identified peaks across replicates and samples using quantile normalization\n\
  3. for each sample, checks replicate reproducibility and removes outliers by first fitting a smooth line between two replicates and then applying a cutoff at the residual values whose distribution is modeled as a normal distribution (see option -outliers)\n\
  4. computes significant fold changes between samples for the identified peaks (see options -pseudo and -fdr)\n\
\n\
REQUIREMENTS: \n\
  * Input files: comma-separated, no spaces!\n\
  * Input formats: REG, GFF, BED, SAM\n\
  * Operand: interval\n\
  * Region requirements: single-interval\n\
  * Region-set requirements: sorted by chromosome/strand/start\n\
\n\
EXAMPLES: \n\
  1. genomic_apps peakdiff -v -o Notch1.vs.Input -g genome.bed Notch1.r1.sorted.bed,Notch1.r2.sorted.bed Input.r1.sorted.bed,Input.r2.sorted.bed\n\
  2. genomic_apps peakdiff -v -o Notch1_cancer.vs.Notch1_normal -g genome.bed Notch1_cancer.r1.sorted.bed,Notch1_cancer.r2.sorted.bed Notch1_normal.r1.sorted.bed,Notch1_normal.r2.sorted.bed\n\
  3. genomic_apps peakdiff -v -o Notch1_cancer.vs.Notch1_normal.with_controls -g genome.bed Notch1_cancer.r1.sorted.bed,Notch1_cancer.r2.sorted.bed Notch1_normal.r1.sorted.bed,Notch1_normal.r2.sorted.bed Input_cancer.r1.sorted.bed,Input_cancer.r2.sorted.bed Input_normal.r1.sorted.bed,Input_normal.r2.sorted.bed" \
  );

  cmd_line->AddOperation("profile", "[OPTIONS] COMMA-SEPARATED-SIGNAL-REG-FILES COMMA-SEPARATED-REFERENCE-REG-FILES", \
  "Create profile(s) of signal regions in reference regions.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operand: interval\n\
  * Region requirements: single-interval\n\
  * Region-set requirements: none", \
  "1. genomic_apps profile -v -o profile1 -i -colors red,green -legend 'Notch1 in TSSs','H3K4me3 in TSSs' Notch1.bed,H3K4me3.bed tss.bed\n\
  2. genomic_apps profile -v -o profile2 -i -colors red,green -legend 'H3K4me1 in Notch1','H3K4me1 in Myc' H3K4me1.bed Notch1.bed,Myc.bed\n\
  3. genomic_apps profile -v -o profile2 -reuse -i -xlab 'distance from center (nt)' -ylab 'number of reads' -title 'read profiles' -colors red,green -legend 'H3K4me1 in Notch1','H3K4me1 in Myc' H3K4me1.bed Notch1.bed,Myc.bed" \
  );

  if (argc<2) { cmd_line->OperationSummary("OPERATION [OPTIONS] INPUT-FILES","Implements common genomics applications."); exit(1); }

  // current operation
  string op = argv[1];
  cmd_line->SetCurrentOperation(op);
  
  // common options
  cmd_line->AddOption("--help", &HELP, false, "help");
  cmd_line->AddOption("-h", &HELP, false, "help");
  cmd_line->AddOption("-v", &VERBOSE, false, "verbose mode");
  
  // Main options
  int n_args;
  if (cmd_line->current_cmd_operation=="heatmap") {
    n_args = 2;
    cmd_line->AddOption("-reuse", &REUSE, false, "reuse histogram data; update paramaters only");	
    cmd_line->AddOption("-R", &RSCRIPT_INPUT_FILE_NAME, "", "R script file to use (not required)");	
    cmd_line->AddOption("-o", &OUT_PREFIX, "", "prefix for output files");
    cmd_line->AddOption("-i", &IGNORE_STRAND, false, "ignore strand while finding overlaps");
    cmd_line->AddOption("-shift", &SHIFT, "5000,5000", "comma-separated upstream/downstream distances from reference center");
    cmd_line->AddOption("-colors", &COLORS, "", "comma-separated colors for heatmap pixels");
    cmd_line->AddOption("-title", &TITLE, "", "heatmap comma-separated titles");
    cmd_line->AddOption("-xlab", &XLABEL, "", "heatmap x-axis label");
    cmd_line->AddOption("-ylab", &YLABEL, "", "heatmap y-axis label");
    cmd_line->AddOption("-itype", &IMAGE_TYPE, "pdf", "image format type {pdf,tif}");
    cmd_line->AddOption("-isize", &IMAGE_SIZE, "2000,4000", "comma-separated image dimensions (for tif format only)");
    cmd_line->AddOption("-ires", &IMAGE_RESOLUTION, 600, "image resolution in dpi (for tif format only)");
  }
  else if (cmd_line->current_cmd_operation=="peakdiff") {
    n_args = 2;
    cmd_line->AddOption("-reuse", &REUSE, false, "reuse histogram data; update paramaters only");	
    cmd_line->AddOption("-R", &RSCRIPT_INPUT_FILE_NAME, "", "R script file to use (not required)");	
    cmd_line->AddOption("-o", &OUT_PREFIX, "", "prefix for output files");
    cmd_line->AddOption("-i", &IGNORE_STRAND, false, "ignore strand while finding overlaps");
    cmd_line->AddOption("-g", &GENOME_REG_FILE, "", "genome region file");
    cmd_line->AddOption("-w", &WIN_SIZE, 500, "window size (must be a multiple of window distance)");
    cmd_line->AddOption("-d", &WIN_DIST, 100, "window distance");
    cmd_line->AddOption("-pval", &PVALUE_CUTOFF, 1.0e-05, "p-value cutoff for calling significant windows");
    cmd_line->AddOption("-outliers", &OUTLIER_PROB, 0.01, "probability cutoff for residuals in outlier detection between replicates");
    cmd_line->AddOption("-pseudo", &PSEUDOCOUNT, 1, "pseudocount to be added to window count for fold-change computations");
    cmd_line->AddOption("-fdr", &FDR, 0.05, "false discover rate for differential peak discovery");
    cmd_line->AddOption("-labels", &LABELS, "", "comma-separated sample labels");
    cmd_line->AddOption("-itype", &IMAGE_TYPE, "pdf", "image format type {pdf,tif}");
    cmd_line->AddOption("-isize", &IMAGE_SIZE, "2000,2000", "comma-separated image dimensions (for tif format only)");
    cmd_line->AddOption("-ires", &IMAGE_RESOLUTION, 300, "image resolution in dpi (for tif format only)");
  }
  else if (cmd_line->current_cmd_operation=="profile") {
    n_args = 2;
    cmd_line->AddOption("-reuse", &REUSE, false, "reuse histogram data; update paramaters only");	
    cmd_line->AddOption("-R", &RSCRIPT_INPUT_FILE_NAME, "", "R script file to use (not required)");	
    cmd_line->AddOption("-o", &OUT_PREFIX, "", "prefix for output files");
    cmd_line->AddOption("-i", &IGNORE_STRAND, false, "ignore strand while finding overlaps");
    cmd_line->AddOption("-norm", &NORMALIZE, false, "normalize against total number of reads and number of reference regions");	
    cmd_line->AddOption("-shift", &SHIFT, "5000,5000", "comma-separated upstream/downstream distances from reference center");
    cmd_line->AddOption("-legend", &LEGEND, "", "comma-separated legend labels for line plot");
    cmd_line->AddOption("-colors", &COLORS, "", "comma-separated colors for line plot");
    cmd_line->AddOption("-title", &TITLE, "", "plot title");
    cmd_line->AddOption("-xlab", &XLABEL, "", "plot x-axis label");
    cmd_line->AddOption("-ylab", &YLABEL, "", "plot y-axis label");
    cmd_line->AddOption("-itype", &IMAGE_TYPE, "pdf", "image format type {pdf,tif}");
    cmd_line->AddOption("-isize", &IMAGE_SIZE, "2000,2000", "comma-separated image dimensions (for tif format only)");
    cmd_line->AddOption("-ires", &IMAGE_RESOLUTION, 300, "image resolution in dpi (for tif format only)");
  }
  else {
    cerr << "Unknown operation '" << op << "'!\n";
    delete cmd_line;
    exit(1);
  }

  // process command line
  *next_arg = cmd_line->Read(argv+1,argc-1) + 1;
  if (HELP||(argc-*next_arg<n_args)) { cmd_line->OperationUsage(); delete cmd_line; exit(1); }
  
  if ((strcmp(IMAGE_TYPE,"pdf")!=0)&&(strcmp(IMAGE_TYPE,"tif")!=0)) { fprintf(stderr, "Error: unsupported image format '%s'!\n", IMAGE_TYPE); exit(1); }
  
  return cmd_line;
}








//-------Tokenize-----------
//
char **Tokenize(char *inp, char delim, int *n_tokens)
{
  if ((inp==NULL)||(strlen(inp)==0)) { *n_tokens = 0; return NULL; }
  *n_tokens = CountTokens(inp,delim);
  char **token = new char*[*n_tokens];
  for (int k=0; k<*n_tokens; k++) token[k] = StrCopy(GetNextToken(&inp,delim)); 
  return token;
}



//-------CreateRscript-----------
//
void CreateRscript(const char *rscript_input_file_name, const char *rscript_template, string rscript_output_file_name)
{
  if (strlen(rscript_input_file_name)==0) {
    FILE *rscript_out_file = fopen(rscript_output_file_name.c_str(),"w");
	if (rscript_out_file==NULL) { fprintf(stderr, "Error: cannot open file '%s' for writing!\n", rscript_output_file_name.c_str()); exit(1); }
	fprintf(rscript_out_file, "%s", rscript_template);
    fclose(rscript_out_file);
  }
  else {
    FILE *rscript_in_file = fopen(rscript_input_file_name,"r");
	if (rscript_in_file==NULL) { fprintf(stderr, "Error: R script file '%s' not found!\n", rscript_input_file_name); exit(1); }
	char *rscript = LoadFile(rscript_in_file);
	fclose(rscript_in_file);
    FILE *rscript_out_file = fopen(rscript_output_file_name.c_str(),"w");
	if (rscript_out_file==NULL) { fprintf(stderr, "Error: cannot open file '%s' for writing!\n", rscript_output_file_name.c_str()); exit(1); }
	fprintf(rscript_out_file, "%s", rscript);
    fclose(rscript_out_file);
  }
}




//-------PrintLogFile-----------
//
void PrintLogFile(string log_file_name)
{
  fprintf(stderr, "\n***********************************************\n"); 
  fprintf(stderr, "\n  R   S C R I P T   L O G   F I L E            \n"); 
  fprintf(stderr, "\n***********************************************\n"); 
  FileBufferText buffer(log_file_name.c_str());
  for (char *inp=buffer.Next(); inp!=NULL; inp=buffer.Next()) fprintf(stderr, "%s\n", inp);
}













//-----ScanReadFiles----------
//
void ScanReadFiles(char **signal_reg_file, int n_signal_files, char **ref_reg_file, int n_ref_files, \
                   char **signal_control_reg_file, int n_signal_control_files, char **ref_control_reg_file, int n_ref_control_files, \
				   char *genome_reg_file, const char *out_file_name)
{
  bool USE_COUNTS = false;
  char PREPROCESS = 'c';

  StringLIntMap *bounds = ReadBounds(genome_reg_file);
  unsigned long int effective_genome_size = CalcBoundSize(bounds);
  if (VERBOSE) fprintf(stderr, "* Effective genome size = %lu\n", effective_genome_size);

  // open signal files
  GenomicRegionSet **SignalRegSet = new GenomicRegionSet*[n_signal_files];
  GenomicRegionSetScanner **SignalRegScanner = new GenomicRegionSetScanner*[n_signal_files];
  double *p_signal = new double[n_signal_files];
  for (int s=0; s<n_signal_files; s++) {
    SignalRegSet[s] = new GenomicRegionSet(signal_reg_file[s],BUFFER_SIZE,VERBOSE,false,true);
    p_signal[s] = (double)CountLines(signal_reg_file[s],BUFFER_SIZE)/effective_genome_size;
    SignalRegSet[s]->Reset();
    if (VERBOSE) fprintf(stderr, "* Signal input file = %s; background probability = %.2e\n", signal_reg_file[s], p_signal[s]);
    SignalRegScanner[s] = new GenomicRegionSetScanner(SignalRegSet[s],bounds,WIN_DIST,WIN_SIZE,USE_COUNTS,IGNORE_STRAND,PREPROCESS);
  }
  
  // open signal control files
  GenomicRegionSet **SignalControlRegSet = n_signal_control_files>0?new GenomicRegionSet*[n_signal_control_files]:NULL;
  GenomicRegionSetScanner **SignalControlRegScanner = n_signal_control_files>0?new GenomicRegionSetScanner*[n_signal_control_files]:NULL;
  double *p_signal_control = n_signal_control_files>0?new double[n_signal_control_files]:NULL;
  for (int s=0; s<n_signal_control_files; s++) {
    SignalControlRegSet[s] = new GenomicRegionSet(signal_control_reg_file[s],BUFFER_SIZE,VERBOSE,false,true);
    p_signal_control[s] = (double)CountLines(signal_control_reg_file[s],BUFFER_SIZE)/effective_genome_size;
    SignalControlRegSet[s]->Reset();
    if (VERBOSE) fprintf(stderr, "* Signal control input file = %s; background probability = %.2e\n", signal_control_reg_file[s], p_signal_control[s]);
    SignalControlRegScanner[s] = new GenomicRegionSetScanner(SignalControlRegSet[s],bounds,WIN_DIST,WIN_SIZE,USE_COUNTS,IGNORE_STRAND,PREPROCESS);
  }
  
  // open reference files
  GenomicRegionSet **RefRegSet = new GenomicRegionSet*[n_ref_files];
  GenomicRegionSetScanner **RefRegScanner = new GenomicRegionSetScanner*[n_ref_files];
  double *p_ref = new double[n_ref_files];
  for (int r=0; r<n_ref_files; r++) {
    RefRegSet[r] = new GenomicRegionSet(ref_reg_file[r],BUFFER_SIZE,VERBOSE,false,true);
    p_ref[r] = (double)CountLines(ref_reg_file[r],BUFFER_SIZE)/effective_genome_size;
    RefRegSet[r]->Reset();
    if (VERBOSE) fprintf(stderr, "* Reference input file = %s; background probability = %.2e\n", ref_reg_file[r], p_ref[r]);
    RefRegScanner[r] = new GenomicRegionSetScanner(RefRegSet[r],bounds,WIN_DIST,WIN_SIZE,USE_COUNTS,IGNORE_STRAND,PREPROCESS);
  }
  
  // open reference control files
  GenomicRegionSet **RefControlRegSet = n_ref_control_files>0?new GenomicRegionSet*[n_ref_control_files]:NULL;
  GenomicRegionSetScanner **RefControlRegScanner = n_ref_control_files>0?new GenomicRegionSetScanner*[n_ref_control_files]:NULL;
  double *p_ref_control = n_ref_control_files>0?new double[n_ref_control_files]:NULL;
  for (int r=0; r<n_ref_control_files; r++) {
    RefControlRegSet[r] = new GenomicRegionSet(ref_control_reg_file[r],BUFFER_SIZE,VERBOSE,false,true);
    p_ref_control[r] = (double)CountLines(ref_control_reg_file[r],BUFFER_SIZE)/effective_genome_size;
    RefControlRegSet[r]->Reset();
    if (VERBOSE) fprintf(stderr, "* Reference control input file = %s; background probability = %.2e\n", ref_control_reg_file[r], p_ref_control[r]);
    RefControlRegScanner[r] = new GenomicRegionSetScanner(RefControlRegSet[r],bounds,WIN_DIST,WIN_SIZE,USE_COUNTS,IGNORE_STRAND,PREPROCESS);
  }
  
  
  // scanning   
  long int *signal_val = new long int[n_signal_files];
  long int *signal_control_val = n_signal_control_files>0?new long int[n_signal_control_files]:NULL;
  long int *ref_val = new long int[n_ref_files];
  long int *ref_control_val = n_ref_control_files>0?new long int[n_ref_control_files]:NULL;
  FILE *out_file = fopen(out_file_name,"w");
  if (out_file==NULL) { fprintf(stderr, "Error: cannot open file '%s' for writing!\n", out_file_name); exit(1); }
  Progress PRG("Scanning...",1);
  while (true) {
    for (int s=0; s<n_signal_files; s++) signal_val[s] = SignalRegScanner[s]->Next();
    for (int s=0; s<n_signal_control_files; s++) signal_control_val[s] = SignalControlRegScanner[s]->Next();
    for (int r=0; r<n_ref_files; r++) ref_val[r] = RefRegScanner[r]->Next();
    for (int r=0; r<n_ref_control_files; r++) ref_control_val[r] = RefControlRegScanner[r]->Next();
	if (signal_val[0]==-1) break;
	
	bool print_interval = false;
    if (n_signal_control_files>0) { 
      for (int s=0; s<n_signal_files; s++) 
	    if ((double)gsl_cdf_binomial_Q(signal_val[s],max(p_signal[s],(double)signal_control_val[s]/WIN_SIZE),WIN_SIZE)<=PVALUE_CUTOFF) { print_interval = true; break; }
	  if (print_interval==false) 
	    for (int r=0; r<n_ref_files; r++) 
		  if ((double)gsl_cdf_binomial_Q(ref_val[r],max(p_ref[r],(double)ref_control_val[r]/WIN_SIZE),WIN_SIZE)<=PVALUE_CUTOFF) { print_interval = true; break; }
	}
	else {
      for (int s=0; s<n_signal_files; s++) if ((double)gsl_cdf_binomial_Q(signal_val[s],p_signal[s],WIN_SIZE)<=PVALUE_CUTOFF) { print_interval = true; break; }
	  if (print_interval==false) for (int r=0; r<n_ref_files; r++) if ((double)gsl_cdf_binomial_Q(ref_val[r],p_ref[r],WIN_SIZE)<=PVALUE_CUTOFF) { print_interval = true; break; }
	}

    if (print_interval==true) {
      SignalRegScanner[0]->PrintInterval(out_file);
      for (int s=0; s<n_signal_files; s++) fprintf(out_file, "\t%ld", signal_val[s]);
      for (int r=0; r<n_ref_files; r++) fprintf(out_file, "\t%ld", ref_val[r]);
      for (int s=0; s<n_signal_control_files; s++) fprintf(out_file, "\t%ld", signal_control_val[s]);
      for (int r=0; r<n_ref_control_files; r++) fprintf(out_file, "\t%ld", ref_control_val[r]);
	  fprintf(out_file, "\n");
    }

    PRG.Check();
  }
  PRG.Done();
  fclose(out_file);
  
  // cleanup
  delete [] SignalRegSet;
  delete [] SignalRegScanner;
  delete p_signal;
  delete [] RefRegSet;
  delete [] RefRegScanner;
  delete p_ref;
  if (SignalControlRegSet) delete [] SignalControlRegSet;
  if (SignalControlRegScanner) delete [] SignalControlRegScanner;
  if (p_signal_control) delete p_signal_control;
  if (RefControlRegSet) delete [] RefControlRegSet;
  if (RefControlRegScanner) delete [] RefControlRegScanner;
  if (p_ref_control) delete p_ref_control;
  delete bounds;
  delete signal_val;
  if (signal_control_val) delete signal_control_val;
  delete ref_val;
  if (ref_control_val) delete ref_control_val;
}







// ************************************************************************************
//
//   M  A  I  N 
//
// ************************************************************************************
int main(int argc, char* argv[]) 
{
  // process command-line arguments
  int next_arg;
  CmdLineWithOperations *cmd_line = InitCmdLine(argc,argv,&next_arg); 
  _MESSAGES_ = VERBOSE;
  
  // initilize random generator
  RANDOM_GENERATOR = InitRandomGenerator(getpid()+time(NULL));

  
  
  // ************************************************************************************
  //
  //     H  E  A  T  M  A  P
  //
  // ************************************************************************************
  
  if (cmd_line->current_cmd_operation=="heatmap") {
    if (strlen(OUT_PREFIX)==0) { fprintf(stderr, "Error: prefix for output files must be specified using the -o option!\n"); exit(1); }
	
    char *SIGNAL_REG_FILES = argv[next_arg];
    char *REF_REG_FILES = argv[next_arg+1];
	
	// process input parameters
	char *shift = SHIFT;
	long int shift_upstream = atol(GetNextToken(&shift,','));
	long int shift_downstream = atol(GetNextToken(&shift,','));
	int n_signal_files;
	char **signal_reg_file = Tokenize(SIGNAL_REG_FILES,',',&n_signal_files);
	int n_ref_files;
	char **ref_reg_file = Tokenize(REF_REG_FILES,',',&n_ref_files);
	if (n_ref_files>1) { fprintf(stderr, "Error: only one reference file allowed!\n"); exit(1); }
	int n_colors = CountTokens(COLORS,',');
	if (n_colors!=n_signal_files) { fprintf(stderr, "Error: number of colors must match number of signal files!\n"); exit(1); }
    int n_titles = CountTokens(TITLE,',');;
	if (n_titles!=n_signal_files) { fprintf(stderr, "Error: number of titles must match number of signal files!\n"); exit(1); }
	
    // setup output file names
	string data_file_name = (string)OUT_PREFIX + (string)".dat";
	string param_file_name = (string)OUT_PREFIX + (string)".params";
	string rscript_file_name = (string)OUT_PREFIX + (string)".r";
	string image_file_name = (string)OUT_PREFIX + "." + (string)IMAGE_TYPE;
	string log_file_name = (string)OUT_PREFIX + (string)".log";

	// print parameters to file
	FILE *param_file = fopen(param_file_name.c_str(),"w");
    if (param_file==NULL) { fprintf(stderr, "Error: cannot open file '%s' for writing!\n", param_file_name.c_str()); exit(1); }
	fprintf(param_file, "%ld\n", shift_upstream);
	fprintf(param_file, "%ld\n", shift_downstream);
	fprintf(param_file, "%s\n", COLORS);
	fprintf(param_file, "%s\n", TITLE);
	fprintf(param_file, "%s\n", XLABEL);
	fprintf(param_file, "%s\n", YLABEL);
	fprintf(param_file, "%s\n", IMAGE_SIZE);
	fprintf(param_file, "%d\n", IMAGE_RESOLUTION);
	fprintf(param_file, "%d\n", n_signal_files);
	for (int i=0; i<argc; i++) fprintf(param_file, "%s%c", argv[i], i<argc-1?' ':'\n');
    fclose(param_file);

	// create R script
	CreateRscript(RSCRIPT_INPUT_FILE_NAME,RSCRIPT_TEMPLATE_HEATMAP,rscript_file_name);
	
    if (REUSE==false) {
      // set parameters
	  char *BIN_BITS = StrCopy("17,20,23,26");
	  bool MATCH_GAPS = false;
	  char *OFFSET_OP = StrCopy("5p");
      long int n_bins = (shift_upstream+shift_downstream)/50;
	  int n_bins_combine = 10; 
      long int bin_min = -shift_upstream;
      long int bin_max = shift_downstream;

      // open reference region set and shift 5prime position upstream/downstream
      GenomicRegionSet RefRegSet(ref_reg_file[0],BUFFER_SIZE,VERBOSE,true,true);
	  long int n_ref = RefRegSet.n_regions;
      GenomicRegion *r = RefRegSet.Get();
	  long int n_ref1 = r->n_line;
      while (r!=NULL) { 
        r->Connect();
	    r->Center();
        r->I.front()->ShiftPos(-shift_upstream,shift_downstream,true);
		r = RefRegSet.Next();
	  }

      // initialize bin array	  
	  int ***bins = new int**[n_signal_files];
	  for (int s=0; s<n_signal_files; s++) {
	    bins[s] = new int*[n_ref];
	    for (int r=0; r<n_ref; r++) { bins[s][r] = new int[n_bins]; for (int q=0; q<n_bins; q++) bins[s][r][q] = 0; }
	  }
	  
	  // process signal files
	  for (int n=0; n<n_signal_files; n++) {
	    if (VERBOSE) fprintf(stderr, "Creating heatmap of '%s' in '%s'...\n", signal_reg_file[n], ref_reg_file[0]); 
        GenomicRegionSet TestRegSet(signal_reg_file[n],BUFFER_SIZE,VERBOSE,false,true);
        UnsortedGenomicRegionSetOverlaps Overlaps(&TestRegSet,&RefRegSet,BIN_BITS);
        Progress PRG("Processing queries...",1);
        for (GenomicRegion *qreg=Overlaps.GetQuery(); Overlaps.Done()==false; qreg=Overlaps.NextQuery()) {
          for (GenomicRegion *ireg=Overlaps.GetOverlap(MATCH_GAPS,IGNORE_STRAND); ireg!=NULL; ireg=Overlaps.NextOverlap(MATCH_GAPS,IGNORE_STRAND)) {
            long int start_offset, stop_offset;
            qreg->I.front()->GetOffsetFrom(ireg->I.front(),OFFSET_OP,IGNORE_STRAND,&start_offset,&stop_offset);
            if (start_offset>stop_offset) { fprintf(stderr, "Error: start offset is greater than stop offest (this must be a bug)!\n"); exit(1); }
            long int x = (start_offset+stop_offset)/2-shift_upstream;
            double z = (double)(x-bin_min)/(bin_max-bin_min);
            if ((z>=0)&&(z<1)) bins[n][ireg->n_line-n_ref1][(int)(n_bins*z)]++;
          }
          PRG.Check();
        }
        PRG.Done();

	  }
	  
	  // store bin data in output file
	  FILE *data_file = fopen(data_file_name.c_str(),"w");
      if (data_file==NULL) { fprintf(stderr, "Error: cannot open file '%s' for writing!\n", data_file_name.c_str()); exit(1); }
      for (int r=0; r<n_ref; r++) {
        fprintf(data_file, "%d\t", r);
	    for (int s=0; s<n_signal_files; s++) {
		  int val = 0;
		  for (int q=0; q<n_bins_combine-1; q++) val += bins[s][r][q];
		  for (int q=0,qq=n_bins_combine-1; qq<n_bins; q++,qq++) {
		    val += bins[s][r][qq];
		    fprintf(data_file, "%d%s", val, qq!=n_bins-1?"\t":"");
		    val -= bins[s][r][q];
		  }
		  fprintf(data_file, "%c", s!=n_signal_files-1?'\t':'\n');
	    }
      }
	  fclose(data_file);
	  
	  // cleanup
      for (int s=0; s<n_signal_files; s++) delete [] bins[s];
      delete bins;	  
      delete BIN_BITS;
	  delete OFFSET_OP;
	}
	
	  
	// execute R script
	fprintf(stderr, "Running R script...\n");
	stringstream s;
	s << "R CMD BATCH '--args " << data_file_name << " " << param_file_name << " " << image_file_name << "' " << rscript_file_name << " " << log_file_name;
	system(s.str().c_str()); 
	PrintLogFile(log_file_name);
	
    // cleanup
	delete [] signal_reg_file;
	delete [] ref_reg_file;
	

  }




  // ************************************************************************************
  //
  //     P  E  A  K  D  I  F  F
  //
  // ************************************************************************************
  
  else if (cmd_line->current_cmd_operation=="peakdiff") {
    if (strlen(OUT_PREFIX)==0) { fprintf(stderr, "Error: prefix for output files must be specified using the -o option!\n"); exit(1); }
	
    char *SIGNAL_REG_FILES = StrCopy(argv[next_arg]);
    char *REF_REG_FILES = StrCopy(argv[next_arg+1]);
	char *SIGNAL_CONTROL_REG_FILES = (argc>next_arg+2)?StrCopy(argv[next_arg+2]):NULL;
	char *REF_CONTROL_REG_FILES = (argc>next_arg+3)?StrCopy(argv[next_arg+3]):NULL;

	// process input parameters
	int n_signal_files;
	char **signal_reg_file = Tokenize(SIGNAL_REG_FILES,',',&n_signal_files);
	int n_ref_files;
	char **ref_reg_file = Tokenize(REF_REG_FILES,',',&n_ref_files);
    if ((n_signal_files!=2)||(n_ref_files!=2)) { fprintf(stderr, "Error: this method requires exactly two replicates per sample!\n"); exit(1); }
	int n_signal_control_files;
	char **signal_control_reg_file = Tokenize(SIGNAL_CONTROL_REG_FILES,',',&n_signal_control_files);
	int n_ref_control_files;
	char **ref_control_reg_file = Tokenize(REF_CONTROL_REG_FILES,',',&n_ref_control_files);
    if ((n_signal_control_files>0)&&((n_signal_control_files!=n_signal_files)||(n_ref_control_files!=n_ref_files))) { fprintf(stderr, "Error: number of control files should match the number of signal files for each sample!\n"); exit(1); }
	
    // setup output file names
	string data_file_name = (string)OUT_PREFIX + (string)".dat";
	string param_file_name = (string)OUT_PREFIX + (string)".params";
	string rscript_file_name = (string)OUT_PREFIX + (string)".r";
	string image_file_name = (string)OUT_PREFIX + "." + (string)IMAGE_TYPE;
	string log_file_name = (string)OUT_PREFIX + (string)".log";

	// print parameters to file
	FILE *param_file = fopen(param_file_name.c_str(),"w");
    if (param_file==NULL) { fprintf(stderr, "Error: cannot open file '%s' for writing!\n", param_file_name.c_str()); exit(1); }
	fprintf(param_file, "%d\n", n_signal_files);
	fprintf(param_file, "%d\n", n_ref_files);
	fprintf(param_file, "%ld\n", WIN_SIZE);
	fprintf(param_file, "%.6e\n", PVALUE_CUTOFF);
	fprintf(param_file, "%lu\n", PSEUDOCOUNT);
	fprintf(param_file, "%.6e\n", OUTLIER_PROB);
	fprintf(param_file, "%.6e\n", FDR);
	fprintf(param_file, "%s\n", LABELS);
	fprintf(param_file, "%s\n", IMAGE_SIZE);
	fprintf(param_file, "%d\n", IMAGE_RESOLUTION);
	for (int i=0; i<argc; i++) fprintf(param_file, "%s%c", argv[i], i<argc-1?' ':'\n');
    fclose(param_file);

	// scan read files for preliminary peak identification
    if (REUSE==false) ScanReadFiles(signal_reg_file,n_signal_files,ref_reg_file,n_ref_files,signal_control_reg_file,n_signal_control_files,ref_control_reg_file,n_ref_control_files,GENOME_REG_FILE,data_file_name.c_str());

	// execute R script
	CreateRscript(RSCRIPT_INPUT_FILE_NAME,RSCRIPT_TEMPLATE_PEAKDIFF,rscript_file_name);
	fprintf(stderr, "Running R script...\n");
	stringstream s;
	s << "R CMD BATCH '--args " << data_file_name << " " << param_file_name << " " << image_file_name << "' " << rscript_file_name << " " << log_file_name;
	system(s.str().c_str()); 
	PrintLogFile(log_file_name);
	
    // cleanup
    delete [] signal_reg_file;	
	delete [] ref_reg_file;
    if (signal_control_reg_file!=NULL) delete [] signal_control_reg_file;	
    if (ref_control_reg_file!=NULL) delete [] ref_control_reg_file;	
	
    delete SIGNAL_REG_FILES;
    delete REF_REG_FILES;
	if (SIGNAL_CONTROL_REG_FILES) delete SIGNAL_CONTROL_REG_FILES;
	if (REF_CONTROL_REG_FILES) delete REF_CONTROL_REG_FILES;
  }


  // ************************************************************************************
  //
  //     P  R  O  F  I  L  E
  //
  // ************************************************************************************
  
  else if (cmd_line->current_cmd_operation=="profile") {
    if (strlen(OUT_PREFIX)==0) { fprintf(stderr, "Error: prefix for output files must be specified using the -o option!\n"); exit(1); }
	
    char *SIGNAL_REG_FILES = argv[next_arg];
    char *REF_REG_FILES = argv[next_arg+1];
	
	// process input parameters
	char *shift = SHIFT;
	long int shift_upstream = atol(GetNextToken(&shift,','));
	long int shift_downstream = atol(GetNextToken(&shift,','));
	int n_signal_files;
	char **signal_reg_file = Tokenize(SIGNAL_REG_FILES,',',&n_signal_files);
	int n_ref_files;
	char **ref_reg_file = Tokenize(REF_REG_FILES,',',&n_ref_files);
	int n_colors = CountTokens(COLORS,',');
	if (n_colors!=n_signal_files*n_ref_files) { fprintf(stderr, "Error: number of colors must match total number of lines in the plot!\n"); exit(1); }
	int n_legend = CountTokens(LEGEND,',');
	if (n_legend!=n_signal_files*n_ref_files) { fprintf(stderr, "Error: number of legend labels must match total number of lines in the plot!\n"); exit(1); }
	
    // setup output file names
	string data_file_name = (string)OUT_PREFIX + (string)".dat";
	string param_file_name = (string)OUT_PREFIX + (string)".params";
	string rscript_file_name = (string)OUT_PREFIX + (string)".r";
	string image_file_name = (string)OUT_PREFIX + "." + (string)IMAGE_TYPE;
	string log_file_name = (string)OUT_PREFIX + (string)".log";

	// print parameters to file
	FILE *param_file = fopen(param_file_name.c_str(),"w");
    if (param_file==NULL) { fprintf(stderr, "Error: cannot open file '%s' for writing!\n", param_file_name.c_str()); exit(1); }
	fprintf(param_file, "%ld\n", shift_upstream);
	fprintf(param_file, "%ld\n", shift_downstream);
	fprintf(param_file, "%s\n", LEGEND);
	fprintf(param_file, "%s\n", COLORS);
	fprintf(param_file, "%s\n", TITLE);
	fprintf(param_file, "%s\n", XLABEL);
	fprintf(param_file, "%s\n", YLABEL);
	fprintf(param_file, "%s\n", IMAGE_SIZE);
	fprintf(param_file, "%d\n", IMAGE_RESOLUTION);
	for (int i=0; i<argc; i++) fprintf(param_file, "%s%c", argv[i], i<argc-1?' ':'\n');
    fclose(param_file);

	// create R script
	CreateRscript(RSCRIPT_INPUT_FILE_NAME,RSCRIPT_TEMPLATE_PROFILE,rscript_file_name);
	
    if (REUSE==false) {
      // set parameters
	  char *BIN_BITS = StrCopy("17,20,23,26");
	  bool MATCH_GAPS = false;
	  char *OFFSET_OP = StrCopy("5p");
      long int n_bins = (shift_upstream+shift_downstream)/100;
      long int bin_min = -shift_upstream;
      long int bin_max = shift_downstream;

	  // compute profiles!
	  FILE *data_file = fopen(data_file_name.c_str(),"w");
      if (data_file==NULL) { fprintf(stderr, "Error: cannot open file '%s' for writing!\n", data_file_name.c_str()); exit(1); }
	  for (int m=0,k=0; m<n_ref_files; m++) { 
        // open reference region set and shift 5prime position upstream/downstream
        GenomicRegionSet RefRegSet(ref_reg_file[m],BUFFER_SIZE,VERBOSE,true,true);
        for (GenomicRegion *r=RefRegSet.Get(); r!=NULL; r=RefRegSet.Next()) { 
	      r->Connect();
	      r->Center();
          r->I.front()->ShiftPos(-shift_upstream,shift_downstream,true);
	    }
	  
	    // process signal files
	    for (int n=0; n<n_signal_files; n++,k++) {
	      if (VERBOSE) fprintf(stderr, "Creating profile of '%s' in '%s'...\n", signal_reg_file[n], ref_reg_file[m]); 
          GenomicRegionSet TestRegSet(signal_reg_file[n],BUFFER_SIZE,VERBOSE,false,true);
          UnsortedGenomicRegionSetOverlaps Overlaps(&TestRegSet,&RefRegSet,BIN_BITS);
          Progress PRG("Processing queries...",1);
          unsigned long int *bins;
          ALLOCATE1D_INIT(bins,n_bins,unsigned long int,0);
		  unsigned long int n_signal_reg = 0;
          for (GenomicRegion *qreg=Overlaps.GetQuery(); Overlaps.Done()==false; qreg=Overlaps.NextQuery()) {
		    n_signal_reg++;
            for (GenomicRegion *ireg=Overlaps.GetOverlap(MATCH_GAPS,IGNORE_STRAND); ireg!=NULL; ireg=Overlaps.NextOverlap(MATCH_GAPS,IGNORE_STRAND)) {
              long int start_offset, stop_offset;
              qreg->I.front()->GetOffsetFrom(ireg->I.front(),OFFSET_OP,IGNORE_STRAND,&start_offset,&stop_offset);
              if (start_offset>stop_offset) { fprintf(stderr, "Error: start offset is greater than stop offest (this must be a bug)!\n"); exit(1); }
              long int x = (start_offset+stop_offset)/2-shift_upstream;
              double z = (double)(x-bin_min)/(bin_max-bin_min);
              if ((z>=0)&&(z<1)) bins[(int)(n_bins*z)]++;
            }
            PRG.Check();
          }
          PRG.Done();
	      fprintf(data_file, "%s in %s\t", signal_reg_file[n], ref_reg_file[m]);
          if (NORMALIZE) for (long int b=0; b<n_bins; b++) fprintf(data_file, "%.6e%c", (double)bins[b]/n_signal_reg/RefRegSet.n_regions, b!=n_bins-1?'\t':'\n');
          else for (long int b=0; b<n_bins; b++) fprintf(data_file, "%lu%c", bins[b], b!=n_bins-1?'\t':'\n');
          delete bins;
	    }
		
	  }
	  fclose(data_file);
      delete BIN_BITS;
	  delete OFFSET_OP;
	}
	
	// execute R script
	fprintf(stderr, "Running R script...\n");
	stringstream s;
	s << "R CMD BATCH '--args " << data_file_name << " " << param_file_name << " " << image_file_name << "' " << rscript_file_name << " " << log_file_name;
	system(s.str().c_str()); 
	PrintLogFile(log_file_name);

	
    // cleanup
	delete [] signal_reg_file;
	delete [] ref_reg_file;
	
  }
  
  
  // clean up
  delete cmd_line;
  delete RANDOM_GENERATOR;
  
  return 0;
}



