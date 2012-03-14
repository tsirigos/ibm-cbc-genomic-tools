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


using namespace std;


//---------------------------------------------------------------------------------//
// Global variables                                                                //
//---------------------------------------------------------------------------------//

const string PROGRAM = "genomic_apps";
const long int BUFFER_SIZE = 10000;
gsl_rng *RANDOM_GENERATOR;



const char *RSCRIPT_PROFILE = "\n##\n## USAGE: genomic_apps.profile.r DATA-FILE PARAMETER-FILE OUTPUT-IMAGE-FILE\n##\n\nargs <- commandArgs(trailingOnly=T);\nprofile_data_file <- args[1];\nprofile_param_file <- args[2];\nprofile_image_file <- args[3];\n\nparams <- readLines(profile_param_file);\nprofile_upstream <- as.numeric(params[1]);\nprofile_downstream <- as.numeric(params[2]);\nprofile_legend <- strsplit(params[3],',')[[1]];\nprofile_colors <- strsplit(params[4],',')[[1]];\nprofile_title <- params[5];\nprofile_xlab <- params[6];\nprofile_ylab <- params[7];\n  \ntiff(profile_image_file,width=3500,height=3500,res=600,compression='lzw',antialias='none');\nY <- as.matrix(read.table(profile_data_file,header=FALSE,row.names=1,sep='\t'));\n\nd <- (profile_upstream+profile_downstream)/ncol(Y);\nx <- seq(-profile_upstream+d/2,+profile_downstream-d/2,d);\nplot(x,Y[1,],type='l',col=profile_colors[1],main=profile_title,xlab=profile_xlab,ylab=profile_ylab,xlim=c(-profile_upstream,profile_upstream),ylim=c(0,max(Y)));\ni <- 2;\nwhile (i <= nrow(Y)) \n{\n  lines(x,Y[i,],col=profile_colors[i]);\n  i <- i + 1;\n}\nlegend('topright',profile_legend,pch=16,col=profile_colors,inset=0.05);\ndev.off();\n\n";




  
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



//-------InitCmdLine-----------
//
CmdLineWithOperations *InitCmdLine(int argc, char *argv[], int *next_arg)
{
  // initialize
  CmdLineWithOperations *cmd_line = new CmdLineWithOperations(); 
  cmd_line->SetProgramName(PROGRAM,VERSION);

  // set operations
/*  cmd_line->AddOperation("heatmap", "[OPTIONS] <REG-FILE>", \
  "Determines input read counts in sliding windows of reference regions.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operand: interval\n\
  * Region requirements: single-interval\n\
  * Region-set requirements: sorted by chromosome/strand/start"\
  );

  cmd_line->AddOperation("peakdiff", "[OPTIONS] SIGNAL-REG-FILE [CONTROL-REG-FILE [GENOME-UNIQ-REG-FILE]]", \
  "Scans input reads to identify peaks.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operand: interval\n\
  * Region requirements: single-interval\n\
  * Region-set requirements: sorted by chromosome/strand/start"\
  );
*/

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
    //n_args = 0;
    //cmd_line->AddOption("-g", &GENOME_REG_FILE, "genome.reg+", "genome region file");
  }
  else if (cmd_line->current_cmd_operation=="peakdiff") {
    //n_args = 0;
    //cmd_line->AddOption("-g", &GENOME_REG_FILE, "genome.reg+", "genome region file");
  }
  else if (cmd_line->current_cmd_operation=="profile") {
    n_args = 2;
    cmd_line->AddOption("-reuse", &REUSE, false, "reuse histogram data; update paramaters only");	
    cmd_line->AddOption("-R", &RSCRIPT_INPUT_FILE_NAME, "", "R script file to use (not required)");	
    cmd_line->AddOption("-o", &OUT_PREFIX, "", "prefix for output files");
    cmd_line->AddOption("-i", &IGNORE_STRAND, false, "ignore strand while finding overlaps");
    cmd_line->AddOption("-shift", &SHIFT, "5000,5000", "comma-separated upstream/downstream distances from reference center");
    cmd_line->AddOption("-legend", &LEGEND, "", "comma-separated legend labels for line plot");
    cmd_line->AddOption("-colors", &COLORS, "", "comma-separated colors for line plot");
    cmd_line->AddOption("-title", &TITLE, "", "plot title");
    cmd_line->AddOption("-xlab", &XLABEL, "", "plot x-axis label");
    cmd_line->AddOption("-ylab", &YLABEL, "", "plot y-axis label");
  }
  else {
    cerr << "Unknown operation '" << op << "'!\n";
    delete cmd_line;
    exit(1);
  }

  // process command line
  *next_arg = cmd_line->Read(argv+1,argc-1) + 1;
  if (HELP||(argc-*next_arg<n_args)) { cmd_line->OperationUsage(); delete cmd_line; exit(1); }
  
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




//---------------------------------------------------------------------------------//
// MAIN	                                                                           //
//---------------------------------------------------------------------------------//
int main(int argc, char* argv[]) 
{
  // process command-line arguments
  int next_arg;
  CmdLineWithOperations *cmd_line = InitCmdLine(argc,argv,&next_arg); 
  _MESSAGES_ = VERBOSE;
  
  // initilize random generator
  RANDOM_GENERATOR = InitRandomGenerator(getpid()+time(NULL));

  
  if (cmd_line->current_cmd_operation=="profile") {
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
	
    // setup output file names
	string data_file_name = (string)OUT_PREFIX + (string)".dat";
	string param_file_name = (string)OUT_PREFIX + (string)".params";
	string rscript_file_name = (string)OUT_PREFIX + (string)".r";
	string image_file_name = (string)OUT_PREFIX + (string)".tif";
	string log_file_name = (string)OUT_PREFIX + (string)".log";

	// print parameters to file
	FILE *param_file = fopen(param_file_name.c_str(),"w");
	fprintf(param_file, "%ld\n", shift_upstream);
	fprintf(param_file, "%ld\n", shift_downstream);
	fprintf(param_file, "%s\n", LEGEND);
	fprintf(param_file, "%s\n", COLORS);
	fprintf(param_file, "%s\n", TITLE);
	fprintf(param_file, "%s\n", XLABEL);
	fprintf(param_file, "%s\n", YLABEL);
    fclose(param_file);

	// create R script
	if (strlen(RSCRIPT_INPUT_FILE_NAME)==0) {
      FILE *rscript_out_file = fopen(rscript_file_name.c_str(),"w");
	  fprintf(rscript_out_file, "%s", RSCRIPT_PROFILE);
      fclose(rscript_out_file);
	}
	else {
	  FILE *rscript_in_file = fopen(RSCRIPT_INPUT_FILE_NAME,"r");
	  if (rscript_in_file==NULL) { fprintf(stderr, "Error: R script file '%s' not found!\n", RSCRIPT_INPUT_FILE_NAME); exit(1); }
	  char *rscript = LoadFile(rscript_in_file);
	  fclose(rscript_in_file);
      FILE *rscript_out_file = fopen(rscript_file_name.c_str(),"w");
	  fprintf(rscript_out_file, "%s", rscript);
      fclose(rscript_out_file);
	}

    if (REUSE==false) {
      // set parameters
	  char *BIN_BITS = StrCopy("17,20,23,26");
	  bool MATCH_GAPS = false;
	  char *OFFSET_OP = StrCopy("5p");
      long int n_bins = 100;
      long int bin_min = -shift_upstream;
      long int bin_max = shift_downstream;

	  // compute profiles!
	  FILE *data_file = fopen(data_file_name.c_str(),"w");
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
          double *bins;
          ALLOCATE1D_INIT(bins,n_bins,double,0);
          for (GenomicRegion *qreg=Overlaps.GetQuery(); Overlaps.Done()==false; qreg=Overlaps.NextQuery()) {
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
          for (long int b=0; b<n_bins; b++) fprintf(data_file, "%.0f%c", bins[b], b!=n_bins-1?'\t':'\n');
          delete bins;
	    }
		
	  }
	  fclose(data_file);
      delete BIN_BITS;
	  delete OFFSET_OP;
	}
	
	// execute R script
	stringstream s;
	s << "R CMD BATCH '--args " << data_file_name << " " << param_file_name << " " << image_file_name << "' " << rscript_file_name << " " << log_file_name;
	system(s.str().c_str()); 

	
    // cleanup
	delete [] signal_reg_file;
	delete [] ref_reg_file;
	
  }
  
  
  // clean up
  delete cmd_line;
  delete RANDOM_GENERATOR;
  
  return 0;
}



