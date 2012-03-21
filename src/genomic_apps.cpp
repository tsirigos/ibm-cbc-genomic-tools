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



const char *RSCRIPT_TEMPLATE_PROFILE = \
"\n##\n## USAGE: genomic_apps.profile.r DATA-FILE PARAMETER-FILE OUTPUT-IMAGE-FILE\n##\n\nargs <- commandArgs(trailingOnly=T);\nprofile_data_file <- args[1];\nprofile_param_file <- args[2];\nprofile_image_file <- args[3];\n\nparams <- readLines(profile_param_file);\nprofile_upstream <- as.numeric(params[1]);\nprofile_downstream <- as.numeric(params[2]);\nprofile_legend <- strsplit(params[3],',')[[1]];\nprofile_colors <- strsplit(params[4],',')[[1]];\nprofile_title <- params[5];\nprofile_xlab <- params[6];\nprofile_ylab <- params[7];\n  \ntiff(profile_image_file,width=3500,height=3500,res=600,compression='lzw',antialias='none');\nY <- as.matrix(read.table(profile_data_file,header=FALSE,row.names=1,sep='\t'));\n\nd <- (profile_upstream+profile_downstream)/ncol(Y);\nx <- seq(-profile_upstream+d/2,+profile_downstream-d/2,d);\nplot(x,Y[1,],type='l',col=profile_colors[1],main=profile_title,xlab=profile_xlab,ylab=profile_ylab,xlim=c(-profile_upstream,profile_upstream),ylim=c(0,max(Y)));\ni <- 2;\nwhile (i <= nrow(Y)) \n{\n  lines(x,Y[i,],col=profile_colors[i]);\n  i <- i + 1;\n}\nlegend('topright',profile_legend,pch=16,col=profile_colors,inset=0.05);\ndev.off();\n\n";


const char *RSCRIPT_TEMPLATE_HEATMAP = \
"\n##\n## USAGE: genomic_apps.heatmap.r DATA-FILE PARAMETER-FILE OUTPUT-IMAGE-FILE\n##\n\nargs <- commandArgs(trailingOnly=T);\ndata_file <- args[1];\nparam_file <- args[2];\nimage_file <- args[3];\n\nparams <- readLines(param_file);\nshift_upstream <- as.numeric(params[1]);\nshift_downstream <- as.numeric(params[2]);\nheatmap_colors <- strsplit(params[3],',')[[1]];\nheatmap_title <- strsplit(params[4],',')[[1]];\nheatmap_xlab <- params[5];\nheatmap_ylab <- params[6];\nheatmap_size <- as.numeric(strsplit(params[7],',')[[1]]);\nheatmap_resolution <- as.numeric(params[8]);\nn_heatmaps <- as.numeric(params[9]);\n\n\nnorm_rows <- function(X) { for (i in 1:nrow(X)) X[i,] <- (X[i,]-mean(X[i,]))/sd(X[i,]); return(X); }\n\nlibrary('MASS');\nlibrary('preprocessCore');           # from Bioconductor\nlibrary('gplots');                   # from Bioconductor\n\n# load data\nD <- as.matrix(read.table(data_file,row.names=1,sep='\t'));\n  \n# create combined heatmap (main version)\ntiff(image_file,width=heatmap_size[1],height=heatmap_size[2],res=heatmap_resolution,compression='lzw');\n\npar(fig=c(0,1,0,1),mar=c(2,2,0,0)); \nplot.new();\nmtext(heatmap_xlab,side=1);\nmtext(heatmap_ylab,side=2);\n\n\nd <- ncol(D)/n_heatmaps;\nI <- 1:d;\ndj <- 0.9/n_heatmaps;\nfor (j in 1:n_heatmaps) {\n  par(fig=c(0.1+(j-1)*dj,0.1+j*dj,0.1,1),mar=c(0.5,0.5,3,0.5),new=TRUE);\n  colorscale <- c(colorpanel(20,low='white',high='white'),colorpanel(50,low='white',high=heatmap_colors[j]),colorpanel(30,low=heatmap_colors[j],high=heatmap_colors[j]));\n  image(z=t(norm_rows(D[,I])),col=colorscale,main=heatmap_title[j],xlab=heatmap_xlab,ylab=heatmap_ylab,xaxt='n',yaxt='n');\n  I <- I+d;\n}\n\ndev.off();\n\n";





  
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
char *IMAGE_SIZE;
int IMAGE_RESOLUTION;
bool NORMALIZE;





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

/*
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
    cmd_line->AddOption("-isize", &IMAGE_SIZE, "2000,4000", "comma-separated image dimensions");
    cmd_line->AddOption("-ires", &IMAGE_RESOLUTION, 600, "image resolution in dpi");
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
    cmd_line->AddOption("-norm", &NORMALIZE, false, "normalize against total number of reads and number of reference regions");	
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



//-------CreateRscript-----------
//
void CreateRscript(const char *rscript_input_file_name, const char *rscript_template, string rscript_output_file_name)
{
  if (strlen(rscript_input_file_name)==0) {
    FILE *rscript_out_file = fopen(rscript_output_file_name.c_str(),"w");
	fprintf(rscript_out_file, "%s", rscript_template);
    fclose(rscript_out_file);
  }
  else {
    FILE *rscript_in_file = fopen(rscript_input_file_name,"r");
	if (rscript_in_file==NULL) { fprintf(stderr, "Error: R script file '%s' not found!\n", rscript_input_file_name); exit(1); }
	char *rscript = LoadFile(rscript_in_file);
	fclose(rscript_in_file);
    FILE *rscript_out_file = fopen(rscript_output_file_name.c_str(),"w");
	fprintf(rscript_out_file, "%s", rscript);
    fclose(rscript_out_file);
  }
}




//-------CreateRscript-----------
//
void PrintLogFile(string log_file_name)
{
  fprintf(stderr, "\n***********************************************\n"); 
  fprintf(stderr, "\n  R   S C R I P T   L O G   F I L E            \n"); 
  fprintf(stderr, "\n***********************************************\n"); 
  FileBuffer buffer(log_file_name.c_str());
  for (char *inp=buffer.Next(); inp!=NULL; inp=buffer.Next()) fprintf(stderr, "%s\n", inp);
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
	string image_file_name = (string)OUT_PREFIX + (string)".tif";
	string log_file_name = (string)OUT_PREFIX + (string)".log";

	// print parameters to file
	FILE *param_file = fopen(param_file_name.c_str(),"w");
	fprintf(param_file, "%ld\n", shift_upstream);
	fprintf(param_file, "%ld\n", shift_downstream);
	fprintf(param_file, "%s\n", COLORS);
	fprintf(param_file, "%s\n", TITLE);
	fprintf(param_file, "%s\n", XLABEL);
	fprintf(param_file, "%s\n", YLABEL);
	fprintf(param_file, "%s\n", IMAGE_SIZE);
	fprintf(param_file, "%d\n", IMAGE_RESOLUTION);
	fprintf(param_file, "%d\n", n_signal_files);
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



