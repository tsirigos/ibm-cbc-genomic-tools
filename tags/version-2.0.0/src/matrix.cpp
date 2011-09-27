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
#include <sys/types.h>
#include <sys/stat.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <string>
#include <iostream>
#include "core.h"
using namespace std;



//---------------------------------------------------------------------------------//
// COMMAND-LINE OPTIONS                                                            //
//---------------------------------------------------------------------------------//


// Generic options
bool VERBOSE;
int DECIMALS;
char *FMT_STR = NULL;
bool FMT_INT;
bool FMT_LINT;
bool FMT_LLINT;
bool RANGE;
bool UPPER;
bool EQUAL;
bool GREATER;
bool HELP;
int OFFSET;
int NBINS;
char *AVG_STR = NULL;
char *STD_STR = NULL;
char *COEFF_STR = NULL;
bool FIRST_COLUMN;
bool ALL_PAIRS;


gsl_rng *RANDOM_GENERATOR;



//-----Usage----------
//
void Usage()
{
  fprintf(stderr, "\n");
  fprintf(stderr, "USAGE: \n");
  fprintf(stderr, "  matrix OPERATION [OPTIONS] <MATRIX>|MATRIX_FILE\n"); 
  fprintf(stderr, "\n");
  fprintf(stderr, "OPERATIONS: \n");
  fprintf(stderr, "  * -T          matrix transpose\n");
  fprintf(stderr, "  * -norm       normalization\n");
  fprintf(stderr, "  * -cnorm      column normalization\n");
  fprintf(stderr, "  * -rnorm      row normalization\n");
  fprintf(stderr, "  * -cstat      column statistics\n");
  fprintf(stderr, "  * -rstat      row statistics\n");
  fprintf(stderr, "  * -csum       column sums\n");
  fprintf(stderr, "  * -cutoff     apply cutoff\n");
  fprintf(stderr, "  * -test       test equality and/or inequality\n");
  fprintf(stderr, "  * -pairs      print all row pairs\n");
  fprintf(stderr, "  * -sparse     convert matrix to sparse format\n");
  fprintf(stderr, "  * -encode     print vectors as encoded sets\n");
  fprintf(stderr, "  * -shuffle    permute the elements of each row separately\n");
  fprintf(stderr, "  * -format     format matrix\n");
  fprintf(stderr, "  * -shrink     remove zero/empty columns\n");
  fprintf(stderr, "  * -mult       multiply matrix elements with coefficient\n");
  fprintf(stderr, "  * -rel        compute relative changes between element pairs in each row\n");
  fprintf(stderr, "  * -del        delete redundant rows\n");
  fprintf(stderr, "\n");
}



//-----Usage----------
//
void Usage(CmdLine *cmdLine, char *method)
{
  if (strcmp(method,"-T")==0) cmdLine->Usage("matrix -T [OPTIONS] MATRIX_FILE => MATRIX"); 
  else if (strcmp(method,"-norm")==0) cmdLine->Usage("matrix -norm [OPTIONS] MATRIX => MATRIX"); 
  else if (strcmp(method,"-cnorm")==0) cmdLine->Usage("matrix -cnorm [OPTIONS] MATRIX => MATRIX"); 
  else if (strcmp(method,"-rnorm")==0) cmdLine->Usage("matrix -rnorm [OPTIONS] MATRIX => MATRIX"); 
  else if (strcmp(method,"-cstat")==0) cmdLine->Usage("matrix -cstat [OPTIONS] MATRIX => STATISTICS"); 
  else if (strcmp(method,"-rstat")==0) cmdLine->Usage("matrix -rstat [OPTIONS] MATRIX => STATISTICS"); 
  else if (strcmp(method,"-csum")==0) cmdLine->Usage("matrix -csum [OPTIONS] MATRIX => VECTOR"); 
  else if (strcmp(method,"-cutoff")==0) cmdLine->Usage("matrix -cutoff [OPTIONS] MATRIX => MATRIX"); 
  else if (strcmp(method,"-test")==0) cmdLine->Usage("matrix -test [OPTIONS] MATRIX => MATRIX"); 
  else if (strcmp(method,"-pairs")==0) cmdLine->Usage("matrix -pairs [OPTIONS] MATRIX => VECTOR_PAIRS"); 
  else if (strcmp(method,"-sparse")==0) cmdLine->Usage("matrix -sparse [OPTIONS] MATRIX => SPARSE_VECTORS"); 
  else if (strcmp(method,"-encode")==0) cmdLine->Usage("matrix -encode [OPTIONS] MATRIX => SETS"); 
  else if (strcmp(method,"-shuffle")==0) cmdLine->Usage("matrix -shuffle [OPTIONS] MATRIX => MATRIX"); 
  else if (strcmp(method,"-format")==0) cmdLine->Usage("matrix -format [OPTIONS] MATRIX => MATRIX"); 
  else if (strcmp(method,"-shrink")==0) cmdLine->Usage("matrix -shrink [OPTIONS] MATRIX => MATRIX"); 
  else if (strcmp(method,"-mult")==0) cmdLine->Usage("matrix -mult [OPTIONS] MATRIX => MATRIX"); 
  else if (strcmp(method,"-rel")==0) cmdLine->Usage("matrix -rel [OPTIONS] MATRIX => MATRIX"); 
  else if (strcmp(method,"-del")==0) cmdLine->Usage("matrix -del [OPTIONS] MATRIX => MATRIX"); 
  else fprintf(stderr, "Unknown method '%s'!\n", method);
  exit(1); 
}



//--------InitCmdLine-----------
//
CmdLine *InitCmdLine(char *method)
{
  CmdLine *cmd_line = new CmdLine(); 

  // Processing
  cmd_line->AddOption("-v", &VERBOSE, false, "verbose mode");
  cmd_line->AddOption("-h", &HELP, false, "help");
  cmd_line->AddOption("-help", &HELP, false, "help");
  cmd_line->AddOption("-f", &FMT_STR, "", "format string");
  cmd_line->AddOption("-i", &FMT_INT, false, "integer format");
  cmd_line->AddOption("-l", &FMT_LINT, false, "long integer format");
  cmd_line->AddOption("-ll", &FMT_LLINT, false, "long long integer format");
  cmd_line->AddOption("-n", &DECIMALS, 2, "double format decimal numbers");

  if (strcmp(method,"-cutoff")==0) {
    cmd_line->AddOption("-c", &COEFF_STR, "0", "cutoff value");
    cmd_line->AddOption("-u", &UPPER, false, "use upper bound for cutoff");
  }
  else if (strcmp(method,"-test")==0) {
    cmd_line->AddOption("-c", &COEFF_STR, "0", "cutoff value");
    cmd_line->AddOption("-e", &EQUAL, false, "test equality");
    cmd_line->AddOption("-g", &GREATER, false, "test greater than");
  }
  else if (strcmp(method,"-norm")==0) {
    cmd_line->AddOption("-m", &AVG_STR, "0", "mean");
    cmd_line->AddOption("-s", &STD_STR, "1", "standard deviation");
  }
  else if (strcmp(method,"-rnorm")==0) {
    cmd_line->AddOption("-r", &RANGE, false, "normalize vector range (default is mean/std)");
  }
  else if (strcmp(method,"-sparse")==0) {
    cmd_line->AddOption("-d", &OFFSET, false, "starting feature number");
  }
  else if (strcmp(method,"-encode")==0) {
    cmd_line->AddOption("-b", &NBINS, 2, "number of bins");
  }
  else if (strcmp(method,"-mult")==0) {
    cmd_line->AddOption("-c", &COEFF_STR, "0", "cutoff value");
  }
  else if (strcmp(method,"-rel")==0) {
    cmd_line->AddOption("-1", &FIRST_COLUMN, false, "use first column as reference");
    cmd_line->AddOption("-2", &ALL_PAIRS, false, "compute all pairs");
  }
  else if (strcmp(method,"-del")==0) {
    cmd_line->AddOption("-c", &COEFF_STR, "0", "cutoff value");
    cmd_line->AddOption("-e", &EQUAL, false, "test equality");
    cmd_line->AddOption("-g", &GREATER, false, "test greater than");
  }
  
  return cmd_line;
}








//-------RunMatrixCommand-------
//
template <class T> void RunMatrixCommand(char *file, char *method, char *fmt, T coeff, T avg, T std)
{
  MatrixTemplate<T> M(file,VERBOSE);
  if (strcmp(method,"-norm")==0) { M.Norm(avg,std); M.Print(fmt);  }			// OPTION -norm: normalize matrix
  else if (strcmp(method,"-rnorm")==0) { M.NormRows(RANGE); M.Print(fmt);  }		// OPTION -rnorm: normalize vectors
  else if (strcmp(method,"-cnorm")==0) { M.NormCols(); M.Print(fmt); }			// OPTION -cnorm: normalize columns
  else if (strcmp(method,"-cutoff")==0) { M.ApplyCutoff(coeff,UPPER); M.Print(fmt); }	// OPTION -cutoff: apply cutoff filter
  else if (strcmp(method,"-test")==0) { M.Test(coeff,EQUAL,GREATER); M.Print(fmt); }	// OPTION -test: test equality and inequality
  else if (strcmp(method,"-cstat")==0) M.PrintColStats(fmt);				// OPTION -cstat: matrix column statistics
  else if (strcmp(method,"-rstat")==0) M.PrintRowStats(fmt);				// OPTION -rstat: matrix row statistics
  else if (strcmp(method,"-csum")==0) M.PrintColSums(fmt);				// OPTION -csum: matrix column sums
  else if (strcmp(method,"-T")==0) M.PrintTranspose(fmt);				// OPTION -T: print matrix transpose
  else if (strcmp(method,"-pairs")==0) M.PrintRowPairs(fmt); 				// OPTION -pairs: print all row pairs
  else if (strcmp(method,"-sparse")==0) M.PrintSparse(fmt,OFFSET);                      // OPTION -sparse: convert matrix to sparse format
  else if (strcmp(method,"-encode")==0) M.PrintEncoded(NBINS);   			// OPTION -encoded: print encoded sets
  else if (strcmp(method,"-shuffle")==0) { M.Shuffle(RANDOM_GENERATOR); M.Print(fmt); }	// OPTION -shuffle: permute the elements of each row separately
  else if (strcmp(method,"-format")==0) M.Print(fmt);					// OPTION -format: format matrix
  else if (strcmp(method,"-shrink")==0) M.Shrink(fmt);					// OPTION -shrink: remove zero/empty columns
  else if (strcmp(method,"-mult")==0) { M.Multiply(coeff); M.Print(fmt); }		// OPTION -mult: multiply matrix elements with coefficient
  else if (strcmp(method,"-rel")==0) { 							// OPTION -rel: compute relative changes between element pairs in each row
    if (ALL_PAIRS==false) { M.Rel(FIRST_COLUMN); M.Print(fmt); }
    else M.PrintRel2(fmt);
  }
  else if (strcmp(method,"-del")==0) M.Delete(fmt,coeff,EQUAL,GREATER);			// OPTION -del: delete redundant rows
  else { fprintf(stderr, "Error: '%s' is an invalid method!\n", method); exit(1); }	// INVALID OPTION
}







//---------------------------------------------------------------------------------//
// MAIN	                                                                           //
//---------------------------------------------------------------------------------//
int main(int argc, char* argv[]) 
{
  // process command-line arguments
  if (argc<2) { Usage(); exit(1); }
  char *method = argv[1];
  CmdLine *cmdLine = InitCmdLine(method);
  int next_arg = cmdLine->Read(argv+1,argc-1) + 1;
  if (HELP) { Usage(cmdLine,method); exit(1); }
  MESSAGES(VERBOSE);
  if (COEFF_STR==NULL) COEFF_STR = StrCopy("0");
  if (AVG_STR==NULL) AVG_STR = StrCopy("0");
  if (STD_STR==NULL) STD_STR = StrCopy("0");

  // initialize random generator
  RANDOM_GENERATOR = InitRandomGenerator(time(NULL));

  // setup output format
  char *MATRIX_FILE = next_arg>=argc ? NULL : argv[next_arg];
  char fmt[100];
  if (FMT_INT) { strcpy(fmt,"%d"); RunMatrixCommand(MATRIX_FILE,method,fmt,atoi(COEFF_STR),atoi(AVG_STR),atoi(STD_STR)); }
  else if (FMT_LINT) { strcpy(fmt,"%ld"); RunMatrixCommand(MATRIX_FILE,method,fmt,atol(COEFF_STR),atol(AVG_STR),atol(STD_STR)); }
  else if (FMT_LLINT) { strcpy(fmt,"%lld"); RunMatrixCommand(MATRIX_FILE,method,fmt,atoll(COEFF_STR),atoll(AVG_STR),atoll(STD_STR)); }
  else { 
    if (strlen(FMT_STR)==0) sprintf(fmt, "%%.%if", DECIMALS); 
    else sprintf(fmt, "%s", FMT_STR);
    RunMatrixCommand(MATRIX_FILE,method,fmt,(double)atof(COEFF_STR),(double)atof(AVG_STR),(double)atof(STD_STR)); 
  }

  // clean up
  delete cmdLine;
  gsl_rng_free(RANDOM_GENERATOR);

  return 0;
}



