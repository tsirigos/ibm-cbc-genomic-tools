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
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <iostream>
#include "core.h"
using namespace std;



extern "C" pid_t getpid(); 




//---------------------------------------------------------------------------------//
// COMMAND-LINE OPTIONS                                                            //
//---------------------------------------------------------------------------------//


// Generic options
bool VERBOSE;
bool HELP;
char *FMT_STR;
int DECIMALS;
long int VEC_SIZE;
char *DELIMS;
bool GREATER;
bool EQUAL;
float CUTOFF;
float CUTOFF_REPLACE;
bool UPPER;
long int OFFSET;
int NBINS;
bool RANGE;
long int RADIUS;
long int DISTANCE;
bool INFO;
long int BUFFER_SIZE;
float ALPHA;
bool ZEROES;
float LOG_BASE;
float EXP_BASE;
float POW;
int HIST_BIN_SIZE;
int HIST_NBINS;
int HIST_NBINS_COMBINE;
double HIST_MIN;
double HIST_MAX;
bool HIST_AUTOMIN;
bool HIST_AUTOMAX;
char *SEPARATOR;
bool RANK_NORM;
double CONST;






//-----Usage----------
//
void Usage()
{
  fprintf(stderr, "\n");
  fprintf(stderr, "USAGE: \n");
  fprintf(stderr, "  vectors OPERATION [OPTIONS] <VECTORS>\n"); 
  fprintf(stderr, "\n");
  fprintf(stderr, "FUNCTION: \n");
  fprintf(stderr, "  Perform vector operations using double-precision floating-point arithmetic.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "OPERATIONS: \n");
  fprintf(stderr, "  absmax     absolute maximum\n");
  fprintf(stderr, "  bins       create bins for vector values\n");
  fprintf(stderr, "  cutoff     apply cutoff\n");
  fprintf(stderr, "  div        divide\n");
  fprintf(stderr, "  exp        compute the exponential (inverse of log)\n");
  fprintf(stderr, "  fold       compute fold changes\n");
  fprintf(stderr, "  format     format vectors\n");
  fprintf(stderr, "  hist       histogram (all vectors)\n");
  fprintf(stderr, "  imax       maximum index\n");
  fprintf(stderr, "  imin       minimum index\n");
  fprintf(stderr, "  items      convert to itemset format\n");
  fprintf(stderr, "  log        compute the logarithm\n");
  fprintf(stderr, "  m          mean\n");
  fprintf(stderr, "  max        maximum\n");
  fprintf(stderr, "  med        median\n");
  fprintf(stderr, "  merge      merge consecutive vectors with identical labels\n");
  fprintf(stderr, "  min        minimum\n");
  fprintf(stderr, "  n          size\n");
  fprintf(stderr, "  norm       normalize\n");
  fprintf(stderr, "  pairs      create all vector pairs\n");
  fprintf(stderr, "  permute    permute order of vectors\n");
  fprintf(stderr, "  pow        compute the power\n");
  fprintf(stderr, "  q          quantiles\n");
  fprintf(stderr, "  rank       ranks\n");
  fprintf(stderr, "  rev        reverse\n");
  fprintf(stderr, "  sd         standard deviation\n");
  fprintf(stderr, "  shuffle    shuffle vector elements\n");
  fprintf(stderr, "  slide      create subvectors by a sliding window\n");
  fprintf(stderr, "  sort       sort\n");
  fprintf(stderr, "  sparse     convert to sparse format\n");
  fprintf(stderr, "  stat       statistics\n");
  fprintf(stderr, "  sum        sum\n");
  fprintf(stderr, "  test       test if greater than cutoff\n");
  fprintf(stderr, "\n");
}



//--------InitCmdLine-----------
//
CmdLine *InitCmdLine(char *method)
{
  CmdLine *cmd_line = new CmdLine(); 

  // Processing
  cmd_line->AddOption("-v", &VERBOSE, false, "verbose mode");
  cmd_line->AddOption("--help", &HELP, false, "help");
  cmd_line->AddOption("-h", &HELP, false, "help");
  cmd_line->AddOption("-f", &FMT_STR, "", "format string");
  cmd_line->AddOption("-n", &DECIMALS, 2, "number of decimal points");
  cmd_line->AddOption("-B", &BUFFER_SIZE, 100000000, "input line buffer size");
  cmd_line->AddOption("-k", &VEC_SIZE, 0, "vector size (0 = auto)");

  if (strcmp(method,"cutoff")==0) {
    cmd_line->AddOption("-u", &UPPER, false, "use upper bound for cutoff");
    cmd_line->AddOption("-c", &CUTOFF, 0, "cutoff value");
    cmd_line->AddOption("-r", &CUTOFF_REPLACE, 0, "cutoff replacement value");
  }
  else if (strcmp(method,"exp")==0) {
    cmd_line->AddOption("-b", &EXP_BASE, 10.0, "exp base");
  }
  else if (strcmp(method,"format")==0) {
    cmd_line->AddOption("-d", &DELIMS, "", "vector delimiters");
  }
  else if (strcmp(method,"hist")==0) {
    cmd_line->AddOption("-b", &HIST_NBINS, 10, "number of bins");
    cmd_line->AddOption("-min", &HIST_MIN, 0.0, "minimum");
    cmd_line->AddOption("-max", &HIST_MAX, 1.0, "maximum");
  }
  else if (strcmp(method,"bins")==0) {
    cmd_line->AddOption("--bin-size", &HIST_BIN_SIZE, 0, "bin size (if > 0, overrides -b)");
    cmd_line->AddOption("-b", &HIST_NBINS, 10, "number of bins");
    cmd_line->AddOption("-m", &HIST_NBINS_COMBINE, 1, "number of sucessive bins to combine");
    cmd_line->AddOption("-min", &HIST_MIN, 0.0, "minimum");
    cmd_line->AddOption("-max", &HIST_MAX, 1.0, "maximum");
    cmd_line->AddOption("--auto-min", &HIST_AUTOMIN, false, "determine minimum automatically");
    cmd_line->AddOption("--auto-max", &HIST_AUTOMAX, false, "determine maximum automatically");
  }
  else if (strcmp(method,"items")==0) {
    cmd_line->AddOption("-b", &NBINS, 3, "number of bins");
  }
  else if (strcmp(method,"log")==0) {
    cmd_line->AddOption("-b", &LOG_BASE, 10.0, "log base");
    cmd_line->AddOption("-c", &CONST, 0.0, "add constant before computing logarithm");
  }
  else if (strcmp(method,"m")==0) {
    cmd_line->AddOption("-a", &ALPHA, 0.0, "trimmed vector parameter");
  }
  else if (strcmp(method,"merge")==0) {
    cmd_line->AddOption("-t", &SEPARATOR, " ", "separator");
  }
  else if (strcmp(method,"norm")==0) {
    cmd_line->AddOption("-r", &RANGE, false, "normalize vector range (default is mean/std)");
    cmd_line->AddOption("-a", &ALPHA, 0.0, "trimmed vector parameter");
  }
  else if (strcmp(method,"pow")==0) {
    cmd_line->AddOption("-b", &POW, 10.0, "power of");
  }
  else if (strcmp(method,"rank")==0) {
    cmd_line->AddOption("-norm", &RANK_NORM, false, "normalized ranks");
  }
  else if (strcmp(method,"sd")==0) {
    cmd_line->AddOption("-a", &ALPHA, 0.0, "trimmed vector parameter");
  }
  else if (strcmp(method,"slide")==0) {
    cmd_line->AddOption("-r", &RADIUS, 1, "radius");
    cmd_line->AddOption("-d", &DISTANCE, 1, "distance");
  }
  else if (strcmp(method,"sparse")==0) {
    cmd_line->AddOption("-d", &OFFSET, 0, "feature offset");
    cmd_line->AddOption("-0", &ZEROES, false, "do not omit zeros, only NaN values");
  }
  else if (strcmp(method,"test")==0) {
    cmd_line->AddOption("-g", &GREATER, false, "test if greater than");
    cmd_line->AddOption("-e", &EQUAL, false, "test if equal");
    cmd_line->AddOption("-c", &CUTOFF, 0, "cutoff value");
  }

  return cmd_line;
}



   
	

//---------------------------------------------------------------------------------//
// CLASS VectorBuffer                                                              //
//---------------------------------------------------------------------------------//

class VectorBuffer
{
 public:
  VectorBuffer(char *file, long int buffer_size);
  VectorBuffer(long int buffer_size, bool load_stdin=false);
  ~VectorBuffer();
  
  // func
  bool ReadNext();
  void PrintLabel();
  
  // data
  FileBufferText *buffer;
  long int n_vectors;
  long int index;
  bool labels;
  string vec_label;
  long int vec_size;
  double *vec;
};

    
	
//-----Constructor-----
//
VectorBuffer::VectorBuffer(char *file, long int buffer_size)
{
  buffer = new FileBufferText(file,buffer_size);
  n_vectors = buffer->CountLines();
  vec = NULL;
  index = 0;
}


//-----Constructor-----
//
VectorBuffer::VectorBuffer(long int buffer_size, bool load_stdin)
{
  if (load_stdin==true) { 
    FILE *fp = LoadStdIn(&n_vectors,buffer_size);
    buffer = new FileBufferText(fp,buffer_size);
  }
  else {
    n_vectors = 1;
    buffer = new FileBufferText((char *)NULL,buffer_size);
  }
  vec = NULL;
  index = 0;
}



//-----Destructor-----
//
VectorBuffer::~VectorBuffer()
{
  if (buffer!=NULL) delete buffer;
  if (vec!=NULL) delete vec;
}



//-----ReadNext-----
//
bool VectorBuffer::ReadNext()
{
  if (vec!=NULL) delete vec;
  vec = NULL;
  char *inp = buffer->Next();
  if (inp==NULL) return false;
  if (strchr(inp,'\t')!=NULL) {
    if (index==0) labels = true;
    else if (labels==false) { fprintf(stderr, "Line %ld: no label please!\n", index+1); exit(1); }
    vec_label = GetNextToken(&inp,'\t');
  }
  else {
    if (index==0) labels = false;
    else if (labels==true) { fprintf(stderr, "Line %ld: expected label here!\n", index+1); exit(1); }
  }
  vec_size = VEC_SIZE;
  vec = ReadDoubleVec(inp,&vec_size);
  ++index;
  return true;
}



//-----PrintLabel-----
//
void VectorBuffer::PrintLabel()
{
  if (labels==true) cout << vec_label << '\t';
}



//---------------------------------------------------------------------------------//
// END CLASS VectorBuffer                                                          //
//---------------------------------------------------------------------------------//







//---------------------------------------------------------------------------------//
// CLASS Vectors                                                                   //
//---------------------------------------------------------------------------------//

class Vectors
{
 public:
  Vectors(VectorBuffer *VB);
  ~Vectors();
  
  // func
  
  // data
  long int n_vectors;
  bool labels;
  string *vec_label;
  long int *vec_size;
  double **vec;
};

    
	
//-----Constructor-----
//
Vectors::Vectors(VectorBuffer *VB)
{
  n_vectors = VB->n_vectors;
  vec_label = new string[n_vectors];
  vec_size = new long int[n_vectors];
  ALLOCATE1D(vec,n_vectors,double *);
  Progress PRG("Reading vectors...",n_vectors);
  for (long int n=0; VB->ReadNext()==true; n++) {
    vec_label[n] = VB->vec_label;
    vec_size[n] = VB->vec_size;
    ALLOCATE1D(vec[n],vec_size[n],double);
    for (long int k=0; k<vec_size[n]; k++) vec[n][k] = VB->vec[k]; 
    PRG.Check();
  }
  PRG.Done();
  labels = VB->labels;
}


//-----Destructor-----
//
Vectors::~Vectors()
{
  delete [] vec_label;
  FREE1D(vec_size);
  FREE2D(vec,n_vectors);
}



//---------------------------------------------------------------------------------//
// END CLASS Vectors                                                               //
//---------------------------------------------------------------------------------//











//---------------------------------------------------------------------------------//
// MAIN	                                                                           //
//---------------------------------------------------------------------------------//
int main(int argc, char* argv[]) 
{
  // process command-line arguments
  if (argc<2) { Usage(); exit(1); }
  char *method = argv[1];
  if (method[0]=='-') ++method;   // compatibility with previous versions
  CmdLine *cmdLine = InitCmdLine(method);
  cmdLine->Read(argv+1,argc-1);
  if (HELP) { 
    char usage_str[200];
    sprintf(usage_str, "%s %s [OPTIONS] <VECTORS>\n", argv[0], method); 
    cmdLine->Usage(usage_str); 
    exit(1); 
  }
  MESSAGES(VERBOSE);

  // load vectors
  bool load_stdin = (strcmp(method,"pairs")==0)||(strcmp(method,"permute")==0);
  VectorBuffer *VB = new VectorBuffer(BUFFER_SIZE,load_stdin);
  long int n_vectors = VB->n_vectors;

  // setup output format
  char fmt[100];
  if (strlen(FMT_STR)>0) sprintf(fmt, "%s", FMT_STR);
  else sprintf(fmt, "%%.%if", DECIMALS);
  if (VERBOSE) fprintf(stderr, "* Found %ld vectors. Using format '%s'.\n", n_vectors, fmt);

  
  //--------------------------------------------------------------------------------------//
  // OPTION -norm: normalize vectors                                                      //
  //--------------------------------------------------------------------------------------//
  if (strcmp(method,"norm")==0) { 
    Progress PRG("Normalizing rows...",n_vectors);
    while (VB->ReadNext()==true) {
      if (ALPHA>0) {
        double *v = VectorCopy(VB->vec,VB->vec_size);
        VectorSort(v,VB->vec_size);
	long int k = (long int)floor(ALPHA*VB->vec_size);
	double mean = VectorAvg(&v[k],VB->vec_size-2*k);
	double sd = VectorStd(&v[k],VB->vec_size-2*k);
        delete v;
        for (unsigned long int k=0; k<(unsigned long int)VB->vec_size; k++) VB->vec[k] = (VB->vec[k]-mean)/sd;
      }
      else {
        if (RANGE) VectorNormRange(VB->vec,VB->vec_size);
        else VectorNorm(VB->vec,VB->vec_size);
      }
      VB->PrintLabel();
      VectorPrint(VB->vec,VB->vec_size,fmt);
      printf("\n");
      PRG.Check();
    }
    PRG.Done();
  }


  //--------------------------------------------------------------------------------------//
  // OPTION -exp: compute the exponential (inverse of log)                                //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"exp")==0) { 
    Progress PRG("Computing powers...",n_vectors);
    while (VB->ReadNext()==true) {
      VB->PrintLabel();
      for (long int j=0; j<VB->vec_size; j++) {
        if (VB->vec[j]!=VB->vec[j]) printf("NaN");
        else printf(fmt, pow((double)EXP_BASE,VB->vec[j]));
	printf(" ");
      }
      printf("\n");
      PRG.Check();
    }
    PRG.Done();
  }


  //--------------------------------------------------------------------------------------//
  // OPERATION format: format vectors                                                       //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"format")==0) { 
    Progress PRG("Formatting rows...",n_vectors);
    while (VB->ReadNext()==true) {
      VB->PrintLabel();
      for (long int j=0; j<VB->vec_size; j++) {
        if (VB->vec[j]!=VB->vec[j]) printf("NaN");
        else printf(fmt, VB->vec[j]);
	printf("%c", (unsigned long int)j>=strlen(DELIMS) ? ' ' : DELIMS[j]);
      }
      printf("\n");
      PRG.Check();
    }
    PRG.Done();
  }


  //--------------------------------------------------------------------------------------//
  // OPERATION rank: compute ranks                                                        //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"rank")==0) { 
    Progress PRG("Ranking vectors...",n_vectors);
    while (VB->ReadNext()==true) {
      VB->PrintLabel();
      int *r = GetVectorRanks(VB->vec,VB->vec_size);
      if (RANK_NORM) for (long int j=0; j<VB->vec_size; j++) { printf(fmt, (double)(r[j]+1)/VB->vec_size); printf(" "); }
      else for (long int j=0; j<VB->vec_size; j++) printf("%d ", r[j]+1);
      printf("\n");
	  delete r;
      PRG.Check();
    }
    PRG.Done();
  }


  //--------------------------------------------------------------------------------------//
  // OPERATION sort: sort vectors                                                           //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"sort")==0) { 
    Progress PRG("Sorting vectors...",n_vectors);
    while (VB->ReadNext()==true) {
      VB->PrintLabel();
      VectorSort(VB->vec,VB->vec_size);
      VectorPrint(VB->vec,VB->vec_size,fmt);
      printf("\n");
      PRG.Check();
    }
    PRG.Done();
  }


  //--------------------------------------------------------------------------------------//
  // OPERATION shuffle: shuffle vectors                                                     //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"shuffle")==0) { 
    gsl_rng *RANDOM_GENERATOR = InitRandomGenerator(getpid()+time(NULL));
    Progress PRG("Shuffling vectors...",n_vectors);
    while (VB->ReadNext()==true) {
      VB->PrintLabel();
      gsl_ran_shuffle(RANDOM_GENERATOR,VB->vec,VB->vec_size,sizeof(double));  
      VectorPrint(VB->vec,VB->vec_size,fmt);
      printf("\n");
      PRG.Check();
    }
    PRG.Done();
    gsl_rng_free(RANDOM_GENERATOR);
  }


  //--------------------------------------------------------------------------------------//
  // OPTION -rev: reverse vector order                                                    //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"rev")==0) { 
    Progress PRG("Reversing...",n_vectors);
    while (VB->ReadNext()==true) {
      VB->PrintLabel();
      for (long int j=VB->vec_size-1; j>=0; j--) {
        if (VB->vec[j]!=VB->vec[j]) printf("NaN");
        else printf(fmt, VB->vec[j]);
	printf(" ");
      }
      printf("\n");
      PRG.Check();
    }
    PRG.Done();
  }


  //--------------------------------------------------------------------------------------//
  // OPTION -log: compute the logarithms                                                  //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"log")==0) { 
    Progress PRG("Computing logarithms...",n_vectors);
    while (VB->ReadNext()==true) {
      VB->PrintLabel();
      for (long int j=0; j<VB->vec_size; j++) {
        if (VB->vec[j]!=VB->vec[j]) printf("NaN");
        else printf(fmt, log(VB->vec[j]+CONST)/log(LOG_BASE));
	printf(" ");
      }
      printf("\n");
      PRG.Check();
    }
    PRG.Done();
  }


  //--------------------------------------------------------------------------------------//
  // OPTION -pow: compute the powers                                                      //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"pow")==0) { 
    Progress PRG("Computing powers...",n_vectors);
    while (VB->ReadNext()==true) {
      VB->PrintLabel();
      for (long int j=0; j<VB->vec_size; j++) {
        if (VB->vec[j]!=VB->vec[j]) printf("NaN");
        else printf(fmt, pow(VB->vec[j],(double)POW));
	printf(" ");
      }
      printf("\n");
      PRG.Check();
    }
    PRG.Done();
  }


  //--------------------------------------------------------------------------------------//
  // OPTION -fold: compute fold changes                                                   //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"fold")==0) { 
    Progress PRG("Computing fold changes...",n_vectors);
    while (VB->ReadNext()==true) {
      VB->PrintLabel();
      for (long int j=0; j<VB->vec_size; j++) { printf(fmt, VB->vec[j]); printf(" "); }
      printf("\t");
      for (long int j=1; j<VB->vec_size; j++) { printf(fmt, VB->vec[j]/VB->vec[0]-1); printf(" "); }
      for (long int j=1; j<VB->vec_size; j++) { printf(fmt, VB->vec[j]/VB->vec[j-1]-1); printf(" "); }
      printf("\n");
      PRG.Check();
    }
    PRG.Done();
  }


  //--------------------------------------------------------------------------------------//
  // OPTION -cutoff: apply cutoff filter                                                  //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"cutoff")==0) {
    Progress PRG("Applying cutoff...",n_vectors);
    while (VB->ReadNext()==true) {
      VB->PrintLabel();
      for (long int j=0; j<VB->vec_size; j++) {
        if (VB->vec[j]!=VB->vec[j]) printf("NaN");
	else if (UPPER) printf(fmt, VB->vec[j]>CUTOFF ? CUTOFF_REPLACE:VB->vec[j]);
        else printf(fmt, VB->vec[j]<CUTOFF ? CUTOFF_REPLACE:VB->vec[j]);
	printf(" ");
      }
      printf("\n");
      PRG.Check();
    }
    PRG.Done();
  }


  
  //--------------------------------------------------------------------------------------//
  // OPTION -test: test greater than                                                      //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"test")==0) {
    Progress PRG("Applying test...",n_vectors);
    while (VB->ReadNext()==true) {
      VB->PrintLabel();
      for (long int j=0; j<VB->vec_size; j++) {
        if (VB->vec[j]!=VB->vec[j]) printf("NaN");
	else printf("%i", (EQUAL&&(VB->vec[j]==CUTOFF))||(GREATER&&(VB->vec[j]>CUTOFF))||(!GREATER&&(VB->vec[j]<CUTOFF)));
	printf(" ");
      }
      printf("\n");
      PRG.Check();
    }
    PRG.Done();
  }


  //--------------------------------------------------------------------------------------//
  // OPTION -slide: sum over sliding window                                               //
  //--------------------------------------------------------------------------------------//

  else if (strcmp(method,"slide")==0) {
    Progress PRG("Printing sliding windows...",n_vectors);
    while (VB->ReadNext()==true) {
      for (long int j=RADIUS; j+RADIUS<VB->vec_size; j+=DISTANCE) {
        VB->PrintLabel();
        for (long int z=j-RADIUS; z<=j+RADIUS; z++) { printf(fmt, VB->vec[z]); printf(" "); }
        printf("\n");
      }
      PRG.Check();
    }
    PRG.Done();
  }

  
  //--------------------------------------------------------------------------------------//
  // OPTION -n: vector size                                                               //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"n")==0) {  
    Progress PRG("Computing sizes...",n_vectors);
    while (VB->ReadNext()==true) {
      VB->PrintLabel();
      printf("%ld\n", VB->vec_size); 
      PRG.Check();
    }
    PRG.Done();
  }
  

  //--------------------------------------------------------------------------------------//
  // OPTION -stat: vector stat                                                            //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"stat")==0) {  
    Progress PRG("Computing statistics...",n_vectors);
    while (VB->ReadNext()==true) {
      VB->PrintLabel();
      printf("%ld", VB->vec_size); printf("\t");
      printf(fmt, VectorSum(VB->vec,VB->vec_size)); printf("\t");
      printf(fmt, VectorMin(VB->vec,VB->vec_size)); printf("\t"); 
      printf(fmt, VectorMax(VB->vec,VB->vec_size)); printf("\t"); 
      printf(fmt, VectorAvg(VB->vec,VB->vec_size)); printf("\t"); 
      printf(fmt, VectorStd(VB->vec,VB->vec_size)); printf("\t");
      VectorSort(VB->vec,VB->vec_size);
      printf(fmt, gsl_stats_median_from_sorted_data(VB->vec,1,VB->vec_size)); 
      printf("\n");
      PRG.Check();
    }
    PRG.Done();
  }
  

  //--------------------------------------------------------------------------------------//
  // OPTION -sum: vector sum                                                              //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"sum")==0) {  
    Progress PRG("Computing sum...",n_vectors);
    while (VB->ReadNext()==true) {
      VB->PrintLabel();
      printf(fmt, VectorSum(VB->vec,VB->vec_size)); 
      printf("\n");
      PRG.Check();
    }
    PRG.Done();
  }
  

  //--------------------------------------------------------------------------------------//
  // OPTION -div: divide                                                                  //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"div")==0) {  
    Progress PRG("Dividing...",n_vectors);
    while (VB->ReadNext()==true) {
      VB->PrintLabel();
      for (long int k=1; k<VB->vec_size; k++) { printf(fmt, VB->vec[k]/VB->vec[0]); printf(" "); }
      printf("\n");
      PRG.Check();
    }
    PRG.Done();
  }
  

  //--------------------------------------------------------------------------------------//
  // OPTION -min: vector min                                                              //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"min")==0) {  
    Progress PRG("Computing min...",n_vectors);
    while (VB->ReadNext()==true) {
      VB->PrintLabel();
      printf(fmt, VectorMin(VB->vec,VB->vec_size)); 
      printf("\n");
      PRG.Check();
    }
    PRG.Done();
  }
  
  
  //--------------------------------------------------------------------------------------//
  // OPTION -imin: vector min index                                                       //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"imin")==0) {  
    Progress PRG("Computing min index...",n_vectors);
    while (VB->ReadNext()==true) {
      VB->PrintLabel();
      printf("%ld\n", (long int)gsl_stats_min_index(VB->vec,1,VB->vec_size)); 
      PRG.Check();
    }
    PRG.Done();
  }
  
  
  //--------------------------------------------------------------------------------------//
  // OPTION -max: vector max                                                              //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"max")==0) {  
    Progress PRG("Computing max...",n_vectors);
    while (VB->ReadNext()==true) {
      VB->PrintLabel();
      printf(fmt, VectorMax(VB->vec,VB->vec_size)); 
      printf("\n");
      PRG.Check();
    }
    PRG.Done();
  }
  

  //--------------------------------------------------------------------------------------//
  // OPTION -absmax: vector absolute max                                                  //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"absmax")==0) {  
    Progress PRG("Computing absmax...",n_vectors);
    while (VB->ReadNext()==true) {
      VB->PrintLabel();
      double max = VectorMax(VB->vec,VB->vec_size);
      double min = VectorMin(VB->vec,VB->vec_size);
      printf(fmt, fabs(max)>fabs(min)?max:min); 
      printf("\n");
      PRG.Check();
    }
    PRG.Done();
  }
  

  //--------------------------------------------------------------------------------------//
  // OPTION -imax: vector max index                                                       //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"imax")==0) {  
    Progress PRG("Computing max index...",n_vectors);
    while (VB->ReadNext()==true) {
      VB->PrintLabel();
      printf("%ld\n", (long int)gsl_stats_max_index(VB->vec,1,VB->vec_size)); 
      PRG.Check();
    }
    PRG.Done();
  }
  
  
  //--------------------------------------------------------------------------------------//
  // OPTION -med: vector median                                                           //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"med")==0) {  
    Progress PRG("Computing median...",n_vectors);
    while (VB->ReadNext()==true) {
      VB->PrintLabel();
      VectorSort(VB->vec,VB->vec_size);
      printf(fmt, gsl_stats_median_from_sorted_data(VB->vec,1,VB->vec_size)); 
      printf("\n");
      PRG.Check();
    }
    PRG.Done();
  }
  
  
  //--------------------------------------------------------------------------------------//
  // OPTION -q: vector quantiles                                                          //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"q")==0) {  
    Progress PRG("Computing quantiles...",n_vectors);
    while (VB->ReadNext()==true) {
      VB->PrintLabel();
      VectorSort(VB->vec,VB->vec_size);
      printf(fmt, gsl_stats_quantile_from_sorted_data(VB->vec,1,VB->vec_size,0.25)); printf(" ");
      printf(fmt, gsl_stats_quantile_from_sorted_data(VB->vec,1,VB->vec_size,0.50)); printf(" ");
      printf(fmt, gsl_stats_quantile_from_sorted_data(VB->vec,1,VB->vec_size,0.75)); printf("\n");
      PRG.Check();
    }
    PRG.Done();
  }
  
  
  //--------------------------------------------------------------------------------------//
  // OPTION -mean: vector mean                                                            //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"m")==0) {  
    Progress PRG("Computing mean...",n_vectors);
    while (VB->ReadNext()==true) {
      VB->PrintLabel();
      if (ALPHA>0) {
        VectorSort(VB->vec,VB->vec_size);
	long int k = (long int)floor(ALPHA*VB->vec_size);
	printf(fmt, VectorAvg(&VB->vec[k],VB->vec_size-2*k));
      }
      else printf(fmt, VectorAvg(VB->vec,VB->vec_size)); 
      printf("\n");
      PRG.Check();
    }
    PRG.Done();
  }
  
 
  //--------------------------------------------------------------------------------------//
  // OPTION -sd: vector sd                                                                //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"sd")==0) {  
    Progress PRG("Computing sd...",n_vectors);
    while (VB->ReadNext()==true) {
      VB->PrintLabel();
      if (ALPHA>0) {
        VectorSort(VB->vec,VB->vec_size);
	long int k = (long int)floor(ALPHA*VB->vec_size);
	printf(fmt, VectorStd(&VB->vec[k],VB->vec_size-2*k));
      }
      else printf(fmt, VectorStd(VB->vec,VB->vec_size)); 
      printf("\n");
      PRG.Check();
    }
    PRG.Done();
  }
    

  //--------------------------------------------------------------------------------------//
  // OPTION -sparse: convert to sparse                                                    //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"sparse")==0) { 
  
    Progress PRG("Converting to sparse format...",n_vectors);
    while (VB->ReadNext()==true) {
      VB->PrintLabel();
      for (long int j=0; j<VB->vec_size; j++) 
        if ((VB->vec[j]==VB->vec[j])&&(ZEROES||(VB->vec[j]!=0))) { printf("%ld:", OFFSET+j); printf(fmt, VB->vec[j]); printf(" "); }
      printf("\n");
      PRG.Check();
    }
    PRG.Done();

  }
  
  
  //--------------------------------------------------------------------------------------//
  // OPTION -items: convert to itemsets                                                   //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"items")==0) { 
    while (VB->ReadNext()==true) {
      VB->PrintLabel();
      for (long int j=0; j<VB->vec_size; j++) 
        if ((VB->vec[j]==VB->vec[j])&&(VB->vec[j]!=0)) printf("%ld ", (long int)NBINS*j+(long int)VB->vec[j]);
      printf("\n");
    }
  }
  
  
  //--------------------------------------------------------------------------------------//
  // OPTION -pairs: print all vector pairs                                                //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"pairs")==0) { 
    Vectors V(VB);
    Progress PRG("Print vectors pairs...",n_vectors);
    for (long int i1=0; i1<n_vectors; i1++) {
      for (long int i2=i1+1; i2<n_vectors; i2++) {
        if (V.labels==true) cout << V.vec_label[i1] << '|' << V.vec_label[i2] << '\t';
        VectorPrint(V.vec[i1],V.vec_size[i1],fmt);
        printf("\t");
        VectorPrint(V.vec[i2],V.vec_size[i2],fmt);
        printf("\n");
      }
      PRG.Check();
    }
    PRG.Done();
  }

  
  //--------------------------------------------------------------------------------------//
  // OPTION -permute: permute vectors' order                                              //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"permute")==0) { 
    Vectors V(VB);
    gsl_rng *RANDOM_GENERATOR = InitRandomGenerator(getpid()+time(NULL));
    long int *order;
    ALLOCATE1D(order,n_vectors,long int);
    for (long int i=0; i<n_vectors; i++) order[i] = i;
    gsl_ran_shuffle(RANDOM_GENERATOR,order,n_vectors,sizeof(long int));  
    Progress PRG("Printing vectors in permuted order...",n_vectors);
    for (long int i=0; i<n_vectors; i++) {
      if (V.labels==true) cout << V.vec_label[order[i]] << '\t';
      VectorPrint(V.vec[order[i]],V.vec_size[order[i]],fmt);
      printf("\n");
      PRG.Check();
    }
    PRG.Done();
    gsl_rng_free(RANDOM_GENERATOR);
  }

  
  //--------------------------------------------------------------------------------------//
  // OPTION -hist: histogram                                                              //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"hist")==0) { 
    Progress PRG("Computing histogram...",n_vectors);
    long int n_bins = HIST_NBINS;
    double min = HIST_MIN;
    double max = HIST_MAX;
    double *bins;
    ALLOCATE1D_INIT(bins,n_bins,double,0);
    while (VB->ReadNext()==true) {
      double *X = VB->vec;
      long int n = VB->vec_size;
      long int N = 0;
      for (long int k=0; k<n; k++) {
        double val = (X[k]-min) / (max-min);
        if ((val<0)||(val>=1)) continue;
        int b = (int)(n_bins*val);
        if ((b<0)||(b>=n_bins)) { fprintf(stderr, "Error: overflow!\n"); exit(1); }
        bins[b]++;
        N++;
      }  
      PRG.Check();
    }
    PRG.Done();
    double *freq;
    ALLOCATE1D(freq,n_bins,double);
    for (long int b=0; b<n_bins; b++) freq[b] = bins[b];
    VectorNormSum(freq,n_bins);
    for (long int b=0; b<n_bins; b++) { printf(fmt, min+b*(max-min)/n_bins); cout << '\t'; printf(fmt, freq[b]); printf("\t%.0f\n", bins[b]); }
    FREE1D(bins);
    FREE1D(freq);
  }


  //--------------------------------------------------------------------------------------//
  // OPTION -bins: binning vector values                                                  //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"bins")==0) { 
    Progress PRG("Computing vector histogram...",n_vectors);
    while (VB->ReadNext()==true) {
      double min, max;
      long int n_bins;
      min = HIST_AUTOMIN?VectorMin(VB->vec,VB->vec_size):HIST_MIN;
      max = HIST_AUTOMAX?VectorMax(VB->vec,VB->vec_size)+HIST_BIN_SIZE:HIST_MAX;
      n_bins = (HIST_BIN_SIZE>0)?(max-min)/HIST_BIN_SIZE:HIST_NBINS;
      double *bins;
      ALLOCATE1D_INIT(bins,n_bins,double,0);
      double *X = VB->vec;
      long int n = VB->vec_size;
      long int N = 0;
      for (long int k=0; k<n; k++) {
        double val = (X[k]-min) / (max-min);
        if ((val<0)||(val>=1)) continue;
        int b = (int)(n_bins*val);
        if (b<0) { fprintf(stderr, "Error (line %ld): underflow (val=%f, min=%f, max=%f, n_bins=%ld)!\n", VB->buffer->n_line, val, min, max, n_bins); exit(1); }
        if (b>=n_bins) { fprintf(stderr, "Error (line %ld): overflow (val=%f, min=%f, max=%f, n_bins=%ld)!\n", VB->buffer->n_line, val, min, max, n_bins); exit(1); }
        bins[b]++;
        N++;
      }  
      VB->PrintLabel();
      if (HIST_NBINS_COMBINE==1) for (long int b=0; b<n_bins; b++) printf("%.0f%c", bins[b], b==n_bins-1?'\n':' ');
      else {
        if (n_bins<HIST_NBINS_COMBINE) printf("\n");
        else {
          for (long int b=0; b<=n_bins-HIST_NBINS_COMBINE; b++) {
            long int sum = 0;
            for (int bb=0; bb<HIST_NBINS_COMBINE; bb++) sum += (long int)floor(bins[b+bb]);
            printf("%ld%c", sum, b==n_bins-HIST_NBINS_COMBINE?'\n':' ');
          }
        }
      }
      FREE1D(bins);
      PRG.Check();
    }
    PRG.Done();
  }


  //--------------------------------------------------------------------------------------//
  // OPTION -merge: merge consecutive vectors                                             //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"merge")==0) {
    Progress PRG("Merging vectors...",n_vectors);
    string label;
    bool first = true;
    while (VB->ReadNext()==true) {
      if (first) { VB->PrintLabel(); label = VB->vec_label; first = false; }
      else if (VB->vec_label!=label) { printf("\n"); VB->PrintLabel(); label = VB->vec_label; }
      else printf("%s", SEPARATOR);
      VectorPrint(VB->vec,VB->vec_size,fmt);
      PRG.Check();
    }
    printf("\n");
    PRG.Done();
  }


  else {
    fprintf(stderr, "Error: '%s' is an invalid method!\n", method);
    exit(1);
  }


  // clean up
  delete cmdLine;
  delete VB;

  return 0;
}



