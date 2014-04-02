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
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <string>
#include <map>
#include <iostream>
#include "core.h"


using namespace std;
typedef map<string,pair<long int,long int> > KeyMap;
class StringSets;
typedef double* (*CalcStatisticFuncPtr)(StringSets*,bool);
extern "C" pid_t getpid();




// options
bool GAUSSIAN = false;
bool VERBOSE;
long int N_PERMUTATIONS;
bool DETAILS;
bool HEADER;
float QVAL_CUTOFF;
bool NORMALIZE;
bool UNDER;
bool APPROX;
bool PRINT_FDR;
char *STATISTIC;
long int MIN_SUPPORT;
long int MAX_SUPPORT;



//-------InitCmdLine-----------
//
CmdLine *InitCmdLine()
{
  CmdLine *cmd_line = new CmdLine(); 

  // common options
  cmd_line->AddOption("-v", &VERBOSE, false, "verbose mode");
  cmd_line->AddOption("-kmin", &MIN_SUPPORT, 10, "minimum support per category");
  cmd_line->AddOption("-kmax", &MAX_SUPPORT, 0, "maximum support per category (default = no maximum)");
  cmd_line->AddOption("-norm", &NORMALIZE, false, "normalize row values (if applicable)");
  cmd_line->AddOption("-S", &STATISTIC, "sum", "choose statistic [sum|n|sens|spec|ratio|t|corr]");
  cmd_line->AddOption("-a", &APPROX, false, "use a distribution for p-value approximation (not applicable to all statistics)");
  cmd_line->AddOption("-u", &UNDER, false, "find depleted categories (default = enriched)");
  //cmd_line->AddOption("-g", &GAUSSIAN, false, "use gaussian distribution for p-value approximation");
  cmd_line->AddOption("-p", &N_PERMUTATIONS, 100, "number of random permutations");
  cmd_line->AddOption("-q", &QVAL_CUTOFF, 1.0, "FDR cutoff");
  cmd_line->AddOption("-f", &PRINT_FDR, false, "print FDR instead of adjusted p-values");
  cmd_line->AddOption("-h", &HEADER, false, "print header");
  cmd_line->AddOption("-d", &DETAILS, false, "print details");

  return cmd_line;
}







//---------------------------------------------------------------------------------//
// CLASS StringSets                                                                //
//---------------------------------------------------------------------------------//

class StringSets
{
 public:
  StringSets(char *file, char *vec_file=NULL);
  ~StringSets();

  // functions
  void Print();
  void PrintGOGenes(long int c);
  void Permute(gsl_rng *rnd_generator);
  double *CalcSumStatistic(bool approx=false);
  double *CalcRatioStatistic(bool approx=false);
  double *CalcTStatistic(bool approx=false);
  double *CalcSensitivityStatistic(bool approx=false);
  double *CalcSpecificityStatistic(bool approx=false);
  double *CalcHyperGeomStatistic(bool approx=false);
  double *CalcCorrStatistic(bool approx=false);
  double *RunPermutations(double *Y, long int n_permutations, CalcStatisticFuncPtr calc_statistic);			// run permutations to determine p-values
  double *RunPermutations(double *Y, long int n_permutations, CalcStatisticFuncPtr calc_statistic, bool gaussian);	// run permutations to determine p-values
  double *RunApproxPermutations(double *Y, long int n_permutations, CalcStatisticFuncPtr calc_statistic); 		// run permutations to determine q-values only

  // data
  long int n_rows, n_cols;
  char **ROW_LABELS, **COL_LABELS;
  long int *COL_STATS;
  KeyMap key_map;
  int **B;
  bool use_totals;
  long int n_values;
  float *V, *Vtotal;
  double Vsum, VsumZ, Vsum2, Vtotal_sum;
};


//----- Constructor -------
//
StringSets::StringSets(char *file, char *vec_file)
{
  // initialize
  FileBufferText buffer(file);
  n_rows = buffer.CountLines();
  ALLOCATE1D(ROW_LABELS,n_rows,char *);
  n_cols = 0;
  n_values = 0;
  use_totals = false;
    
  // parse keys
  Progress PRG("Reading labels...",n_rows);
  KeyMap key_map0;
  for (long int r=0; r<n_rows; r++) {
    char *inp = buffer.Next();
    ROW_LABELS[r] = StrCopy(GetNextToken(&inp,'\t'));
    if (vec_file==NULL) GetNextToken(&inp,'\t');
    while (inp[0]!=0) {
      char *key = GetNextToken(&inp," ");
      if (key_map0.find(key)==key_map0.end()) { key_map0[key].first = n_cols++; key_map0[key].second = 1; }
      else key_map0[key].second += 1;
    }
    PRG.Check();
  }
  PRG.Done();

  // enforce minimum support criteria
  n_cols = 0;
  if (MAX_SUPPORT==0) MAX_SUPPORT = n_rows;
  for (KeyMap::iterator x=key_map0.begin(); x!=key_map0.end(); x++)
    if ((x->second.second>=MIN_SUPPORT)&&(x->second.second<=MAX_SUPPORT)) { key_map[x->first].first = n_cols++; key_map[x->first].second = x->second.second; }
  if (VERBOSE) fprintf(stderr, "* Found %ld rows and %ld columns.\n", n_rows, n_cols);

  // create table
  ALLOCATE2D(B,n_cols,n_rows+1,int);
  for (long int c=0; c<n_cols; c++) B[c][0] = 0;
  ALLOCATE1D(V,n_rows,float);
  ALLOCATE1D(Vtotal,n_rows,float);
  ALLOCATE1D_INIT(COL_LABELS,n_cols,char *,NULL);
  ALLOCATE1D_INIT(COL_STATS,n_cols,long int,0);
  buffer.Reset();
  Progress PRG2("Creating table...",n_rows);
  for (long int r=0; r<n_rows; r++) {
    char *inp = buffer.Next();
    GetNextToken(&inp,'\t');
    if (vec_file==NULL) { 
      char *v_str = GetNextToken(&inp,'\t');
      int n_tokens = CountTokens(v_str,' '); 
      if ((n_tokens==0)||(n_tokens>3)) { fprintf(stderr, "Line %ld: 2nd column should contain 1 or 2 values!\n", r+1); exit(1); }
      if (r==0) n_values = n_tokens;
      else if (n_tokens!=n_values) { fprintf(stderr, "Line %ld: expected %ld instead of %d tokens in 2nd column!\n", r+1, n_values, n_tokens); exit(1); }
      V[r] = atof(GetNextToken(&v_str,' '));
      Vtotal[r] = n_values==2 ? atof(GetNextToken(&v_str,' ')) : 1;
    }
    while (inp[0]!=0) {
      char *key = GetNextToken(&inp," ");
      if (key_map.find(key)!=key_map.end()) {
        long int c = key_map[key].first;
        if (COL_LABELS[c]==NULL) COL_LABELS[c] = StrCopy(key);
	B[c][0]++;
	B[c][B[c][0]] = r;
        COL_STATS[c]++;
      }
    }
    PRG2.Check();
  }
  PRG2.Done();

  // load vector
  if (vec_file!=NULL) {
    long int n_vec_rows;
    float **X = LoadMatrix(vec_file,&n_vec_rows,&n_values);
    if (VERBOSE) fprintf(stderr, "* Found a %ldx%ld matrix.\n", n_vec_rows, n_values); 
    if ((n_vec_rows!=n_rows)||(n_values>2)) { fprintf(stderr, "Wrong dimensions!\n"); exit(1); }
    for (long int r=0; r<n_rows; r++) V[r] = X[r][0];
    if (n_values==2) for (long int r=0; r<n_rows; r++) Vtotal[r] = X[r][1];
    FREE2D(X,n_vec_rows);
  }

  // compute normalized value if applicable
  if ((n_values==2)&&(NORMALIZE==true)) {
    use_totals = false;
    for (long int r=0; r<n_rows; r++) { V[r] /= Vtotal[r]; Vtotal[r] = 1; }
  }
  else use_totals = true;

  // compute vector sums
  Vsum = VsumZ = Vsum2 = Vtotal_sum = 0;
  if (use_totals==false) for (long int r=0; r<n_rows; r++) { Vsum += V[r]; Vsum2 += V[r]*V[r]; Vtotal_sum += Vtotal[r]; }
  else for (long int r=0; r<n_rows; r++) { Vsum += V[r]; VsumZ += V[r]/Vtotal[r]; Vsum2 += pow((double)(V[r]/Vtotal[r]),2.0); Vtotal_sum += Vtotal[r]; }
 
  if (VERBOSE) fprintf(stderr, "* using normalized values = %s\n", use_totals?"NO":"YES");
}



//----- Destructor -------
//
StringSets::~StringSets()
{
  FREE2D(B,n_cols);
  FREE1D(V);
  FREE1D(Vtotal);
  FREE2D(ROW_LABELS,n_rows);
  FREE2D(COL_LABELS,n_cols);
  FREE1D(COL_STATS);
}



//----- Print -------
//
void StringSets::Print()
{
  // print table
  Progress PRG("Printing...",n_rows);
  printf("\t"); for (long int c=0; c<n_cols; c++) printf("%s ", COL_LABELS[c]); printf("\n");
/*  for (long int r=0; r<n_rows; r++) {
    printf("%s\t", ROW_LABELS[r]);
    for (long int c=0; c<n_cols; c++) printf("%i ", A[c][r]);
    printf("\n");
    PRG.Check();
  }*/
  PRG.Done();
}



//----- PrintGOGenes -------
//
void StringSets::PrintGOGenes(long int c)
{
  // print gene list for a given GO
  for (long int k=1; k<=B[c][0]; k++) printf("%s ", ROW_LABELS[B[c][k]]);
}



//----- Permute -------
//
void StringSets::Permute(gsl_rng *rnd_generator)
{
  if (use_totals==false) gsl_ran_shuffle(rnd_generator,V,n_rows,sizeof(float));
  else {
    long int *x = new long int[n_rows];
    for (long int r=0; r<n_rows; r++) x[r] = r;
    gsl_ran_shuffle(rnd_generator,x,n_rows,sizeof(long int));
    float *Q = new float[n_rows];
    float *Qtotal = new float[n_rows];
    for (long int r=0; r<n_rows; r++) { Q[r] = V[x[r]]; Qtotal[r] = Vtotal[x[r]]; }
    delete V;
    delete Vtotal;
    V = Q;
    Vtotal = Qtotal;
    delete x;
  }
}




//---- CalcTStatistic ------
//
double *StringSets::CalcTStatistic(bool approx)
{
  double *Y;
  ALLOCATE1D(Y,n_cols,double);

  if (use_totals==false) {   
    for (long int c=0; c<n_cols; c++) {
      double mean[2] = {0,0};
      double var[2] = {0,0};
      long int n[2] = {0,0};
      for (long int z=1; z<=B[c][0]; z++) {
        long int r = B[c][z];
        n[1]++;
        mean[1] += V[r];
        var[1] += V[r]*V[r];
      }
      n[0] = n_rows - n[1];
      mean[0] = Vsum - mean[1];
      var[0] = Vsum2 - var[1];
      for (int k=0; k<=1; k++) { mean[k] /= n[k]; var[k] = var[k]/n[k] - mean[k]*mean[k]; }     
      Y[c] = (mean[1]-mean[0])/sqrt(var[1]/n[1]+var[0]/n[0]); 
      if (UNDER) Y[c] = -Y[c];
      if (approx==true) {
        long int df = (long int)floor(pow(var[0]/n[0]+var[1]/n[1],2.0)/(pow(var[0]/n[0],2.0)/(n[0]-1)+pow(var[1]/n[1],2.0)/(n[1]-1)));
        Y[c] = (df<0)?1.0:gsl_cdf_tdist_Q(Y[c],df);
      }
    }
  }
  else {
    for (long int c=0; c<n_cols; c++) {
      long int n[2] = {0,0};
      double sum[2] = {0,0};
      double total[2] = {0,0};
      double mean[2] = {0,0};
      double sumZ[2] = {0,0};
      double sumqZ[2] = {0,0};
      double varZ[2] = {0,0};
      for (long int z=1; z<=B[c][0]; z++) {
        long int r = B[c][z];
        n[1]++;
        sum[1] += V[r];
        total[1] += Vtotal[r];
        sumZ[1] += V[r]/Vtotal[r];
        sumqZ[1] += pow((double)V[r]/Vtotal[r],2.0);
      }
      n[0] = n_rows - n[1];
      total[0] = Vtotal_sum - total[1];
      sum[0] = Vsum - sum[1];
      sumZ[0] = VsumZ - sumZ[1];
      sumqZ[0] = Vsum2 - sumqZ[1];
      for (int k=0; k<=1; k++) { mean[k] = sum[k]/total[k]; varZ[k] = sumqZ[k]/n[k] - pow((double)sumZ[k]/n[k],2.0); }     
      Y[c] = (mean[1]-mean[0])/sqrt(varZ[1]/n[1]+varZ[0]/n[0]); 
      if (UNDER) Y[c] = -Y[c];
      if (approx==true) {
        long int df = (long int)floor(pow(varZ[0]/n[0]+varZ[1]/n[1],2.0)/(pow(varZ[0]/n[0],2.0)/(n[0]-1)+pow(varZ[1]/n[1],2.0)/(n[1]-1)));
        Y[c] = (df<0)?1.0:gsl_cdf_tdist_Q(Y[c],df);
      }
    }
  }

  return Y;
}




//---- CalcSpecificityStatistic ------
//
double *StringSets::CalcSpecificityStatistic(bool approx)
{
  double *Y;
  ALLOCATE1D(Y,n_cols,double);
  long int t = 0;
  for (long int r=0; r<n_rows; r++) t += UNDER?V[r]<0:V[r]>0;
  for (long int c=0; c<n_cols; c++) {
    long int k = 0;
    for (long int z=1; z<=B[c][0]; z++) {
      long int r = B[c][z];
      k += UNDER?V[r]<0:V[r]>0;
    }
    Y[c] = (double)k/t;
    if (approx==true) { cerr << "Error: not implemented yet!\n"; exit(1); }
  }

  return Y;
}




//---- CalcSensitivityStatistic ------
//
double *StringSets::CalcSensitivityStatistic(bool approx)
{
  double *Y;
  ALLOCATE1D(Y,n_cols,double);
  for (long int c=0; c<n_cols; c++) {
    long int n1 = COL_STATS[c];
    long int k = 0;
    for (long int z=1; z<=B[c][0]; z++) {
      long int r = B[c][z];
      k += UNDER?V[r]<0:V[r]>0;
    }
    Y[c] = (double)k/n1;
    if (approx==true) { cerr << "Error: not implemented yet!\n"; exit(1); }
  }

  return Y;
}




//---- CalcHyperGeomStatistic ------
//
double *StringSets::CalcHyperGeomStatistic(bool approx)
{
  double *Y;
  ALLOCATE1D(Y,n_cols,double);
  long int t = 0;
  for (long int r=0; r<n_rows; r++) t += UNDER?V[r]<0:V[r]>0;
  for (long int c=0; c<n_cols; c++) {
    long int n1 = COL_STATS[c];
    long int k = 0;
    for (long int z=1; z<=B[c][0]; z++) {
      long int r = B[c][z];
      k += UNDER?V[r]<0:V[r]>0;
    }
    Y[c] = k;
    if (approx==true) Y[c] = k==0?1.0:(double)gsl_cdf_hypergeometric_Q(k-1,n1,n_rows-n1,t);
  }

  return Y;
}




//---- CalcRatioStatistic ------
//
double *StringSets::CalcRatioStatistic(bool approx)
{
  double *Y;
  ALLOCATE1D(Y,n_cols,double);

  if (use_totals==false) {   
    for (long int c=0; c<n_cols; c++) {
      double mean[2] = {0,0};
      double var[2] = {0,0};
      long int n[2] = {0,0};
      for (long int z=1; z<=B[c][0]; z++) {
        long int r = B[c][z];
        n[1]++;
        mean[1] += V[r];
        var[1] += V[r]*V[r];
      }
      n[0] = n_rows - n[1];
      mean[0] = Vsum - mean[1];
      var[0] = Vsum2 - var[1];
      //n[0] = n_rows;
      //mean[0] = Vsum;
      //var[0] = Vsum2;
      for (int k=0; k<=1; k++) { mean[k] /= n[k]; var[k] = var[k]/n[k] - mean[k]*mean[k]; }     
      Y[c] = UNDER ? mean[0]/mean[1] : mean[1]/mean[0];
      if (approx==true) {
        double m = Vsum/n_rows;
        double v = Vsum2/n_rows;
        Y[c] = gsl_cdf_ugaussian_Q((m*Y[c]-m)/sqrt(v*pow(Y[c],2.0)+v));
      }
    }
  }
  else {
    for (long int c=0; c<n_cols; c++) {
      long int n[2] = {0,0};
      double sum[2] = {0,0};
      double total[2] = {0,0};
      double mean[2] = {0,0};
      for (long int z=1; z<=B[c][0]; z++) {
        long int r = B[c][z];
        n[1]++;
        sum[1] += V[r];
        total[1] += Vtotal[r];
      }
      n[0] = n_rows - n[1];
      total[0] = Vtotal_sum - total[1];
      sum[0] = Vsum - sum[1];
      //n[0] = n_rows;
      //total[0] = Vtotal_sum;
      //sum[0] = Vsum;
      for (int k=0; k<=1; k++) mean[k] = sum[k]/total[k];
      Y[c] = UNDER ? mean[0]/mean[1] : mean[1]/mean[0];
      if (approx==true) {
        cerr << "Error: not implemented yet!\n"; exit(1);
      }
    }
  }

  return Y;
}




//---- CalcSumStatistic ------
//
double *StringSets::CalcSumStatistic(bool approx)
{
  double *Y;
  ALLOCATE1D(Y,n_cols,double);

  if (use_totals==true) { 
    for (long int c=0; c<n_cols; c++) {
      Y[c] = 0;
      long int n = B[c][0];
      double Ytotal = 0;
      for (long int r=1; r<=n; r++) { Y[c] += V[B[c][r]]; Ytotal += Vtotal[B[c][r]]; }
      Y[c] /= Ytotal;
      if (UNDER) Y[c] = -Y[c];
      if (approx==true) {
        cerr << "Error: not implemented yet!\n"; exit(1);
      }
    }
  }
  else {
    for (long int c=0; c<n_cols; c++) {
      Y[c] = 0;
      long int n = B[c][0];
      for (long int r=1; r<=n; r++) Y[c] += V[B[c][r]];
      Y[c] /= n;
      if (UNDER) Y[c] = -Y[c];
      if (approx==true) {
        cerr << "Error: not implemented yet!\n"; exit(1);
      }
    }
  }

  return Y;
}



//---- CalcCorrStatistic ------
//
double *StringSets::CalcCorrStatistic(bool approx)
{
  double *Y;
  ALLOCATE1D(Y,n_cols,double);

  if (use_totals==false) { 
    cerr << "Error: this operation is not permitted!\n"; exit(1);
  }
  else { 
    for (long int c=0; c<n_cols; c++) {
      long int n = B[c][0];
      double *v1 = new double[n];
      double *v2 = new double[n];
      for (long int r=1; r<=n; r++) { v1[r-1] = V[B[c][r]]; v2[r-1] = Vtotal[B[c][r]]; }
      Y[c] = fabs(VectorCorr(v1,v2,n));
      if (UNDER) Y[c] = 1.0-Y[c];
      if (approx==true) Y[c] = gsl_cdf_tdist_Q(Y[c]*sqrt((n-2)/(1-pow(Y[c],2))),n-2);
      delete v1;
      delete v2;
    }
  }

  return Y;
}



//---- RunPermutations ------
//
double *StringSets::RunPermutations(double *Y, long int n_permutations, CalcStatisticFuncPtr calc_statistic)
{
  gsl_rng *rnd_generator = InitRandomGenerator(getpid()+time(NULL));
  double *pval;
  ALLOCATE1D_INIT(pval,n_cols,double,0.0); 
  Progress PRG("Running random permutations...",n_permutations);
  for (long int p=0; p<n_permutations; p++) {
    Permute(rnd_generator);
    double *Y_random = (*calc_statistic)(this,false);
    for (long int c=0; c<n_cols; c++) pval[c] += (Y_random[c]>=Y[c]);
    FREE1D(Y_random);
    PRG.Check();
  }
  PRG.Done();
  for (long int c=0; c<n_cols; c++) pval[c] /= n_permutations; 
  free(rnd_generator);
  return pval;
}



//---- RunPermutationsGaussian ------
//
double *StringSets::RunPermutations(double *Y, long int n_permutations, CalcStatisticFuncPtr calc_statistic, bool gaussian)
{
  if (gaussian==false) return RunPermutations(Y,n_permutations,calc_statistic);

  gsl_rng *rnd_generator = InitRandomGenerator(getpid()+time(NULL));
  double **Y_random = new double*[n_permutations];
  Progress PRG("Running random permutations...",n_permutations);
  for (long int p=0; p<n_permutations; p++) {
    Permute(rnd_generator);
    Y_random[p] = (*calc_statistic)(this,false);
    PRG.Check();
  }
  PRG.Done();
  double *Y_random_avg = ColumnAvg(Y_random,n_permutations,n_cols);
  double *Y_random_std = ColumnStd(Y_random,n_permutations,n_cols);
  delete [] Y_random;
  double *pval = new double[n_cols];
  Progress PRG2("Computing p-values based on gaussian approximation...",n_cols);
  for (long int c=0; c<n_cols; c++) {
    pval[c] = gsl_cdf_gaussian_Q(Y[c]-Y_random_avg[c],Y_random_std[c]);
    PRG2.Check();
  }
  PRG2.Done();
  delete Y_random_avg;
  delete Y_random_std;
  free(rnd_generator);
  return pval;
}



//---- RunApproxPermutations ------
//
double *StringSets::RunApproxPermutations(double *Y, long int n_permutations, CalcStatisticFuncPtr calc_statistic)
{
  // run permutations
  gsl_rng *rnd_generator = InitRandomGenerator(getpid()+time(NULL));
  int *counts;
  ALLOCATE1D_INIT(counts,n_cols,int,0);
  Progress PRG("Running random permutations...",n_permutations);
  for (long int p=0; p<n_permutations; p++) {
    Permute(rnd_generator);
    double *Y_random = (*calc_statistic)(this,true);
    VectorSort(Y_random,n_cols);    
    for (long int z=0,c=0; (z<n_cols)&&(c<n_cols); c++) {
      while ((z<n_cols)&&(Y[z]<Y_random[c])) z++;
      if (z<n_cols) counts[z]++;
    }
    FREE1D(Y_random);
    PRG.Check();
  }
  PRG.Done();
  free(rnd_generator);

  // compute FDR values
  double *FDR = new double[n_cols];
  for (long int k=1,c=0; c<n_cols; c++,k++) {
    FDR[c] = (double)counts[c]/n_permutations/k;
    if (c+1<n_cols) counts[c+1] += counts[c];
  }
  double min_q = FDR[n_cols-1];
  for (long int c=n_cols-1; c>=0; c--) { if (FDR[c]>min_q) FDR[c] = min_q; else min_q = FDR[c]; }
  delete counts;
 
  return FDR;
}



//---------------------------------------------------------------------------------//
// END CLASS StringSets                                                            //
//---------------------------------------------------------------------------------//




//---- calc_sum_statistic ------
//
double *calc_sum_statistic(StringSets *S, bool approx)
{
  return S->CalcSumStatistic(approx);
}


//---- calc_ratio_statistic ------
//
double *calc_ratio_statistic(StringSets *S, bool approx)
{
  return S->CalcRatioStatistic(approx);
}


//---- calc_t_statistic ------
//
double *calc_t_statistic(StringSets *S, bool approx)
{
  return S->CalcTStatistic(approx);
}


//---- calc_hypergeom_statistic ------
//
double *calc_hypergeom_statistic(StringSets *S, bool approx)
{
  return S->CalcHyperGeomStatistic(approx);
}


//---- calc_specificity_statistic ------
//
double *calc_specificity_statistic(StringSets *S, bool approx)
{
  return S->CalcSpecificityStatistic(approx);
}


//---- calc_sensitivity_statistic ------
//
double *calc_sensitivity_statistic(StringSets *S, bool approx)
{
  return S->CalcSensitivityStatistic(approx);
}


//---- calc_corr_statistic ------
//
double *calc_corr_statistic(StringSets *S, bool approx)
{
  return S->CalcCorrStatistic(approx);
}


//---- run_binomial_test ------
//
double *run_binomial_test(StringSets *S)
{
  double *Y;
  ALLOCATE1D(Y,S->n_cols,double);

  if (S->use_totals==true) { 
    double p = (double)S->Vsum/S->Vtotal_sum;
    for (long int c=0; c<S->n_cols; c++) {
      long int k = 0;
      long int n = 0;
      for (long int r=1; r<=S->B[c][0]; r++) { k += (long int)floor(S->V[S->B[c][r]]); n += (long int)floor(S->Vtotal[S->B[c][r]]); }
      Y[c] = gsl_cdf_binomial_Q(k,p,n);
      if (UNDER) Y[c] = (double)1.0 - Y[c];
      Y[c] = -Y[c];
    }
  }
  else {
    cerr << "Error: not implemented yet!\n"; exit(1); 
    for (long int c=0; c<S->n_cols; c++) {
      Y[c] = 0;
      long int n = S->B[c][0];
      for (long int r=1; r<=n; r++) Y[c] += S->V[S->B[c][r]];
      Y[c] /= n;
      if (UNDER) Y[c] = -Y[c];
    }
  }

  return Y;
}







// ------------------------------------------------------------------------------------------------------------------//
// MAIN	                                                                                                             //
// ------------------------------------------------------------------------------------------------------------------//
int main(int argc, char* argv[]) 
{
  // read command-line options 
  CmdLine *cmdLine = InitCmdLine();
  int next_arg = cmdLine->Read(argv,argc);
  if (argc-next_arg<1) { 
    cmdLine->Usage("permutation_test [OPTIONS] vector(LABEL<tab>DATA<tab>CATEGORIES)\n  permutation_test [OPTIONS] vector(LABEL<tab>CATEGORIES) vector(DATA)"); 
    exit(1); 
  }
  MESSAGES(VERBOSE); 
  
  // open input files
  char *MATRIX_FILE = argv[next_arg++];
  char *VECTOR_FILE = next_arg<argc ? argv[next_arg] : NULL;
  
  // load values and annotations
  StringSets INPUT(MATRIX_FILE,VECTOR_FILE);

  // choose statistic
  CalcStatisticFuncPtr calc_statistic;
  if (strcmp(STATISTIC,"sum")==0) calc_statistic = &calc_sum_statistic;
  else if (strcmp(STATISTIC,"n")==0) calc_statistic = &calc_hypergeom_statistic;
  else if (strcmp(STATISTIC,"sens")==0) calc_statistic = &calc_sensitivity_statistic;
  else if (strcmp(STATISTIC,"spec")==0) calc_statistic = &calc_specificity_statistic;
  else if (strcmp(STATISTIC,"ratio")==0) calc_statistic = &calc_ratio_statistic;
  else if (strcmp(STATISTIC,"t")==0) calc_statistic = &calc_t_statistic;
  else if (strcmp(STATISTIC,"corr")==0) calc_statistic = &calc_corr_statistic;
  else { fprintf(stderr, "Error: unknown statistic '%s'!\n", STATISTIC); exit(1); }

  // compute statistic and p-values
  double *VAL = (*calc_statistic)(&INPUT,false);
  double *PVAL = APPROX ? (*calc_statistic)(&INPUT,true) : INPUT.RunPermutations(VAL,N_PERMUTATIONS,calc_statistic,GAUSSIAN);

  // sort p-values
  int *R = VectorRank(PVAL,INPUT.n_cols);
  VectorSort(PVAL,INPUT.n_cols);
 
  // determine FDR values
  double *FDR;
  if (APPROX==true) {
    // run random permutations on approximate p-values to determine FDR
    FDR = INPUT.RunApproxPermutations(PVAL,N_PERMUTATIONS,calc_statistic);
  }
  else {
    // compute FDR using conservative multi-test hypothesis correction
    FDR = new double[INPUT.n_cols];
    for (long int k=1,c=0; c<INPUT.n_cols; c++,k++) FDR[c] = PVAL[c]*INPUT.n_cols/k;
    double min_q = FDR[INPUT.n_cols-1];
    for (long int c=INPUT.n_cols-1; c>=0; c--) { if (FDR[c]>min_q) FDR[c] = min_q; else min_q = FDR[c]; }
  }

  // compute adjusted p-values
  double *QVAL = new double[INPUT.n_cols];
  QVAL[0] = 0;
  for (long int c=1; c<INPUT.n_cols; c++) { QVAL[c] = (c+1)*FDR[c]-c*FDR[c-1]; if (QVAL[c]<QVAL[c-1]) QVAL[c] = QVAL[c-1]; if (QVAL[c]>1) QVAL[c] = 1; }

  // print results
  if (HEADER) printf("CATEGORY\tCATEGORY-SIZE\tQ-VALUE\tP-VALUE\tSTATISTIC\n");
  for (long int c=0; c<INPUT.n_cols; c++) {
    if (QVAL[c]>QVAL_CUTOFF) break;
    printf("%s\t%ld\t%.2e\t%.2e\t%f", INPUT.COL_LABELS[R[c]], INPUT.COL_STATS[R[c]], PRINT_FDR?FDR[c]:QVAL[c], PVAL[c], VAL[R[c]]);
    if (DETAILS) { printf("\t"); INPUT.PrintGOGenes(R[c]); }
    printf("\n");
  }

  // cleanup
  FREE1D(VAL);
  FREE1D(PVAL);
  FREE1D(R);
  FREE1D(FDR);
  FREE1D(QVAL);
  delete cmdLine;
  
  return 0;
}



