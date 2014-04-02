//
// Copyright (c) 2011 IBM Corporation. 
// All rights reserved. This program and the accompanying materials are made available under the terms of the GNU General Public License v2.0 
// which accompanies this distribution, and is available at http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <iostream>
#include <limits>
#include <fstream>
#include <string>
#include <map>
#include <list>
#include "core.h"
#include "gzstream.h"
using namespace std;




int _MESSAGES_ = 0;







//------------------------------------------------------------------------------------------------//
// CLASS Progress                                                                                 //
//------------------------------------------------------------------------------------------------//

Progress::Progress(const char *msg, long int max_count, bool print_remaining)
{
  if (_MESSAGES_) fprintf(stderr, "%s\n", msg);
  this->max_count = max_count;
  this->print_remaining = print_remaining;
  this->count = 0;
  TIME0 = TIME = time(NULL);
}


Progress::Progress()
{
  max_count = 0;
}


Progress::~Progress()
{


}


void Progress::Init(const char *msg, long int max_count)
{
  if (max_count<=0) return;
  if (_MESSAGES_) fprintf(stderr, "%s\n", msg);
  this->max_count = max_count;
  this->count = 0;
  TIME0 = TIME = time(NULL);
}


void Progress::Check(long int c)
{
  if (max_count<=0) return;
  count += c;
  if (_MESSAGES_==0) return;
  if (time(NULL)-TIME>=1) { 
    TIME=time(NULL); 
    if (max_count>1) {
	  if (print_remaining) {
	    unsigned long int sec = (unsigned long int)(TIME-TIME0)*(max_count-count)/count;
	    fprintf(stderr, "%10.2f%% [%5lu:%02lu remaining]\r", 100.0*count/max_count, sec/60, sec%60);
	  }
      else fprintf(stderr, "%10.2f%%\r", 100.0*count/max_count); 
	}
    else fprintf(stderr, "%15ld\r", count);
  }
}


void Progress::Done()
{
  if (max_count<=0) return;
  if (_MESSAGES_==0) return;
  if (max_count>1) fprintf(stderr, "%10.2f%%\n", 100.0);
  else fprintf(stderr, "%15ld\n", count);
}



//------------------------------------------------------------------------------------------------//
// CLASS Progress                                                                                 //
//------------------------------------------------------------------------------------------------//














//------------------------------------------------------------------------------------------------//
// CLASS FileBuffer: class for reading files (text, gz, bam)                                      //
//------------------------------------------------------------------------------------------------//

//------Destructor------
//
FileBuffer::~FileBuffer()
{
  if (BUFFER!=NULL) delete BUFFER;
  if (file_name!=NULL) delete file_name;
}
 

//------CountLines-------
//
long int FileBuffer::CountLines()
{
  if (is_stdin==true) return 1;
  long int n = 0;
  Progress PRG("Counting lines in file...",1);
  while (Next()!=NULL) { n++; PRG.Check(); }
  PRG.Done();
  Reset();
  return n;
}



//------Get-------
//
char *FileBuffer::Get()
{
  return n_line>0 ? BUFFER : NULL;
}




//------------------------------------------------------------------------------------------------//
// END CLASS FileBuffer                                                                           //
//------------------------------------------------------------------------------------------------//











//------------------------------------------------------------------------------------------------//
// CLASS FileBufferText : for reading text files                                                  //
//------------------------------------------------------------------------------------------------//

//------Constructor------
//
FileBufferText::FileBufferText(const char *file, unsigned long int buffer_size)
{
  if (file==NULL) {
    is_stdin = true;
    file_ptr = stdin;
    file_name = NULL;
  }
  else {
    is_stdin = false;
    file_name = StrCopy(file);
    file_ptr = fopen(file_name,"r");
    if (file_ptr==0) { fprintf(stderr, "[FileBuffer] Error: cannot open file '%s' for reading!\n", file); exit(1); }
  }
  BUFFER_SIZE = buffer_size;
  ALLOCATE1D(BUFFER,BUFFER_SIZE,char);
  BUFFER[0] = 0;
  n_line = 0;
}
 


//------Constructor------
//
FileBufferText::FileBufferText(FILE *file_ptr, unsigned long int buffer_size)
{
  is_stdin = false;
  file_name = NULL;
  this->file_ptr = file_ptr;
  if (file_ptr==0) { fprintf(stderr, "[FileBuffer] Error: invalid file handler!\n"); exit(1); }
  BUFFER_SIZE = buffer_size;
  ALLOCATE1D(BUFFER,BUFFER_SIZE,char);
  BUFFER[0] = 0;
  n_line = 0;
}
 


//------Destructor------
//
FileBufferText::~FileBufferText()
{
  if (is_stdin==false) if (file_ptr!=NULL) fclose(file_ptr);
}
 

 
//------Reset-------
//
void FileBufferText::Reset()
{
  if (is_stdin==true) { fprintf(stderr, "[FileBufferText] Error: cannot reset standard input!\n"); exit(1); }
  rewind(file_ptr);
  n_line = 0;
}



//------Next-------
//
char *FileBufferText::Next()
{
  if ((fgets(BUFFER,BUFFER_SIZE,file_ptr)==NULL)||(feof(file_ptr)==true)) { n_line = 0; return NULL; }
  n_line++;

  size_t len = strlen(BUFFER);
  while (BUFFER[len-1]!='\n') { 
    BUFFER_SIZE *= 2;
    char *buffer = new char[BUFFER_SIZE];
    memcpy(buffer,BUFFER,len);
    delete BUFFER;
    BUFFER = buffer;
    if ((fgets(BUFFER+len,BUFFER_SIZE/2,file_ptr)==NULL)||(feof(file_ptr)==true)) { n_line = 0; return NULL; }
    len = strlen(BUFFER);
  }
  BUFFER[len-1] = 0;
  
  return BUFFER;
}




//------------------------------------------------------------------------------------------------//
// END CLASS FileBufferText                                                                       //
//------------------------------------------------------------------------------------------------//









//------------------------------------------------------------------------------------------------//
// CLASS FileBufferGZ : for reading GZ files                                                      //
//------------------------------------------------------------------------------------------------//

//------Constructor------
//
FileBufferGZ::FileBufferGZ(const char *file, unsigned long int buffer_size)
{
  if (file==NULL) { fprintf(stderr, "[FileBufferGZ] Error: cannot read gz file from standard input, use 'gunzip' instead!\n"); exit(1); }
  is_stdin = false;
  file_name = StrCopy(file);
  file_stream = new igzstream(file_name);
  BUFFER_SIZE = buffer_size;
  ALLOCATE1D(BUFFER,BUFFER_SIZE,char);
  BUFFER[0] = 0;
  n_line = 0;
}
 

//------Destructor------
//
FileBufferGZ::~FileBufferGZ()
{
  if (file_stream!=NULL) { file_stream->close(); delete file_stream; }
}
 

//------Reset-------
//
void FileBufferGZ::Reset()
{
  file_stream->close(); 
  delete file_stream;
  file_stream = new igzstream(file_name);
  n_line = 0;
}



//------Read-------
//
bool FileBufferGZ::Read(char *buffer, unsigned long int buffer_size)
{
  if (file_stream->eof()) return false;
  file_stream->getline(buffer,buffer_size-1);
  long int len = (long int)strlen(buffer);
  if (file_stream->gcount()==len+1) { buffer[len] = '\n'; buffer[len+1] = 0; }          // add <EOL> to maintain consistency with fgets()
  else if (file_stream->eof()==false) file_stream->clear();
  return true;
}



//------Next-------
//
char *FileBufferGZ::Next()
{
  if (Read(BUFFER,BUFFER_SIZE)==false) { n_line = 0; return NULL; }
  n_line++;

  size_t len = strlen(BUFFER);
  while (BUFFER[len-1]!='\n') { 
    BUFFER_SIZE *= 2;
    char *buffer = new char[BUFFER_SIZE];
    memcpy(buffer,BUFFER,len);
    delete BUFFER;
    BUFFER = buffer;
    if (Read(BUFFER+len,BUFFER_SIZE/2)==false) { n_line = 0; return NULL; }
    len = strlen(BUFFER);
  }
  BUFFER[len-1] = 0;
  
  return BUFFER;
}



//------------------------------------------------------------------------------------------------//
// END CLASS FileBufferGZ                                                                         //
//------------------------------------------------------------------------------------------------//









//------------------------------------------------------------------------------------------------//
// CLASS FileBufferBAM : for reading BAM files                                                    //
//------------------------------------------------------------------------------------------------//

//------Constructor------
//
FileBufferBAM::FileBufferBAM(const char *file, unsigned long int buffer_size)
{
  if (file==NULL) { fprintf(stderr, "[FileBufferBAM] Error: cannot read BAM file from standard input, use 'samtools view' instead!\n"); exit(1); }
  is_stdin = false;
  file_name = StrCopy(file);
  if ((samfile_ptr=samopen(file_name,"rb",NULL))==0) { fprintf(stderr, "[FileBufferBAM] Error: cannot open file '%s'!\n", file_name); exit(1); }
  header = StrCopy(samfile_ptr->header->text);
  next_header_line = header; 
  bam_ptr = bam_init1();
  BUFFER_SIZE = 0;
  BUFFER = NULL;
  n_line = 0;
}
 


//------Destructor------
//
FileBufferBAM::~FileBufferBAM()
{
  delete header;
  bam_destroy1(bam_ptr);
  samclose(samfile_ptr);
}
 

//------Reset-------
//
void FileBufferBAM::Reset()
{
  if (BUFFER!=NULL) delete BUFFER;
  bam_destroy1(bam_ptr);
  samclose(samfile_ptr);
  if ((samfile_ptr=samopen(file_name,"rb",NULL))==0) { fprintf(stderr, "[FileBufferBAM] cannot open file '%s'!\n", file_name); exit(1); }
  next_header_line = header; 
  bam_ptr = bam_init1();
  BUFFER_SIZE = 0;
  BUFFER = NULL;
  n_line = 0;
}


//------Next-------
//
char *FileBufferBAM::Next()
{
  if (next_header_line[0]!=0) {
    n_line++;
	if (BUFFER!=NULL) delete BUFFER;
    BUFFER = StrCopy(GetNextToken(&next_header_line,'\n'));
	return BUFFER;
  }
  else if (samread(samfile_ptr,bam_ptr)>=0) {
    n_line++;
	if (BUFFER!=NULL) delete BUFFER;
    BUFFER = bam_format1_core(samfile_ptr->header,bam_ptr,0);
	return BUFFER;
  }
  else { n_line = 0; return NULL; }
}




//------------------------------------------------------------------------------------------------//
// END CLASS FileBufferBAM                                                                        //
//------------------------------------------------------------------------------------------------//





//--------CreateFileBuffer----------
//
FileBuffer *CreateFileBuffer(const char *file, unsigned long int buffer_size)
{
  if (file==NULL) return new FileBufferText(file,buffer_size);
  int file_type = GetFileType(file);
  switch (file_type) {
   case -1: { fprintf(stderr, "[CreateFileBuffer] Error: cannot open file '%s'!\n", file); exit(1); }
   case 0: return new FileBufferText(file,buffer_size);
   case 1: return new FileBufferGZ(file,buffer_size);
   case 2: return new FileBufferBAM(file,buffer_size);
  }
  return NULL;
}










//------------------------------------------------------------------------------------------------//
// CLASS FuncTable                                                                                //
//------------------------------------------------------------------------------------------------//

void FuncTable::Add(char *name, FuncType *f, char *description)
{
  if (nfunc==MAX_FUNCTIONS) { printf("Maximum number of functions exceeded!\n"); return; }
  strcpy(this->name[nfunc], name);
  strcpy(this->description[nfunc], description);
  this->f[nfunc] = f;
  nfunc++;
}


void FuncTable::Execute(char *name, char *arg)
{
  int i;
  for (i=0; i<nfunc; i++) if (strcmp(this->name[i],name)==0) break;
  if (i<nfunc) f[i](arg);
  else { printf("Error: Cannot find function '%s'!\n", name); return; }
}


void FuncTable::Help()
{
  fprintf(stderr, "  Command         Description                                                         \n");
  fprintf(stderr, "--------------------------------------------------------------------------------------\n");
  for (int i=0; i<nfunc; i++) fprintf(stderr, "  %-15s %s\n", name[i], description[i]);
}


FuncType *FuncTable::Lookup(char *name)
{
  int i;
  for (i=0; i<nfunc; i++) if (strcmp(this->name[i],name)==0) break;
  if (i<nfunc) return f[i];
  else { printf("Error: Cannot find function '%s'!\n", name); exit(1); }
  return NULL;
}


//------------------------------------------------------------------------------------------------//
// CLASS FuncTable                                                                                //
//------------------------------------------------------------------------------------------------//




//------Interpreter-------
//
void Interpreter(FuncTable *FTABLE)
{
  char *inp = (char *) malloc(1000*sizeof(char));
  char *command;

  while (1) {
    printf("\n>> ");
    int i = 0;
    do { inp[i++] = getchar(); } while (inp[i-1] != '\n');
    inp[i-1] = inp[i] = 0;
    char *arg = inp;
    command = GetNextToken(&arg);
    printf("COMMAND = '%s'\nARGUMENTS = '%s'\n\n", command, arg);
    FTABLE->Execute(command, arg);
  }

  free(inp);
  return;
}









//------------------------------------------------------------------------------------------------//
// TOKENIZING                                                                                     //
//------------------------------------------------------------------------------------------------//
// TokenPreprocess  | Preprocess string before tokenizing                                         //
// CountTokens      | Count the number of tokens                                                  //
// GetNextToken     | Get the next command line token                                             //
// SkipBlank        | Skip all spaces at the beginning of the input string                        //
// Read             | Read a single value from an input string                                    //
// ReadVec          | Process a single vector from an input string                                //
//------------------------------------------------------------------------------------------------//

//-----TokenPreprocess------
//
char *TokenPreprocess(char *inp, char *inp_delims, char out_delim)
{
  if (inp==NULL) return NULL;
  char *out = StrCopy(inp);
  out[0] = 0;
  char *p = inp;
  while ((p[0]!=0)&&(strchr(inp_delims,p[0])!=NULL)) p++;				// skip input delimiters
  for (char *q=out; p[0]!=0; q++) {
    while ((p[0]!=0)&&(strchr(inp_delims,p[0])==NULL)) { q[0] = p[0]; q++; p++; }	// write to output
    while ((p[0]!=0)&&(strchr(inp_delims,p[0])!=NULL)) p++;				// skip input delimiters
    if (p[0]==0) q[0] = 0;
    else { q[0] = out_delim; }
  }
  return out;
}


//-----CountTokens------
//
int CountTokens(char *s, char TOKEN_DELIM)
{
  if (s==NULL) return 0;
  int k = 0;
  while ((s[k]!=0)&&(s[k]==' ')) k++;                         // ignore spaces
  int n = 0;
  while (1) {
    if (s[k]==0) return n;
    while ((s[k]!=0)&&(s[k]!=TOKEN_DELIM)) k++;
    if (s[k]==TOKEN_DELIM) k++;
    n++;
    while ((s[k]!=0)&&(s[k]==' ')) k++;                       // ignore spaces
    if (s[k]==0) return n;
  }
}


//-----GetNextToken------
//
char *GetNextToken(char **pBUF, const char *DELIMS)
{
  char *BUF = *pBUF;
  while (BUF&&(BUF[0]==' ')) BUF++;
  int k = 0;
  while ((BUF[k]!=0)&&(strchr(DELIMS,BUF[k])==NULL)) k++;
  if (BUF[k]==0) *pBUF = BUF+k;
  else {
    BUF[k] = 0;
    *pBUF = BUF+k+1;
  }
  return BUF;
}


//-----GetNextToken------
//
char *GetNextToken(char **pBUF, char DELIM)
{
  char *BUF = *pBUF;
  while (BUF&&(BUF[0]==' ')) BUF++;
  int k = 0;
  while ((BUF[k]!=0)&&(BUF[k]!=DELIM)) k++;
  if (BUF[k]==0) *pBUF = BUF+k;
  else {
    BUF[k] = 0;
    *pBUF = BUF+k+1;
  }
  return BUF;
}


//-----GetNextToken------
//
char *GetNextToken(char **inp)
{
  char *p = *inp;
  while (((*inp)[0]!=0)&&((*inp)[0]!=' ')) (*inp)++;
  (*inp)[0] = 0;
  do {(*inp)++;} while (((*inp)[0]!=0)&&((*inp)[0]==' ')) ;
  return p;
}


//-----SkipBlank-----
//
char *SkipBlank(char *inp)
{
  while ((inp[0]!=0)&&(inp[0]==' ')) inp++;
  return inp;
}



//-----Read<T>-----
//
void Read(char *s, float *v) { *v = atof(s); }
void Read(char *s, double *v) { *v = (double)atof(s); }
void Read(char *s, int *v) { *v = atoi(s); }
void Read(char *s, long int *v) { *v = atol(s); }
void Read(char *s, unsigned long int *v) { *v = (unsigned long int)atol(s); }
void Read(char *s, long long int *v) { *v = atoll(s); }



//-----ReadVec------
//
template <class T> void ReadVec(char *s, T **vec, long int *n, char sep)
{
  char *scopy = StrCopy(s);
  *n = CountTokens(scopy,sep);
  T *V;
  ALLOCATE1D(V,*n,T);
  char *p = scopy;
  for (long int k=0; k<*n; k++) {
    char *s = GetNextToken(&p,sep);
    LowerCase(s);
    if (strcmp(s,"nan")!=0) Read(s,&V[k]);
    else V[k] = numeric_limits<T>::quiet_NaN();
  }
  free(scopy);
  *vec = V;
}
template void ReadVec<long int>(char *s, long int **vec, long int *n, char sep);
template void ReadVec<long long int>(char *s, long long int **vec, long int *n, char sep);
template void ReadVec<double>(char *s, double **vec, long int *n, char sep);
template void ReadVec<int>(char *s, int **vec, long int *n, char sep);


//-----ReadSparseVec------
//
template <class T> void ReadSparseVec(char *inp, T **vec, long int *n)
{
  // determine vector size
  if (*n<=0) { 
    char *pp = StrCopy(inp);
    char *p = pp;
    while (p[0]!=0) {
      char *s = GetNextToken(&p,' ');
      long int k = atol(GetNextToken(&s,':'))+1;
      if (k>*n) *n = k;
    }
    free(pp);
  }

  // store vector
  char *scopy = StrCopy(inp);
  T *V;
  ALLOCATE1D(V,*n,T);
  for (long int k=0; k<*n; k++) V[k] = 0;
  char *p = scopy;
  while (p[0]!=0) {
    char *s = GetNextToken(&p,' ');
    long int k = atol(GetNextToken(&s,':'));
    if (k<0) { fprintf(stderr, "Error: negative feature number!\n"); exit(1); }
    else if (k>=*n) { fprintf(stderr, "Error: feature number (%ld) exceeds vector size (%ld)!\n", k, *n); exit(1); }
    LowerCase(s);
    if (strcmp(s,"nan")!=0) Read(s,&V[k]);
    else V[k] = numeric_limits<T>::quiet_NaN();
  }
  free(scopy);
  *vec = V;
}
template void ReadSparseVec<long int>(char *inp, long int **vec, long int *n);
template void ReadSparseVec<long long int>(char *inp, long long int **vec, long int *n);
template void ReadSparseVec<double>(char *inp, double **vec, long int *n);




//------ReadFloatVec------
//
float *ReadFloatVec(char *s, long int *n, char sep)
{
  float *V;
  ReadVec(s,&V,n,sep);
  return V;
}



//------ReadDoubleVec------
//
double *ReadDoubleVec(char *s, long int *n, char sep)
{
  double *V;
  if (strchr(s,':')==NULL) ReadVec(s,&V,n);
  else ReadSparseVec(s,&V,n);
  return V;
}



//------ReadIntVec------
//
int *ReadIntVec(char *s, long int *n, char sep)
{
  int *V;
  ReadVec(s,&V,n);
  return V;
}












//------------------------------------------------------------------------------------------------//
// VECTOR STATISTICS                                                                              //
//------------------------------------------------------------------------------------------------//
// VectorAvg     | Compute the mean value of vector v of size n                                   //
// VectorStd     | Compute the standard deviation of vector v of size n                           //
// VectorSum     | Compute the sum of the vector elements                                         //
// VectorMax     | Compute the max value of vector v of size n                                    //
// VectorMin     | Compute the min value of vector v of size n                                    //
// VectorTest    | Perform a boolean test on each vector element                                  //
// VectorRange   | Compute the range of a vector                                                  //
// VectorVariab  | Compute vector variability measure                                             //
// VectorEntropy | Compute vector entropy                                                         //
// VectorNonZero | Checks if vector <> 0                                                          //
// VectorCount   | Count nonempty values                                                          //
// VectorL2      | L2 distance                                                                    //
// Apply         | Apply a vector function to all rows of a matrix                                //
// Hist          | Histogram                                                                      //
// GetThreshold  | Determine a threshold for an exponentially decreasing distribution             //
//------------------------------------------------------------------------------------------------//

//------VectorAvg--------
//
template <typename type> type VectorAvg(type *v, int n)
{
  type m = 0;
  int C = 0;
  for(int i=0; i<n; i++) 
    if (v[i]==v[i]) {
      m += v[i];
      C++;
    }
  return m/C;
}
// explicit function template instantiations
template int VectorAvg<int>(int *, int); 
template long int VectorAvg<long int>(long int *, int); 
template long long int VectorAvg<long long int>(long long int *, int); 
template float VectorAvg<float>(float *, int);
template double VectorAvg<double>(double *, int);


//------VectorStd--------
//
template <typename type> type VectorStd(type *v, int n)
{
  type m = VectorAvg(v,n);
  type s = 0;
  int C = 0;
  for(int i=0; i<n; i++) 
    if (v[i]==v[i]) {
      s += v[i]*v[i];
      C++;
    }
  return (type)sqrt((double)s/C-m*m);
}
// explicit function template instantiations
template int VectorStd<int>(int *, int); 
template long int VectorStd<long int>(long int *, int); 
template long long int VectorStd<long long int>(long long int *, int); 
template float VectorStd<float>(float *, int); 
template double VectorStd<double>(double *, int); 




//------VectorSum-------------
//
template <typename type> type VectorSum(type *v, int n)
{
  type sum = 0;
  for(int i=0; i<n; i++) if (v[i]==v[i]) sum += v[i];
  return sum;
}
// explicit function template instantiations
template int VectorSum<int>(int*, int);  
template long int VectorSum<long int>(long int *, int);  
template long long int VectorSum<long long int>(long long int *, int);  
template unsigned long int VectorSum<unsigned long int>(unsigned long int *, int);  
template float VectorSum<float>(float *, int); 
template double VectorSum<double>(double *, int); 



//------VectorSumSq-------------
//
template <typename type> type VectorSumSq(type *v, int n)
{
  type sum = 0;
  for(int i=0; i<n; i++) if (v[i]==v[i]) sum += v[i]*v[i];
  return sum;
}
// explicit function template instantiations
template int VectorSumSq<int>(int*, int);  
template long int VectorSumSq<long int>(long int *, int);  
template long long int VectorSumSq<long long int>(long long int *, int);  
template unsigned long int VectorSumSq<unsigned long int>(unsigned long int *, int);  
template float VectorSumSq<float>(float *, int); 
template double VectorSumSq<double>(double *, int); 



//------VectorMax-------------
//
template <typename type> type VectorMax(type *v, int n)
{
  type max = v[0];
  int i;
  for(i=0; i<n; i++) if (v[i]==v[i]) { max = v[i]; break; }
  for( ; i<n; i++) if ((v[i]==v[i])&&(v[i]>max)) max = v[i];
  return max;
}
// explicit function template instantiations
template int VectorMax<int>(int *, int);  
template long int VectorMax<long int>(long int *, int);  
template long long int VectorMax<long long int>(long long int *, int);  
template float VectorMax<float>(float *, int); 
template double VectorMax<double>(double *, int); 



//------VectorMin-------------
//
template <typename type> type VectorMin(type *v, int n)
{
  type min = v[0];
  int i;
  for(i=0; i<n; i++) if (v[i]==v[i]) { min = v[i]; break; }
  for( ; i<n; i++) if ((v[i]==v[i])&&(v[i]<min)) min = v[i];
  return min;
}
// explicit function template instantiations
template int VectorMin<int>(int *, int);  
template long int VectorMin<long int>(long int *, int);  
template long long int VectorMin<long long int>(long long int *, int);  
template float VectorMin<float>(float *, int); 
template double VectorMin<double>(double *, int); 



//------VectorTest-------------
//
template <class type> int VectorTest(type *v, int n, bool (*func)(type))
{
  int C = 0;
  for(int i=0; i<n; i++) C += func(v[i]);
  return C;
}
// explicit function template instantiations
template int VectorTest<int>  (int *, int, bool (*func)(int));  
template int VectorTest<float>(float *, int, bool (*func)(float));  



//------VectorRange-----------
//
template <typename type> type VectorRange(type *v, int n)
{
  return VectorMax(v,n)-VectorMin(v,n);
}
// explicit function template instantiations
template int   VectorRange<int>  (int *, int);  
template float VectorRange<float>(float *, int); 



//------VectorNonZero---------
//
template <typename type> bool VectorNonZero(type *X, int n)
{
  for (int i=0; i<n; i++) if (X[i]!=0) return true;
  return false;
}
// explicit function template instantiations
template bool VectorNonZero<int>  (int *, int);  
template bool VectorNonZero<float>(float *, int); 



//------VectorCount---------
//
template <typename type> int VectorCount(type *X, int n)
{
  int nonempty = 0;
  for (int i=0; i<n; i++) nonempty += X[i]==X[i];
  return nonempty;
}
// explicit function template instantiations
template int VectorCount<int>  (int *, int);  
template int VectorCount<float>(float *, int); 



//------VectorL2---------
//
template <typename type> float VectorL2(type *X, int n)
{
  float y = 0;
  for (int i=0; i<n; i++) if (!IsNaN(X[i])) y += X[i]*X[i];
  return sqrt(y);
}
// explicit function template instantiations
template float VectorL2<int>  (int *, int);  
template float VectorL2<float>(float *, int); 



//----VectorVariab---------
//
float VectorVariab(float *X, int n)
{
  float *Y;
  ALLOCATE1D(Y,n,float);
  float avg = VectorAvg(X,n);
  float min = VectorMin(X,n);
  float max = VectorMax(X,n);
  for (int i=0; i<n; i++) Y[i] = pow((double)(X[i]-avg),2.0)/(max-min);
  VectorNormRange(Y,n);
  VectorSort(Y,n);
  int nbins = 20;
  float *b = Hist(Y,n,nbins);
  float y = 0;
  for (int i=0; i<4; i++) y += b[i];
  return (n-y)/n;
}



//----VectorEntropy------
//
float VectorEntropy(float *X, int n, int nbins)
{
  float *b = Hist(X,n,nbins);
  VectorNormSum(b,nbins);
  float y = 0;
  for (int i=0; i<nbins; i++) if (b[i]!=0) y -= b[i]*log(b[i]);
  return y/log((double)nbins);
}



//-----Hist--------------
//
float *Hist(float *X, int N, int n_bins, float min, float max)
{
  float *bins;
  ALLOCATE1D_INIT(bins,n_bins,float,0);
  if (max<=min) { min = VectorMin(X,N); max = VectorMax(X,N); }
  
  for (int k=0; k<N; k++) {
    float val = (X[k]-min) / (max-min);
    int b = (int)(n_bins*val);
    bins[b<n_bins?b:n_bins-1]++;
  }
  
  return bins;
}



//--------GetThreshold-----------
//
float GetThreshold(float *X, int n, int delta)
{
  // sort and compute derivative and threshold
  int *I = VectorRank(X,n);
  float *sorted;
  ALLOCATE1D(sorted,n,float);
  for (int g=0; g<n; g++) sorted[g] = X[I[n-g-1]];
  float *smooth = VectorSmooth(sorted,n,delta);
  float *diff;
  ALLOCATE1D(diff,n-delta,float);
  for (int g=0; g<n-delta; g++) diff[g] = 1 - smooth[g]/smooth[g+1];
  float *smooth_diff = VectorSmooth(diff,n,delta);
  //for (int g=0; g<n-2*delta; g++) printf("%f\n", smooth_diff[g]);
  
  // determine the derivative cutoff
  float avg = VectorAvg(smooth_diff,n-2*delta); 
  //float max = VectorMax(smooth_diff,n-2*delta); 
  float cutoff = avg; //max + 2*(avg-max);
  //fprintf(stderr, "avg = %f; max = %f; cutoff = %f\n", avg, max, cutoff);

  // determine the variance cutoff
  for (int g=0; g<n-2*delta; g++) if (smooth_diff[g]>cutoff) { cutoff = smooth[g+delta]; break; }
  
  // clean up
  FREE1D(smooth_diff);
  FREE1D(smooth);
  FREE1D(diff);

  return cutoff;
}











//------------------------------------------------------------------------------------------------//
// VECTOR OPERATIONS                                                                              //
//------------------------------------------------------------------------------------------------//
// VectorSmooth    | Return a "smoothed" version of vector x                                      //
// VectorSort      | Sort the vector                                                              //
// VectorRank      | Get the vector element ranking information                                   //
// VectorChop      | Change values below a cutoff to NaN                                          //
// VectorNorm      | Normalize vector mean to zero and standard deviation to one                  //
// VectorNormMean  | Normalize vector mean to zero                                                //
// VectorNormRange | Normalize the range of a vector to 1                                         //
// VectorNormSum   | Normalize the sum to 1                                                       //
// VectorNormSum2  | Normalize the sum of squared values to 1                                     //
// VectorCut       | Return a 'cut' of the vector                                                 //
// VectorPrint     | Print vector on screen                                                       //
// Apply           | Apply a vector function to all rows of a matrix                              //
// Pairwise        | Apply a function to all pairs of vectors in a matrix                         //
// HistPairwise    | Apply a function to all pairs of vectors in a matrix & compute distribution  //
// Int2Float       | Convert a vector of integer to a vector of float                             //
//------------------------------------------------------------------------------------------------//


//_____VectorSmooth_____________
float *VectorSmooth(float *x, int n, int d)
{
  float *y;
  ALLOCATE1D(y,n-d+1,float);
  y[0] = x[0];
  for (int k=1; k<d; k++) y[0] += x[k];
  for (int k=1; k<n-d+1; k++) y[k] = y[k-1]-x[k-1]+x[k+d-1];
  for (int k=0; k<n-d+1; k++) y[k] /= d;
  return y;
}


//___VectorChop___________
void VectorChop(float *x, int n, float cutoff)
{
  for (int k=0; k<n; k++) if (Abs_(x[k])<cutoff) x[k] = numeric_limits<float>::quiet_NaN();     // 0.0/0.0;
}


//----------VectorCopy--------------
//
template <class type> type* VectorCopy(type *x, unsigned long int n)
{
  type* v = new type[n];
  for (unsigned long int k=0; k<n; k++) v[k] = x[k];
  return v;
}
template int* VectorCopy<int>(int *,unsigned long int);
template long int* VectorCopy<long int>(long int *,unsigned long int);
template long long int* VectorCopy<long long int>(long long int *,unsigned long int);
template float* VectorCopy<float>(float *,unsigned long int);
template double* VectorCopy<double>(double *,unsigned long int);



//----------VectorNorm--------------
//
template <class type> void VectorNorm(type *x, unsigned long int n)
{
  type avg = VectorAvg(x,n);
  type std = VectorStd(x,n);
  for (unsigned long int k=0; k<n; k++) x[k] = (x[k]-avg)/std;
}
template void VectorNorm<int>(int *,unsigned long int);
template void VectorNorm<long int>(long int *,unsigned long int);
template void VectorNorm<long long int>(long long int *,unsigned long int);
template void VectorNorm<float>(float *,unsigned long int);
template void VectorNorm<double>(double *,unsigned long int);



//----------VectorNormSum--------------
//
template <class type> void VectorNormSum(type *x, unsigned long int n)
{
  type sum = 0;
  for (unsigned long int k=0; k<n; k++) sum += x[k];
  for (unsigned long int k=0; k<n; k++) x[k] /= sum;
}
template void VectorNormSum<float>(float *,unsigned long int);
template void VectorNormSum<double>(double *,unsigned long int);



//----------VectorNormRange--------------
//
template <class type> void VectorNormRange(type *x, int n)
{
  type min = VectorMin(x,n);
  type max = VectorMax(x,n);
  for (int k=0; k<n; k++) x[k] = (x[k]-min)/(max-min);
}
template void VectorNormRange<int>(int *,int);
template void VectorNormRange<long int>(long int *,int);
template void VectorNormRange<long long int>(long long int *,int);
template void VectorNormRange<float>(float *,int);
template void VectorNormRange<double>(double *,int);



//___VectorNormSum2___________
void VectorNormSum2(float *x, int n)
{
  float sum = 0;
  for (int k=0; k<n; k++) sum += x[k]*x[k];
  for (int k=0; k<n; k++) x[k] /= sqrt(sum);
}


//___VectorNormMean___________
void VectorNormMean(float *x, int n)
{
  float avg = VectorAvg(x,n);
  for (int k=0; k<n; k++) x[k] = x[k]-avg;
}


//-----VectorSort------------
//
int compareSortF(const void *a, const void *b) { return *(float*)a > *(float*)b ? 1 : -1; }
int compareSortI(const void *a, const void *b) { return *(int*)a > *(int*)b ? 1 : -1; }
int compareSortLI(const void *a, const void *b) { return *(long int*)a > *(long int*)b ? 1 : -1; }
int compareSortULI(const void *a, const void *b) { return *(unsigned long int*)a > *(unsigned long int*)b ? 1 : -1; }
int compareSortD(const void *a, const void *b) { return *(double*)a > *(double*)b ? 1 : -1; }
void VectorSort(float *V, int n) { qsort(V,n,sizeof(float),compareSortF); }
void VectorSort(int *V, int n) { qsort(V,n,sizeof(int),compareSortI); }
void VectorSort(long int *V, int n) { qsort(V,n,sizeof(long int),compareSortLI); }
void VectorSort(unsigned long int *V, int n) { qsort(V,n,sizeof(unsigned long int),compareSortULI); }
void VectorSort(double *V, long int n) { qsort(V,n,sizeof(double),compareSortD); }



//-----VectorRank-----------
//
struct SortPair { int rank; float val; };
int compareRank(const void *a, const void *b) 
{ 
  struct SortPair *A = (SortPair *) a;
  struct SortPair *B = (SortPair *) b;
  return A->val>B->val?1:-1;
}
int *VectorRank(float *V, int n)
{
  SortPair *X;
  ALLOCATE1D(X,n,SortPair);
  for (int i=0; i<n; i++) {
    X[i].rank = i;
    X[i].val = V[i];
  }
  qsort(X,n,sizeof(SortPair),compareRank);
  
  int *I;
  ALLOCATE1D(I,n,int);
  for (int i=0; i<n; i++) I[i] = X[i].rank;
  delete [] X; 

  return I;
}



//-----VectorRank-----------
//
struct SortPairDouble { int rank; double val; };
int compareRankDouble(const void *a, const void *b) 
{ 
  struct SortPairDouble *A = (SortPairDouble *) a;
  struct SortPairDouble *B = (SortPairDouble *) b;
  return A->val>B->val?1:-1;
}
int *VectorRank(double *V, int n)
{
  SortPairDouble *X;
  ALLOCATE1D(X,n,SortPairDouble);
  for (int i=0; i<n; i++) {
    X[i].rank = i;
    X[i].val = V[i];
  }
  qsort(X,n,sizeof(SortPairDouble),compareRankDouble);
  
  int *I;
  ALLOCATE1D(I,n,int);
  for (int i=0; i<n; i++) I[i] = X[i].rank;
  delete [] X; 

  return I;
}

int *GetVectorRanks(double *V, int n)
{
  int *order = VectorRank(V,n);
  int *rank = new int[n];
  for (int i=0; i<n; i++) rank[order[i]] = i;
  delete order; 
  return rank;
}


//-----VectorRank-----------
//
struct SortPairLI { int rank; long int val; };
int compareRankLI(const void *a, const void *b) 
{ 
  struct SortPairLI *A = (SortPairLI *) a;
  struct SortPairLI *B = (SortPairLI *) b;
  return A->val>B->val?1:-1;
}
int *VectorRank(long int *V, int n)
{
  SortPairLI *X;
  ALLOCATE1D(X,n,SortPairLI);
  for (int i=0; i<n; i++) {
    X[i].rank = i;
    X[i].val = V[i];
  }
  qsort(X,n,sizeof(SortPairLI),compareRankLI);
  
  int *I;
  ALLOCATE1D(I,n,int);
  for (int i=0; i<n; i++) I[i] = X[i].rank;
  delete [] X; 

  return I;
}



//-----VectorCut----------
//
float *VectorCut(float *X, int n, bool *vert_cut)
{
  int nn = 0;
  if (vert_cut) for (int i=0; i<n; i++) nn += vert_cut[i];
  else nn = n;

  float *XX;
  ALLOCATE1D(XX,nn,float);
  
  for (int i=0,k=0; i<n; i++) 
    if ((vert_cut==NULL) || vert_cut[i]) {
      XX[k] = X[i];
      ++k;
    }

  return XX;
}



//-----VectorCut----------
//
float *VectorCut(float *X, int n, long int mask, int *nn)
{
  *nn = 0;
  long int column = 1;
  for (int c=0; c<n; c++,column<<=1) if (mask&column) (*nn)++;

  float *XX;
  ALLOCATE1D(XX,*nn,float);
  
  column = 1;
  for (int i=0,k=0; i<n; i++,column<<=1) if (mask&column) XX[k++] = X[i];

  return XX;
}



//-----VectorPrint----------
//
template <typename type> void VectorPrint(type *X, int n, const char *format, bool *vert_cut, FILE *output)
{
  for (int i=0; i<n; i++) if ((vert_cut==NULL) || vert_cut[i]) {
    if (X[i]==X[i]) fprintf(output, format, X[i]);
    else fprintf(output, "NaN");
    printf(" ");
  }
}

template void VectorPrint<int>  (int *, int, const char *, bool *, FILE *);  
template void VectorPrint<float>(float *, int, const char *, bool *, FILE *);  
template void VectorPrint<double>(double *, int, const char *, bool *, FILE *);  



//------Apply----------
//
float *Apply(float **A, int N, int M, float (*func)(float *,int))
{
  float *X;
  ALLOCATE1D(X,N,float);
  for (int k=0; k<N; k++) X[k] = (*func)(A[k], M);
  return X;
}



//------Pairwise------------
//
float **Pairwise(float **A, int N, int M, float (*func)(float *,float *,unsigned long int))
{
  float **pairs;
  ALLOCATE2D(pairs,N,N,float);

  Progress PRG("Computing pairwise distance matrix...",N*(N+1)/2);
  for (int a=0; a<N; a++) {
    for(int b=a; b<N; b++) {
      pairs[b][a] = pairs[a][b] = (*func)(A[a], A[b], M);
      PRG.Check();
    }
  }
  PRG.Done();
  return pairs;
}



//------Pairwise------------
//
float **Pairwise(float **X, int Nx, float **Y, int Ny, int M, float (*func)(float *,float *,unsigned long int))
{
  float **pairs;
  ALLOCATE2D(pairs,Nx,Ny,float);

  Progress PRG("Computing pairwise distance matrix...",Nx*Ny);
  for (int i=0; i<Nx; i++) {
    for(int j=0; j<Ny; j++) {
      pairs[i][j] = (*func)(X[i],Y[j],M);
      PRG.Check();
    }
  }
  PRG.Done();
  return pairs;
}




//-----HistPairwise-------------
//
float *HistPairwise(float **A, int N, int M, float (*func)(float *,float *,unsigned long int), int n_bins)
{
  float *bins;
  ALLOCATE1D_INIT(bins,n_bins,float,0);

  Progress PROGR("<HistPairwise>: computing pairwise correlations...",N*(N-1)/2);

  for (int i=0; i<N; i++) {
    for(int j=i+1; j<N; j++) {
      float val = (1 + (*func)(A[i], A[j], M)) / 2;
      int b = (int)(n_bins*val);
      bins[b<n_bins?b:n_bins]++;
      PROGR.Check();
    }
  }
  PROGR.Done();

  return bins;
}



//-----HistPairwise-------------
//
float *HistPairwise(float **A, float **B, int Na, int Nb, int M, float (*func)(float *,float *,unsigned long int), int n_bins)
{
  float *bins;
  ALLOCATE1D_INIT(bins,n_bins,float,0);

  Progress PROGR("<HistPairwise>: computing pairwise correlations...", Na*Nb);

  for (int i=0; i<Na; i++) {
    for(int j=0; j<Nb; j++) {
      float val = (1 + (*func)(A[i], B[j], M)) / 2;
      int b = (int)(n_bins*val);
      bins[b<n_bins?b:n_bins]++;
      PROGR.Check();
    }
  }
  PROGR.Done();

  return bins;
}



//-------Int2Float-------
//
float *Int2Float(int *x, int n)
{
  float *y;
  ALLOCATE1D(y,n,float);
  for (int i=0; i<n; i++) y[i] = x[i];
  return y;
}












//------------------------------------------------------------------------------------------------//
// VECTOR DISTANCES                                                                               //
//------------------------------------------------------------------------------------------------//
// InnerProduct    | The inner product of two vectors                                             //
// Covariance      | The inner product of two vectors                                             //
// Euclidean       | The euclidean distance                                                       //
// VectorCorr      | Calculate the correlation of two vectors                                     //
// RelativeEntropy | Calculate the relative entropy of distribution p with respect to q           //
// Chi2            | Calculate the chi-square distance for vector x relative to expected values   //
// Chi2std         | Calculate the chi-square distance for vector x relative to mean and stdev    //
// Corr            | Calculate the correlation two vectors (handles NaN values)                   //
// Hamming1        | Calculate the common 'ones' in two vectors                                   //
// Jaccard         | Calculate the Jaccard of two vectors                                         //
//------------------------------------------------------------------------------------------------//


//------InnerProduct------
//
template <class type> type InnerProduct(type *A, type *B, unsigned long int n)
{
  type C = 0;
  for (unsigned long int i=0; i<n; i++) C += A[i]*B[i];
  return C; 
}
template float InnerProduct<float>(float *, float *, unsigned long int);
template double InnerProduct<double>(double *, double *, unsigned long int);



//------Covariance------
//
template <class type> type Covariance(type *A, type *B, unsigned long int n)
{
  return InnerProduct(A,B,n)/n;
}
template float Covariance<float>(float *, float *, unsigned long int);
template double Covariance<double>(double *, double *, unsigned long int);



//------Euclidean------
//
template <class type> type Euclidean(type *A, type *B, unsigned long int n)
{
  type d = 0;
  for (unsigned long int k=0; k<n; k++) d += pow((double)(A[k]-B[k]),2.0);
  return sqrt(d);
}
template float Euclidean<float>(float *, float *, unsigned long int);
template double Euclidean<double>(double *, double *, unsigned long int);



//-------VectorCorr-------
//
template <class type> double VectorCorr(type *A, type *B, unsigned long int n)
{
  double Ex = 0, Ey = 0;
  double Ex2 = 0, Ey2 = 0;
  double Exy = 0;
  unsigned long int C = 0;
  for (unsigned long int i=0; i<n; i++) 
    if ((A[i]==A[i])&&(B[i]==B[i])) {
      C++;
      Ex += A[i];
      Ex2 += pow((double)A[i],2.0);
      Ey += B[i];
      Ey2 += pow((double)B[i],2.0);
      Exy += A[i]*B[i];
    }
  Ex = Ex/C;
  Ey = Ey/C;
  Ex2 = Ex2/C;
  Ey2 = Ey2/C;
  Exy = Exy/C;
  double Y = (Exy-Ex*Ey)/sqrt((Ex2-pow(Ex,2.0))*(Ey2-pow(Ey,2.0)));
  if (fabs(Y)>1.00001) fprintf(stderr, "Correlation error: (%f,%f,%f)\n", Y, Ex2-pow(Ex,2.0), Ey2-pow(Ey,2.0));
  return Y; 
}
template double VectorCorr<float>(float *, float *, unsigned long int);
template double VectorCorr<double>(double *, double *, unsigned long int);
template double VectorCorr<int>(int *, int *, unsigned long int);
template double VectorCorr<long int>(long int *, long int *, unsigned long int);
template double VectorCorr<long long int>(long long int *, long long int *, unsigned long int);





//------RelativeEntropy------
//
template <class type> type RelativeEntropy(type *p, type *q, unsigned long int n)
{
  type r = 0;
  type N = 0, M = 0;
  for (unsigned long int k=0; k<n; k++) { N += p[k]; M += q[k]; }
  for (unsigned long int k=0; k<n; k++) if (p[k]>0) r += p[k]/N*log(p[k]/N/(q[k]/M));
  return r/log(2.0);
}
template float RelativeEntropy<float>(float *, float *, unsigned long int);
template double RelativeEntropy<double>(double *, double *, unsigned long int);



//------Chi2------
//
template <class type> type Chi2(type *x, type *avg, unsigned long int n)
{
  type d = 0;
  for (unsigned long int k=0; k<n; k++) d += pow((double)(x[k]-avg[k]),(double)2)/avg[k];
//fprintf(stderr,"x=%f m=%f n=%i => y=%f",x[0],avg[0],n,sqrt(d));getchar();
  return sqrt(d);
}
template float Chi2<float>(float *, float *, unsigned long int);
template double Chi2<double>(double *, double *, unsigned long int);



//------Chi2std------
//
template <class type> type Chi2std(type *x, type *avg, type *std, unsigned long int n)
{
  type d = 0;
  for (unsigned long int k=0; k<n; k++) d += pow((double)(x[k]-avg[k])/std[k],(double)2);
  return sqrt(d);
}
template float Chi2std<float>(float *, float *, float *, unsigned long int);
template double Chi2std<double>(double *, double *, double *, unsigned long int);





//------Corr------
//
float Corr(float *A, float *B, unsigned long int M)
{
  double Ex = 0, Ey = 0;
  double Ex2 = 0, Ey2 = 0;
  double Exy = 0;
  int C = 0; 
  for (unsigned long int i=0; i<M; i++) 
    if ((A[i]==A[i])&&(B[i]==B[i])) {
      Ex += A[i];
      Ex2 += pow((double)A[i],2.0);
      Ey += B[i];
      Ey2 += pow((double)B[i],2.0);
      Exy += A[i]*B[i];
      C++;
    }
  Ex = Ex/C;
  Ey = Ey/C;
  Ex2 = Ex2/C;
  Ey2 = Ey2/C;
  Exy = Exy/C;
  float y = (Exy-Ex*Ey)/sqrt((Ex2-pow(Ex,2.0))*(Ey2-pow(Ey,2.0)));
  if ((y>1)||(C<3)) return 0;
  return y;
}



//------Corr------
//
float Corr(float *A, float *B, unsigned long int *I, unsigned long int M)
{
  double Ex = 0, Ey = 0;
  double Ex2 = 0, Ey2 = 0;
  double Exy = 0;
  int C = 0; 
  for (unsigned long int k=0; k<M; k++)
  {
    unsigned long int i = I[k];
    if ((A[i]==A[i])&&(B[i]==B[i])) {
      Ex += A[i];
      Ex2 += pow((double)A[i],2.0);
      Ey += B[i];
      Ey2 += pow((double)B[i],2.0);
      Exy += A[i]*B[i];
      C++;
    }
  }
  Ex = Ex/C;
  Ey = Ey/C;
  Ex2 = Ex2/C;
  Ey2 = Ey2/C;
  Exy = Exy/C;
  float y = (Exy-Ex*Ey)/sqrt((Ex2-pow(Ex,2.0))*(Ey2-pow(Ey,2.0)));
  return y;
}


//------Hamming1------
//
int Hamming1(int *x, int *y, int n)
{
  int c = 0;  // common number of 'ones'
  for (int i=0; i<n; i++) c += (x[i]==1)&&(y[i]==1);
  return c;
}


//------Jaccard------
//
float Jaccard(int *A, int *B, int n)
{
  float d = 0;
  int C = 0;
  for (int i=0; i<n; i++) { d += (A[i]>0)&&(B[i]>0); C += (A[i]>0)||(B[i]>0); }
  return 1.0-d/C; 
}








//------------------------------------------------------------------------------------------------//
// MISCELLANEOUS                                                                                  //
//------------------------------------------------------------------------------------------------//
// Abs_                | Absolute value                                                           //
// Min_                | Minimum of two numbers                                                   //
// Max_                | Maximum of two numbers                                                   //
// IsNaN               | Check for NaN values                                                     //
// InitRandomGenerator | Initialize a random generator                                            //
//------------------------------------------------------------------------------------------------//



float Abs_(float x) { return x>0?x:-x; }
double Abs_(double x) { return x>0?x:-x; }
bool IsNaN(float x) { return x!=x; }
int Min_(int x, int y) { return x<y?x:y; }
int Max_(int x, int y) { return x>y?x:y; }
long int Min_(long int x, long int y) { return x<y?x:y; }
long int Max_(long int x, long int y) { return x>y?x:y; }
size_t Min_(size_t x, size_t y) { return x<y?x:y; }
size_t Max_(size_t x, size_t y) { return x>y?x:y; }










//------------------------------------------------------------------------------------------------//
// VECTOR I/O                                                                                     //
//------------------------------------------------------------------------------------------------//
// Exists      | Check whether a file exists                                                      //
// LoadStdIn   | Load standard input into a buffer                                                //
// LoadFile    | Load a file into a buffer                                                        //
// LoadVectors | Load an array of vectors from a file into memory                                 //
// SaveVectors | Save an array of vectors from memory into a file                                 //
// ReadVec     | Process a single vector from an input string                                     //
// LoadMatrix  | load a matrix from a file with error checking                                    //
//------------------------------------------------------------------------------------------------//

//-----Exists-----
//
bool Exists(char *file)
{
  FILE *F = fopen(file,"r");
  if (F!=NULL) {
    fclose(F);
    return true;
  }
  return false;
}


//-----GetFileType---------
//
int GetFileType(const char *file)
{
  FILE *F = fopen(file,"r");
  if (F==NULL) return -1;
  int byte1 = fgetc(F);
  int byte2 = fgetc(F);  
  fclose(F);
  if ((byte1==0x1f)&&(byte2==0x8b)) {
    igzstream file_stream(file);
	char c1 = file_stream.get();
	char c2 = file_stream.get();
	char c3 = file_stream.get();
	char c4 = file_stream.get();
	file_stream.close();
    if ((c1=='B')&&(c2=='A')&&(c3=='M')&&(c4=='\1')) return 2;   // bam format
	else return 1;   // gz format
  }
  else return 0;    // text format
}



//---LoadStdIn---
//
char *LoadStdIn() 
{
  FILE *F = tmpfile();
  char *B;
  ALLOCATE1D(B,MAX_BUFFER_SIZE+1,char);
  for (int z=1; ; z++) {
    if (fgets(B,MAX_BUFFER_SIZE,stdin)==NULL) break;             // quit when end of file is reached
    int blen = strlen(B);
    if (blen>=MAX_BUFFER_SIZE-1) { fprintf(stderr, "<LoadStdIn>: line not entirely read!\n"); exit(1); }  
    fwrite(B,1,strlen(B),F);
    //if (z%100000==0) fprintf(stderr, "%i\n", z);
  }
  rewind(F);
  char *buffer = LoadFile(F);
  fclose(F);
  FREE1D(B);
  return buffer;
}



//-----LoadStdIn-----
//
FILE *LoadStdIn(long int *n_lines, long int buffer_size)
{
  FILE *F = tmpfile();
  char *B;
  ALLOCATE1D(B,buffer_size+1,char);
  for (*n_lines=0; ; (*n_lines)++) {
    if (fgets(B,buffer_size,stdin)==NULL) break;             // quit when end of file is reached
    int blen = strlen(B);
    if (blen>=buffer_size-1) { fprintf(stderr, "Line %ld: line not entirely read!\n", *n_lines+1); exit(1); }  
    fwrite(B,1,blen,F);
  }
  rewind(F);
  return F;
}



//------LoadFile-----
//
char *LoadFile(FILE *F) 
{
  fseek(F,0,SEEK_END);
  long int buffer_size = ftell(F);
  //if (_MESSAGES_) fprintf(stderr, "* Read %ld bytes.\n", buffer_size);
  if (buffer_size<0) { fprintf(stderr, "<LoadFile>: error reading file!\n"); exit(1); }
  rewind(F);
  char *buffer;
  ALLOCATE1D(buffer,buffer_size+1,char);
  fread(buffer,1,buffer_size,F);
  buffer[buffer_size] = 0;
  return buffer;
}



//------LoadFile-----
//
char *LoadFile(char *file) 
{
  FILE *F = fopen(file,"r");
  if (F==0) { fprintf(stderr,"<LoadFile>: can't open file '%s'!\n", file); exit(1); }
  char *buffer = LoadFile(F);
  fclose(F);
  return buffer;
}



//-----LoadVectors-----
//
float **LoadVectors(char *file, int *N, int *M)
{
  char *buffer = file==NULL ? LoadStdIn() : LoadFile(file);
  Fields fields(buffer," \n");
  float **V = fields.GetFloatMatrix(N,M);
  free(buffer);
  return V;
}



//------SaveVectors-----
//
void SaveVectors(char *file, char *format, float **data, int *FILTER, int N, int M)
{
  FILE *F;
  fprintf(stderr, "Saving data to file '%s'...\n", file);
  if ((F = fopen(file, "w")) == 0) { printf("File '%s' cannot be created!\n", file); return; };
  Progress PRG("Storing data...",N*M);
  for (int n=0; n<N; n++) if ((FILTER==NULL)||(FILTER[n]==1)) {
    for (int m=0; m<M; m++) { 
      fprintf(F, format, data[n][m]);
      PRG.Check();
    }    
    fprintf(F, "\n");
  }
  PRG.Done();
  fclose(F);
}



//------SaveVectors-----
//
void SaveVectors(char *file, char *format, float **data, int N, int M)
{
  SaveVectors(file, format, data, NULL, N, M);
}



//------SaveVectors-----
//
void SaveVectors(char *file, int **data, int N, int M)
{
  FILE *F;
  if ((F = fopen(file, "w")) == 0) { printf("File '%s' cannot be created!\n", file); return; };
  for (int n=0; n<N; n++) {
    for (int m=0; m<M; m++) fprintf(F, "%i ", data[n][m]);
    fprintf(F, "\n");
  }
  fclose(F);
}




//-----LoadMatrix-----
//
float **LoadMatrix(char *file, long int *n_rows, long int *n_cols, char sep)
{
  // read input
  char *buffer = file==NULL ? LoadStdIn() : LoadFile(file);

  // allocate memory
  *n_rows = CountTokens(buffer,'\n');
  float **MATRIX;
  ALLOCATE1D(MATRIX,*n_rows,float *);
  
  // read vectors
  char *inp = buffer;
  Progress PRG("Reading vectors...",*n_rows);
  for (long int k=0; k<*n_rows; k++) {
    char *vec = GetNextToken(&inp,'\n');
    long int n;
    MATRIX[k] = ReadFloatVec(vec,&n,sep);
    if (k==0) *n_cols = n; 
    else if (*n_cols!=n) { fprintf(stderr, "Line %ld: number of columns (%ld) should be equal to %ld!\n%s\n", k+1, n, *n_cols, vec); exit(1); } 
    PRG.Check();
  } 
  PRG.Done();


  // cleanup
  free(buffer);

  return MATRIX;
}




//-----LoadDoubleMatrix-----
//
double **LoadDoubleMatrix(char *file, long int *n_rows, long int *n_cols)
{
  // read input
  char *buffer = file==NULL ? LoadStdIn() : LoadFile(file);

  // allocate memory
  *n_rows = CountTokens(buffer,'\n');
  double **MATRIX;
  ALLOCATE1D(MATRIX,*n_rows,double *);
  
  // read vectors
  char *inp = buffer;
  Progress PRG("Reading vectors...",*n_rows);
  for (long int k=0; k<*n_rows; k++) {
    char *vec = GetNextToken(&inp,'\n');
    long int n;
    MATRIX[k] = ReadDoubleVec(vec,&n);
    if (k==0) *n_cols = n; 
    else if (*n_cols!=n) { fprintf(stderr, "Line %ld: number of columns (%ld) should be equal to %ld!\n%s\n", k+1, n, *n_cols, vec); exit(1); } 
    PRG.Check();
  } 
  PRG.Done();

  // cleanup
  free(buffer);

  return MATRIX;
}




//-----LoadIntMatrix-----
//
int **LoadIntMatrix(char *file, long int *n_rows, long int *n_cols)
{
  // read input
  char *buffer = file==NULL ? LoadStdIn() : LoadFile(file);

  // allocate memory
  *n_rows = CountTokens(buffer,'\n');
  int **MATRIX;
  ALLOCATE1D(MATRIX,*n_rows,int *);
  
  // read vectors
  char *inp = buffer;
  Progress PRG("Reading vectors...",*n_rows);
  for (long int k=0; k<*n_rows; k++) {
    char *vec = GetNextToken(&inp,'\n');
    long int n;
    MATRIX[k] = ReadIntVec(vec,&n);
    if (k==0) *n_cols = n; 
    else if (*n_cols!=n) { fprintf(stderr, "Line %ld: number of columns (%ld) should be equal to %ld!\n%s\n", k+1, n, *n_cols, vec); exit(1); } 
    PRG.Check();
  } 
  PRG.Done();


  // cleanup
  free(buffer);

  return MATRIX;
}







//------------------------------------------------------------------------------------------------//
// STRINGS                                                                                        //
//------------------------------------------------------------------------------------------------//
// StrCopy        | Make a new copy of a string                                                   //
// StrChop        | Removes trailing newlines                                                     //
// AddSuffix      | Add a suffix to a string                                                      //
// UpperCase      | Convert to uppercase                                                          //
// LowerCase      | Convert to lowercase                                                          //
// IsSuffix       | Test if suffix matches                                                        //
// StrCountLines  | Count newline characters                                                      //
//------------------------------------------------------------------------------------------------//

//------StrCopy--------
//
char *StrCopy(const char *source)
{
  if (source==NULL) return NULL;
  char *target;
  ALLOCATE1D(target,strlen(source)+1,char);
  strcpy(target,source);
  return target;
}



//------StrChop--------
//
void StrChop(char *source, char c)
{
  if (source==NULL) return;
  char *p = source + strlen(source) - 1;
  while (p!=source) if (p[0]==c) p--; else break;
  p[1] = 0;
}



//------AddSuffix--------
//
char *AddSuffix(char *source, char *suffix)
{
  char *out;
  ALLOCATE1D(out,strlen(source)+strlen(suffix)+1,char);
  sprintf(out, "%s%s", source, suffix);
  //fprintf(stderr, "<AddSuffix>: %s%s = %s\n", source, suffix, out);
  return out;
}



//-----UpperCase-----
//
void UpperCase(char *s)
{
  for (char *p=s; p[0]!=0; p++) p[0] = toupper(p[0]);
}



//-----LowerCase-----
//
void LowerCase(char *s)
{
  for (char *p=s; p[0]!=0; p++) p[0] = tolower(p[0]);
}



//-----IsSuffix-----
//
bool IsSuffix(char *str, char *sfx)
{
  if (strlen(str)<strlen(sfx)) return false;
  for (int i=strlen(str)-1,j=strlen(sfx)-1; j>=0; i--,j--) if (str[i]!=sfx[j]) return false;
  return true;
}



//-----StrCountLines-----
//
long int StrCountLines(char *s)
{
  long int n = 0;
  for (; s[0]; s++) if (s[0]=='\n') n++;
  return n;
}




//-----MyRand-----------
//
unsigned long int  MyRand(unsigned long int N)
{
  //double s = (RAND_MAX/N)/((double)RAND_MAX/N);
  //return ((unsigned long int)floor(s*rand()))%N;
  return rand()%N;
}












//-----complement-----------
//
char complement(char c)
{
  switch (c) {
  case 'a': return 't';
  case 'c': return 'g';
  case 'g': return 'c';
  case 't': return 'a';
  case 'r': return 'y';
  case 'y': return 'r';
  case 'w': return 'w';
  case 's': return 's';
  case 'k': return 'm';
  case 'm': return 'k';
  case 'b': return 'v';
  case 'v': return 'b';
  case 'd': return 'h';
  case 'h': return 'd';
  case 'n': return 'n';

  case 'A': return 'T';
  case 'C': return 'G';
  case 'G': return 'C';
  case 'T': return 'A';
  case 'R': return 'Y';
  case 'Y': return 'R';
  case 'W': return 'W';
  case 'S': return 'S';
  case 'K': return 'M';
  case 'M': return 'K';
  case 'B': return 'V';
  case 'V': return 'B';
  case 'D': return 'H';
  case 'H': return 'D';
  case 'N': return 'N';
  }
  return c;
}



//-----ReverseComplement-----------
//
string ReverseComplement(string &x)
{
  size_t len = x.size();
  string y("0",len);
  size_t i, j;
  for (i=0,j=len-1; i<len; i++,j--) y[i] = complement(x[j]);
  return y;
}




//-----ReverseGene-----------
//
void ReverseGene(char *s, unsigned long int len)
{
  if (len==0) return;
  unsigned long int i,j;
  char t;
  for (i=0,j=len-1; i<j; i++,j--) {
    t = s[i];
    s[i] = complement(s[j]);
    s[j] = complement(t);
  }
  if (i==j) s[i] = complement(s[i]);
}





void CmdOptionPrint(bool val) { printf("[%s]", val?"true":"false"); }
void CmdOptionPrint(char val) { printf("[%c]", val); }
void CmdOptionPrint(int val) { printf("[%d]", val); }
void CmdOptionPrint(long int val) { printf("[%ld]", val); }
void CmdOptionPrint(unsigned long int val) { printf("[%lu]", val); }
void CmdOptionPrint(float val) { printf("[%.6f]", val); }
void CmdOptionPrint(double val) { printf("[%.6e]", val); }

void CmdOptionRead(bool *ptr, char *s) { *ptr = !*ptr; }
void CmdOptionRead(char *ptr, char *s) { *ptr = s[0]; }
void CmdOptionRead(int *ptr, char *s) { *ptr = atoi(s); }
void CmdOptionRead(long int *ptr, char *s) { *ptr = atol(s); }
void CmdOptionRead(unsigned long int *ptr, char *s) { *ptr = (unsigned long int)atol(s); }
void CmdOptionRead(float *ptr, char *s) { *ptr = atof(s); }
void CmdOptionRead(double *ptr, char *s) { *ptr = (double)atof(s); }




//---------------------------------------------------------------------------------------------------//
// CLASS CmdOption                                                                                   //
//---------------------------------------------------------------------------------------------------//


CmdOption::CmdOption(const char *opt, const char *description)
{
  this->opt = new char[strlen(opt)+1];
  strcpy(this->opt,opt);
  this->description = new char[strlen(description)+1];
  strcpy(this->description,description);
}

void CmdOption::Init() {}

void CmdOption::Print() {}

void CmdOption::Read(int *argc, char ***argv) {}

CmdOption::~CmdOption() {}





//---------------------------------------------------------------------------------------------------//
// CLASS CmdOptionTemplate                                                                           //
//---------------------------------------------------------------------------------------------------//

//--------Constructor-----------
//
template <class T> CmdOptionTemplate<T>::CmdOptionTemplate(const char *opt, T *ptr, T val, const char *description) : CmdOption(opt,description)
{
  this->ptr = ptr;
  this->val = val;
}

//--------Destructor-----------
//
template <class T> CmdOptionTemplate<T>::~CmdOptionTemplate()
{
}

//--------Init-----------
//
template <class T> void CmdOptionTemplate<T>::Init()
{
  *ptr = val;
}

//--------Print-----------
//
template <class T> void CmdOptionTemplate<T>::Print()
{
  CmdOptionPrint(*ptr);
}

//--------Read-----------
//
template <class T> void CmdOptionTemplate<T>::Read(int *argc, char ***argv)
{
  ++(*argv);
  --(*argc);
  if(*argc == 0) { fprintf(stderr, "Error: could not set option '%s'!\n", opt); exit(1); }
  CmdOptionRead(ptr,*argv[0]);
  (*argc)--;
  (*argv)++;
}



//---------------------------------------------------------------------------------------------------//
// CLASS BoolCmdOption                                                                               //
//---------------------------------------------------------------------------------------------------//

//--------Constructor-----------
//
BoolCmdOption::BoolCmdOption(const char *opt, bool *ptr, bool val, const char *description) : CmdOption(opt,description)
{
  this->ptr = ptr;
  this->val = val;
}

//--------Destructor-----------
//
BoolCmdOption::~BoolCmdOption()
{
}

//--------Init-----------
//
void BoolCmdOption::Init()
{
  *ptr = val;
}

//--------Print-----------
//
void BoolCmdOption::Print()
{
  CmdOptionPrint(*ptr);
}

//--------Read-----------
//
void BoolCmdOption::Read(int *argc, char ***argv)
{
  if(*argc == 0) { fprintf(stderr, "Error: could not set option '%s'!\n", opt); exit(1); }
  CmdOptionRead(ptr,*argv[0]);
  (*argc)--;
  (*argv)++;
}



//---------------------------------------------------------------------------------------------------//
// CLASS StrCmdOption                                                                                //
//---------------------------------------------------------------------------------------------------//


StrCmdOption::StrCmdOption(const char *opt, char **ptr, const char *val, const char *description) : CmdOption(opt,description)
{
  this->val = new char[strlen(val)+1];
  strcpy(this->val,val);
  this->ptr = ptr;
}

void StrCmdOption::Init()
{
  *ptr = val;
}

void StrCmdOption::Print()
{
  printf("[%s]", *ptr);
}


void StrCmdOption::Read(int *argc, char ***argv)
{
  (*argc)--;
  (*argv)++;
  if (*argc == 0) { fprintf(stderr, "Error: could not set option '%s'\n", opt); exit(1); }
  delete val;
  val = new char[strlen(*argv[0])+1];
  strcpy(val,*argv[0]);
  *ptr = val;
  (*argc)--;
  (*argv)++;
}


StrCmdOption::~StrCmdOption()
{
  if (val!=NULL) delete val;
}









//---------------------------------------------------------------------------------------------------//
// CLASS CmdLine                                                                                     //
//---------------------------------------------------------------------------------------------------//


//------Constructor-----------
//
CmdLine::CmdLine()
{

}


//------Destructor-----------
//
CmdLine::~CmdLine()
{
  for (OptionMap::iterator x=cmd_options.begin(); x!=cmd_options.end(); x++) delete x->second;
}


//------SetProgramName-----------
//
void CmdLine::SetProgramName(string program_name, string version)
{
  this->program_name = program_name;
  this->version = version;
}


//------Read----------------
//
int CmdLine::Read(char **argv, int argc)
{
  // initialize
  Init();

  // read options
  int _argc = argc-1;
  char **_argv = argv+1;
  while (_argc>0) {
    if (_argv[0][0]!='-') return argc - _argc;
    OptionMap::iterator x = cmd_options.find(_argv[0]);
    if (x!=cmd_options.end()) x->second->Read(&_argc,&_argv);
    else { fprintf(stderr, "Error: unknown option '%s'!\n", _argv[0]); exit(1); }
  }

  return argc;
}



//------Init---------------
//
void CmdLine::Init()
{
  for (OptionMap::iterator x=cmd_options.begin(); x!=cmd_options.end(); x++) x->second->Init();
}



//------Print---------------
//
void CmdLine::Print()
{
  for (OptionList::iterator x=cmd_option_list.begin(); x!=cmd_option_list.end(); x++) {
    printf("  %-25s %-80s ", (*x)->opt, (*x)->description);
    (*x)->Print();
    printf("\n");
  }
  /*for (OptionMap::iterator x=cmd_options.begin(); x!=cmd_options.end(); x++) {
    printf("%-10s %-80s ", x->second->opt, x->second->description);
    x->second->Print();
    printf("\n");
  }*/
}



//------Usage-----------
//
void CmdLine::Usage(const char *text)
{
  printf("\n");
  printf("USAGE: \n");
  printf("  %s\n", text);
  printf("OPTIONS: \n");
  Print();
  printf("\n");
}


//------Usage-----------
//
void CmdLine::Usage(char *prog_name, const char *usage)
{
  printf("\n");
  printf("USAGE: \n");
  printf("  %s %s\n", prog_name, usage);
  printf("OPTIONS: \n");
  Print();
  printf("\n");
}



//------AddOption-----------
//
void CmdLine::AddOption(CmdOption *option)
{
  cmd_options[option->opt] = option;
  cmd_option_list.push_back(option);
}

void CmdLine::AddOption(const char *opt, bool *ptr, bool val, const char *description)
{
  AddOption(new BoolCmdOption(opt,ptr,val,description));
}

void CmdLine::AddOption(const char *opt, char *ptr, char val, const char *description)
{
  AddOption(new CmdOptionTemplate<char>(opt,ptr,val,description));
}

void CmdLine::AddOption(const char *opt, int *ptr, int val, const char *description)
{
  AddOption(new CmdOptionTemplate<int>(opt,ptr,val,description));
}

void CmdLine::AddOption(const char *opt, long int *ptr, long int val, const char *description)
{
  AddOption(new CmdOptionTemplate<long int>(opt,ptr,val,description));
}

void CmdLine::AddOption(const char *opt, unsigned long int *ptr, unsigned long int val, const char *description)
{
  AddOption(new CmdOptionTemplate<unsigned long int>(opt,ptr,val,description));
}

void CmdLine::AddOption(const char *opt, float *ptr, float val, const char *description)
{
  AddOption(new CmdOptionTemplate<float>(opt,ptr,val,description));
}

void CmdLine::AddOption(const char *opt, double *ptr, double val, const char *description)
{
  AddOption(new CmdOptionTemplate<double>(opt,ptr,val,description));
}

void CmdLine::AddOption(const char *opt, char **ptr, const char *val, const char *description)
{
  AddOption(new StrCmdOption(opt,ptr,val,description));
}



//---------------------------------------------------------------------------------------------------//
// END CLASS CmdLine                                                                                 //
//---------------------------------------------------------------------------------------------------//








//---------------------------------------------------------------------------------------------------//
// CLASS CmdLineWithOperations                                                                       //
//---------------------------------------------------------------------------------------------------//


//------Constructor-----------
//
CmdLineWithOperations::CmdLineWithOperations()
{

}


//------Destructor-----------
//
CmdLineWithOperations::~CmdLineWithOperations()
{
  for (OperationMap::iterator it=cmd_operations.begin(); it!=cmd_operations.end(); it++) delete it->second;
}


//------AddOperation-----------
//
void CmdLineWithOperations::AddOperation(string operation, string usage, string description, string details, string examples)
{
  if (cmd_operations.find(operation)!=cmd_operations.end()) { fprintf(stderr, "Error: [CmdLine::AddOperation] operation is already defined!\n"); exit(1); }
  cmd_operations[operation] = new cmd_info(usage,description,details,examples);
  cmd_operations_list.push_back(operation);
}



//------SetCurrentOperation-----------
//
void CmdLineWithOperations::SetCurrentOperation(string operation)
{
  current_cmd_operation = operation;
}



//-----OperationSummary----------
//
void CmdLineWithOperations::OperationSummary(string usage, string description, bool keep_order)
{
  cout << '\n';
  cout << "USAGE: \n";
  cout << "  " << program_name << " " << usage << '\n';
  cout << '\n';
  if (version!="") {
    cout << "VERSION: \n";
    cout << "  " << version << '\n';
    cout << '\n';
  }
  cout << "DESCRIPTION: \n";
  cout << "  " << description << '\n'; 
  cout << '\n';
  cout << "OPERATION: \n";
  if (keep_order==true) for (OperationList::iterator it=cmd_operations_list.begin(); it!=cmd_operations_list.end(); it++) printf("  %-15s %s\n", it->c_str(), cmd_operations[*it]->description.c_str());
  else for (OperationMap::iterator it=cmd_operations.begin(); it!=cmd_operations.end(); it++) printf("  %-15s %s\n", it->first.c_str(), it->second->description.c_str());
  cout << '\n';
}



//------OperationUsage-----------
//
void CmdLineWithOperations::OperationUsage()
{
  OperationMap::iterator it = cmd_operations.find(current_cmd_operation);
  if (it==cmd_operations.end()) { fprintf(stderr, "Error: [CmdLine::AddOperation] operation not found!\n"); exit(1); }
  cout << '\n';
  cout << "USAGE: \n";
  cout << "  " << program_name << " " << it->first << " " << it->second->usage << '\n';
  cout << '\n';
  cout << "DESCRIPTION: \n";
  cout << "  " << it->second->description << '\n'; 
  cout << '\n';
  if (it->second->details!="") {
    cout << "DETAILS: \n";
    cout << "  " << it->second->details << '\n'; 
    cout << '\n';
  }
  if (it->second->examples!="") {
    cout << "EXAMPLES: \n";
    cout << "  " << it->second->examples << '\n'; 
    cout << '\n';
  }
  cout << "OPTIONS: \n";
  Print();
  cout << '\n';
}



//---------------------------------------------------------------------------------------------------//
// END CLASS CmdLineWithOperations                                                                   //
//---------------------------------------------------------------------------------------------------//







//------------------------------------------------------------------------------------------------//
// CLASS Fields                                                                                   //
//------------------------------------------------------------------------------------------------//
// ResetPointer    | Reset the buffer pointer                                                     //
// GetNextToken    | Return next token and advance the buffer pointer                             //
// GetMatrix       | Convert a buffer to a 2-dimensional array of char[]                          //
// GetFloatMatrix  | Convert a buffer to a 2-dimensional array of float                           //
// GetVector       | Convert a buffer to a vector of char[]                                       //
// GetFloatVector  | Convert a buffer to a vector of float                                        //
// GetIntVector    | Convert a buffer to a vector of int                                          //
// Clean           | Free memory occupied by the buffer                                           //
//------------------------------------------------------------------------------------------------//


//--------Fields--------
//
Fields::Fields(char *buffer, const char *delimiters) 
{
  this->buffer = buffer;
  buffer_size = strlen(buffer);
  
  n_rows = 0;
  n_columns = 0;
  int n_col = 0;
  int k = 0;
  while (1) {
    while ((buffer[k]!=0)&&(buffer[k]==' ')) buffer[k++] = 0;            // skip blanks before token
    if (buffer[k]==0) break;

    while ((buffer[k]!=0)&&(strchr(delimiters,buffer[k])==NULL)) ++k;    // read token

    if (n_rows==0) ++n_columns;                                          // update number of columns
    else ++n_col;

    while ((buffer[k]!=0)&&(buffer[k]==' ')) buffer[k++] = 0;            // skip blanks after token

    if (buffer[k]==0) break;
    else if (buffer[k]=='\n') {                                          // check for newline character
      if ((n_rows>0)&&(n_col!=n_columns)) { 
        //fprintf(stderr, "Warning <Fields>: inconsistent number of columns in row #%i!\n", n_rows); 
        //exit(1); 
      }
      ++n_rows;
      n_col = 0;
      buffer[k++] = 0;
    }
  }
  if ((n_columns>0)&&(n_rows==0)) n_rows = 1;
  pointer = buffer;
}



//-----~Fields-----
//
Fields::~Fields() { }



//-----ResetPointer-----
//
void Fields::ResetPointer()
{
  pointer = buffer;
}



//-----GetNextToken-----
//
char *Fields::GetNextToken()
{
  if (pointer-buffer>buffer_size) return NULL;
  while (pointer[0]==0) ++pointer;
  char *token = pointer;
  pointer += strlen(token);
  return token;
}

 

//-----GetPointers-----
//
char ***Fields::GetPointers()
{
  char ***A;
  ALLOCATE2D(A,n_rows,n_columns,char *)
  for (int r=0,j=0; r<n_rows; r++) {
    for (int c=0; c<n_columns; c++) {
      while (buffer[j]==0) j++;
      A[r][c] = buffer + j;
      j += strlen(A[r][c]);
    }
  }
  return A;
}



//-----GetMatrix-----
//
char ***Fields::GetMatrix(int *n, int *m)
{
  char ***A;
  ALLOCATE2D(A,n_rows,n_columns,char *)
  for (int r=0,j=0; r<n_rows; r++) {
    for (int c=0; c<n_columns; c++) {
      while (buffer[j]==0) j++;
      A[r][c] = StrCopy(buffer + j);
      j += strlen(A[r][c]);
    }
  }
  *n = n_rows;
  *m = n_columns;
  return A;
}



//-----GetFloatMatrix-----
//
float **Fields::GetFloatMatrix(int *n, int *m)
{
  float **A;
  ALLOCATE2D(A,n_rows,n_columns,float)
  Progress PRG("<Fields>: extracting float matrix...",n_rows);
  for (int r=0,j=0; r<n_rows; r++) {
    for (int c=0; c<n_columns; c++) {
      while (buffer[j]==0) j++;
      A[r][c] = atof(buffer + j);
      j += strlen(buffer + j);
    }
    PRG.Check();
  }
  PRG.Done();
  *n = n_rows;
  *m = n_columns;
  return A;
}



//-----GetIntMatrix-----
//
int **Fields::GetIntMatrix(int *n, int *m)
{
  int **A;
  ALLOCATE2D(A,n_rows,n_columns,int)
  int T = 0;
  for (int r=0,j=0; r<n_rows; r++) {
    for (int c=0; c<n_columns; c++) {
      while (buffer[j]==0) j++;
      A[r][c] = atoi(buffer + j);
      j += strlen(buffer + j);
    }
    if (_MESSAGES_) if (--T<0) { T=500; fprintf(stderr, "%8i%%\r", (int)(100*(r+1)/n_rows)); }
  }
  if (_MESSAGES_) fprintf(stderr, "%8i%%\n", 100);
  *n = n_rows;
  *m = n_columns;
  return A;
}



//-----GetVector-----
//
char **Fields::GetVector(int *n_items)
{
  *n_items = n_rows * n_columns;
  char **A;
  ALLOCATE1D(A,*n_items,char *)
  for (int i=0,j=0; i<*n_items; i++) {
    while (buffer[j]==0) j++;
    A[i] = StrCopy(buffer + j);
    j += strlen(A[i]);
  }
  return A;
}



//-----GetFloatVector-----
//
float *Fields::GetFloatVector(int *n_items)
{
  *n_items = n_rows * n_columns;
  float *A;
  ALLOCATE1D(A,*n_items,float)
  for (int i=0,j=0; i<*n_items; i++) {
    while (buffer[j]==0) j++;
    A[i] = atof(buffer + j);
    j += strlen(buffer + j);
  }
  return A;
}



//-----GetIntVector-----
//
int *Fields::GetIntVector(int *n_items)
{
  *n_items = n_rows * n_columns;
  int *A;
  ALLOCATE1D(A,*n_items,int)
  for (int i=0,j=0; i<*n_items; i++) {
    while (buffer[j]==0) j++;
    A[i] = atoi(buffer + j);
    j += strlen(buffer + j);
  }
  return A;
}



//-----Find-----
//
bool Fields::Find(char *s)
{
  ResetPointer();
  char *token;
  while ((token=GetNextToken())!=NULL) if (strcmp(token,s)==0) return true;
  return false;
}


//-----Clean-----
//
void Fields::Clean()
{
  free(buffer);
}


//------------------------------------------------------------------------------------------------//
// CLASS Fields                                                                                   //
//------------------------------------------------------------------------------------------------//












//------------------------------------------------------------------------------------------------//
// CLASS IntList                                                                                  //
//------------------------------------------------------------------------------------------------//
// Print    | Print the list elements                                                             //
// Add      | Add a new element                                                                   //
// SortAdd  | Add a new element keeping the list sorted                                           //
// Match    | Match list to an integer vector                                                     //
//------------------------------------------------------------------------------------------------//

IntList::IntList() 
{ 
  n = 0; 
  list = last = NULL;
}

IntList::~IntList()
{
  if (list) FreeList(list);
}

void IntList::FreeList(IntNode *L)
{
  if (L->next) FreeList(L->next);
  else free(L);
}

void IntList::Print()
{
  IntNode *curr = list;
  while (curr) {
    printf(" %i", curr->val);
    curr = curr->next;
  }
  printf("\n");
}

void IntList::Print(FILE *F)
{
  IntNode *curr = list;
  while (curr) {
    fprintf(F, " %i", curr->val);
    curr = curr->next;
  }
  fprintf(F, "\n");
}

bool IntList::Match(int *X, int m)
{
  IntNode *curr = list;
  for (int j=0; j<m; j++) {
    if (curr==NULL) return false;
    if (X[j]!=curr->val) return false;
    curr = curr->next;
  }
  return true;
}

void IntList::Add(int val)
{
  IntNode *new_node = (IntNode *) malloc(sizeof(IntNode));
  if (list) {
    if (last->val==val) return;
    last->next = new_node;
    new_node->val = val;
    new_node->next = NULL;
    last = new_node;
  }
  else {
    list = last = new_node;
    new_node->val = val;
    new_node->next = NULL;
  }
  n++;
}

void IntList::SortAdd(int val)
{
  IntNode *new_node = (IntNode *) malloc(sizeof(IntNode));
  new_node->val = val;
  
  if (list) {
    if ((val<=list->val)||list->next==NULL) {
      if (val<list->val) {
        new_node->next = list;
        list = new_node;
      }
      else {
        list->next = new_node;
        new_node->next = NULL;
        last = new_node;
      }
    }
    else {
      IntNode *curr = list;
      while (curr->next&&(val>curr->next->val)) curr = curr->next;
      if (curr->next) {
        new_node->next = curr->next;
        curr->next = new_node;
      }
      else {
        last->next = new_node;
        new_node->next = NULL;
        last = new_node;
      }
    }
  }
  else {
    list = last = new_node;
    new_node->next = NULL;
  }
  n++;
}

//------------------------------------------------------------------------------------------------//
// CLASS IntList                                                                                  //
//------------------------------------------------------------------------------------------------//










//------------------------------------------------------------------------------------------------//
// SYSTEM COMMANDS                                                                                //
//------------------------------------------------------------------------------------------------//
// CmdExecute     | Execute unix shell commands (safe version)                                    //
// CmdCountLines  | Count the number of lines                                                     //
// IsDirectory    | Checks if given file is a directory                                           //
//------------------------------------------------------------------------------------------------//

//-----CmdExecute-----
//
char *CmdExecute(char *command)
{
  // first check if temporary file already exists
  char temp_file[1000];
  do { sprintf(temp_file, "_temp_%i_", rand()%10000); } while (Exists(temp_file)); 

  // execute the command and store result in temporary file
  char *cmd;
  ALLOCATE1D(cmd,strlen(command)+strlen(temp_file)+20,char);
  sprintf(cmd,"%s > %s",command,temp_file);
  if (system(cmd)==-1) { fprintf(stderr, "Error executing command '%s'!\n", cmd); exit(1); }

  // store contents of temporary file in buffer
  FILE *F;
  if ((F=fopen(temp_file,"r")) == 0) { fprintf(stderr,"Can't open file '%s'!\n",temp_file); exit(1); }
  fseek(F,0,SEEK_END);
  long int buffer_size = ftell(F);
  rewind(F);
  char *buffer;
  ALLOCATE1D(buffer,buffer_size+1,char);
  fread(buffer,1,buffer_size,F);
  buffer[buffer_size] = 0;
  fclose(F);

  // remove temporary and return buffer
  sprintf(cmd,"rm %s",temp_file);
  if (system(cmd)==-1) { fprintf(stderr, "Error executing command '%s'!\n", cmd); exit(1); }
  return buffer;
}



//-----CmdCountLines-----
//
long int CmdCountLines(char *command) 
{
  char *buffer = CmdExecute(command);
  long int n = StrCountLines(buffer);
  FREE1D(buffer);
  return n;
}



//-----IsDirectory-----
//
int IsDirectory(char *file)
{
  struct stat info;
  stat(file, &info);  
  return (info.st_mode & S_IFMT)==S_IFDIR;
}




//-------CountLines-----------
//
unsigned long int CountLines(char *file, unsigned long int buffer_size)
{
  FileBuffer *buffer = CreateFileBuffer(file,buffer_size);
  unsigned long int n = buffer->CountLines();
  delete buffer;
  return n;
}














//------------------------------------------------------------------------------------------------//
// INDEX VECTORS                                                                                  //
//------------------------------------------------------------------------------------------------//
// IPartition      | Partition vector into groups, return indices                                 //
// IPrint          | Print index vector                                                           //
// IApply          | Index-based apply                                                            //
//------------------------------------------------------------------------------------------------//

//-----IPartition----
//
int **IPartition(int *X, int n, int m)
{
  // NOTE: assumes min(X) = 0, max(X) = m-1
  
  int *N, **K;
  ALLOCATE1D_INIT(N,m,int,0);          // partition sizes
  ALLOCATE1D(K,m,int *);               // partition indices

  // count group sizes
  for (int i=0; i<n; i++) {
    int c = X[i];
    if ((c<0)||(c>=m)) { fprintf(stderr, "<IPartition>: value out of bounds!\n"); exit(1); }
    N[c]++;
  }

  // construct partition
  for (int c=0; c<m; c++) { 
    ALLOCATE1D(K[c],N[c]+1,int); 
    K[c][0] = N[c]; 
  }
  for (int c=0; c<m; c++) N[c] = 1;
  for (int i=0; i<n; i++) {
    int c = X[i];
    K[c][N[c]] = i;
    N[c]++;
  }

  // clean up
  FREE1D(N);

  return K;
}



//----IPrint----
//
void IPrint(int *X)
{
  for (int i=1; i<=X[0]; i++) printf("%i ", X[i]);
  printf("\n");
}



//----IApply----
//
float IApply(float **A, int *I, float (*func)(float *,int))
{
  float *V;
  int n = I[0]*I[0];
  ALLOCATE1D(V,n,float);
  for (int i=0,j=1; j<=I[0]; j++) 
    for (int k=1; k<=I[0]; k++,i++) V[i] = A[I[j]][I[k]];
  float y = (*func)(V,n);
  FREE1D(V);
  return y;
}



//----IApply----
//
float IApply(float **A, int *X, int N, float (*func)(float *,int))
{
  // complement indices of X	
  bool *b;
  ALLOCATE1D_INIT(b,N,bool,true);
  for (int x=1; x<=X[0]; x++) b[X[x]] = false;

  // extract sub-matrix
  float *V;
  int n = X[0]*(N-X[0]);
  ALLOCATE1D(V,n,float);
  for (int i=0,x=1; x<=X[0]; x++) 
    for (int k=0; k<N; k++) if (b[k]) V[i++] = A[X[x]][k];
  float z = (*func)(V,n);

  // clean up
  FREE1D(V);
  FREE1D(b);

  return z;
}











//------------------------------------------------------------------------------------------------//
// SEQUENCES = unsorted                                                                           //
// MULTISETS = sorted SEQUENCES                                                                   //
//------------------------------------------------------------------------------------------------//
// Seq               | Construct a sequence from a string                                         //
// SeqSort           | Sort a sequence, effectively converting it into a multiset                 //
// SeqCopy           | Create a copy of the sequence                                              //
// SeqUnique         | Converts a sequence into the corresponding set                             //
// SeqIntersect      | Returns the intersection of two multisets                                  //
// SeqUnion          | Returns the union of two sequences                                         //
// SeqProduct        | Returns the inner product of two multisets                                 //
// SeqCheckIntersect | Checks if the intersection of two multisets is not empty                   //
// SeqPrint          | Display                                                                    //
//------------------------------------------------------------------------------------------------//

//-----Seq------
//
Sequence Seq(const char *s, const char delimiter)
{
  char *scopy = StrCopy(s);
  int n = CountTokens(scopy,delimiter);
  long int *z;
  ALLOCATE1D(z,n+1,long int);
  z[0] = n;
  char *p = scopy;
  for (int k=1; k<=n; k++) z[k] = atol(GetNextToken(&p,delimiter));
  //SeqPrint(z);
  free(scopy);
  
  return z; 
}



//-----Seq------
//
Sequence Seq(long int n)
{
  long int *z;
  ALLOCATE1D(z,n+1,long int);
  z[0] = n;
  
  return z; 
}



//----------Seq2Bag---------
//
Sequence Seq2Bag(Sequence Q)
{
  // count uniq elements in S
  Sequence S = SeqCopy(Q);
  SeqSort(S);
  long int n = 0;
  for (long int k=1; k<=S[0]; ) {
    long int x = S[k];
    while ((k<=S[0])&&(S[k]==x)) k++;
    n++;
  }
 
  // create bag 
  Sequence B = Seq(2*n);
  B[0] = 2*n;
  for (long int k=1,r=1; k<=S[0]; r+=2) {
    long int x = S[k];
    long int y = 0;
    while ((k<=S[0])&&(S[k]==x)) { k++; y++; }
    B[r] = x;
    B[r+1] = y;
  }

  FREE1D(S);

  return B;
}  




//-----SeqSort------
//
void SeqSort(Sequence X)
{
  VectorSort(&X[1],X[0]);
}



//-----SeqIsSorted------
//
bool SeqIsSorted(Sequence X)
{
  for (long int z=2; z<=X[0]; z++) if (X[z]<X[z-1]) return false;
  return true;
}



//-----SeqCopy------
//
Sequence SeqCopy(Sequence X)
{
  Sequence Y;
  ALLOCATE1D(Y,X[0]+1,long int);
  for (long int z=0; z<=X[0]; z++) Y[z] = X[z];
  return Y;
}



//-----SeqUnique--------
//
void SeqUnique(Sequence X)
{
  int n = (int)X[0];
  if (n==0) return;

  SeqSort(X);
  int k = 1;
  for (int i=2; i<=n; i++) {
    while ((i<=n)&&(X[i]==X[k])) i++;                // if constant, ignore
    if (i<=n) X[++k] = X[i];                         // insert new unique element
  }
  X[0] = k;
}


//-----SeqIntersect--------
//
Sequence SeqIntersect(Sequence X, Sequence Y)
{
  int n = (int)X[0];
  int m = (int)Y[0];
  Sequence Z;
  ALLOCATE1D(Z,n+m+1,long int);
  int k = 0;

  for (int i=1,j=1; (i<=n)&&(j<=m); ) {
    if (X[i]<Y[j]) i++;
    else if (Y[j]<X[i]) j++;
    else {
      long int c = X[i];
      for (; (i<=n)&&(X[i]==c); i++) Z[++k] = X[i];
      for (; (j<=m)&&(Y[j]==c); j++) Z[++k] = Y[j];
    }
  }

  Z[0] = k;
  return Z;
}



//-----SeqUnion--------
//
Sequence SeqUnion(Sequence X, Sequence Y)
{
  int n = (int)X[0];
  int m = (int)Y[0];
  
  Sequence Z;
  ALLOCATE1D(Z,n+m+1,long int);
  Z[0] = n+m;
  int k = 1;
  for (int i=1; i<=n; i++) Z[k++] = X[i];
  for (int j=1; j<=m; j++) Z[k++] = Y[j];
  
  return Z;
}



//-----SeqProduct--------
//
float SeqProduct(Sequence X, Sequence Y)
{
  int n = (int)X[0];
  int m = (int)Y[0];

  float prod = 0;
  for (int i=1,j=1; (i<=n)&&(j<=m); ) {
    if (X[i]<Y[j]) i++;
    else if (Y[j]<X[i]) j++;
    else {
      long int c = X[i];
      int a = 0; for (; (i<=n)&&(X[i]==c); i++) a++;
      int b = 0; for (; (j<=m)&&(Y[j]==c); j++) b++;
      prod += a*b;
    }
  }

  return prod;
}



//-----SeqCheckIntersect--------
//
bool SeqCheckIntersect(Sequence X, Sequence Y)
{
  int n = (int)X[0];
  int m = (int)Y[0];

  for (int i=1,j=1; (i<=n)&&(j<=m); ) {
    if (X[i]<Y[j]) i++;
    else if (Y[j]<X[i]) j++;
    else return true;
  }

  return false;
}



//-----SeqPrint------
//
void SeqPrint(Sequence X, FILE *stream)
{
  for (long int z=1; z<=X[0]; z++) fprintf(stream, "%ld ", X[z]);
}


//-----BagPrint------
//
void BagPrint(Sequence X, FILE *stream)
{
  for (long int z=1; z<=X[0]; z+=2) fprintf(stream, "%ld:%ld ", X[z], X[z+1]);
}














//------------------------------------------------------------------------------------------------//
// VECTORS                                                                                        //
//------------------------------------------------------------------------------------------------//
// Vec               | Construct a vector from a string                                           //
// VecCorr           | Vector correlation                                                         //
// VecPrint          | Display                                                                    //
//------------------------------------------------------------------------------------------------//


//-----Vec------
//
Vector Vec(char *s, char delimiter)
{
  char *scopy = StrCopy(s);
  long int n = CountTokens(scopy,delimiter);
  Vector V;
  ALLOCATE1D(V,n+1,float);
  V[0] = (float)n;
  char *p = scopy;
  for (long int k=1; k<=n; k++) {
    char *s = GetNextToken(&p,delimiter);
    LowerCase(s);
    if (strcmp(s,"nan")!=0) V[k] = atof(s);
    else V[k] = numeric_limits<float>::quiet_NaN();   // 0.0/0.0;
  }
  free(scopy);
  
  return V; 
}



//-----Vec------
//
Vector Vec(long int n)
{
  Vector V;
  ALLOCATE1D(V,n+1,float);
  V[0] = (float)n;
  
  return V; 
}



//-----VecCorr------
//
float VecCorr(Vector V1, Vector V2)
{
  long int n1 = (long int)V1[0];
  long int n2 = (long int)V2[0];
  if (n1!=n2) { fprintf(stderr, "<VecCorr>: vector sizes (%ld,%ld) are not equal!\n", n1, n2); exit(1); }  
  return Corr(&V1[1],&V2[1],n1);
}




//-----VecDiff------
//
Vector VecDiff(Vector V1, Vector V2)
{
  long int n1 = (long int)V1[0];
  long int n2 = (long int)V2[0];
  if (n1!=n2) { fprintf(stderr, "<VecCorr>: vector sizes (%ld,%ld) are not equal!\n", n1, n2); exit(1); }  
  Vector V = Vec(n1);
  for (long int k=1; k<=n1; k++) 
    if ((V1[k]==V1[k])&&(V2[k]==V2[k])) V[k] = V1[k]-V2[k];
    else V[k] = numeric_limits<float>::quiet_NaN();   // 0.0/0.0;
  return V;
}




//----------VecNorm--------------
//
Vector VecNorm(Vector V)
{
  long int n = (long int)V[0];
  float avg = VectorAvg(&V[1],n);
  float std = VectorStd(&V[1],n);
  Vector VV = Vec(n);
  for (long int k=1; k<=n; k++) 
    if (V[k]==V[k]) VV[k] = (V[k]-avg)/std;
    else VV[k] = V[k];
  return VV;
}




//-----VecPrint------
//
void VecPrint(Vector V, FILE *stream)
{
  for (long int z=1; z<=(long int)V[0]; z++) 
    if (V[z]==V[z]) fprintf(stream, "%f ", V[z]);
    else fprintf(stream, "NaN ");
}

















//------------------------------------------------------------------------------------------------//
// MATRIX STATISTICS                                                                              //
//------------------------------------------------------------------------------------------------//
// ColumnConst   | Check whether column vectors are constant                                      //
// ColumnSum     | Compute the sum of the values of the columns of matrix M of size n x m         //
// ColumnAvg     | Compute the mean values of the columns of matrix M of size n x m               //
// ColumnStd     | Compute the standard deviation of the columns of matrix M of size n x m        //
// ColumnMax     | Compute the maximum of the columns of matrix M of size n x m                   //
// ColumnMin     | Compute the minimum of the columns of matrix M of size n x m                   //
// ColumnCount   | Count nonempty values of the columns of matrix M of size n x m                 //
//------------------------------------------------------------------------------------------------//

//---ColumnConst-----
//
template <class type> bool *ColumnConst(type **M, unsigned long int n, unsigned long int m)
{
  bool *b;
  ALLOCATE1D_INIT(b,m,bool,true);
  
  for (unsigned long int j=0; j<m; j++) {
    type c = M[0][j];
    for (unsigned long int i=1; i<n; i++) if (M[i][j]!=c) { b[j] = false; break; }
  }

  return b;
}
template bool *ColumnConst<float>(float **, unsigned long int, unsigned long int);
template bool *ColumnConst<double>(double **, unsigned long int, unsigned long int);



//---ColumnSum-----
//
template <class type> type *ColumnSum(type **M, unsigned long int n, unsigned long int m)
{
  type *sum;
  ALLOCATE1D_INIT(sum,m,type,0);
  for (unsigned long int j=0; j<m; j++) 
    for (unsigned long int i=0; i<n; i++) if (M[i][j]==M[i][j]) sum[j] += M[i][j];
  return sum;
}
template float *ColumnSum<float>(float **, unsigned long int, unsigned long int);
template double *ColumnSum<double>(double **, unsigned long int, unsigned long int);
template int *ColumnSum<int>(int **, unsigned long int, unsigned long int);
template long int *ColumnSum<long int>(long int **, unsigned long int, unsigned long int);
template long long int *ColumnSum<long long int>(long long int **, unsigned long int, unsigned long int);
template unsigned long int *ColumnSum<unsigned long int>(unsigned long int **, unsigned long int, unsigned long int);





//---ColumnAvg-----
//
template <class type> type *ColumnAvg(type **M, unsigned long int n, unsigned long int m)
{
  type *avg;
  unsigned long int *count;
  ALLOCATE1D_INIT(avg,m,type,0);
  ALLOCATE1D_INIT(count,m,unsigned long int,0);
  for (unsigned long int j=0; j<m; j++) 
    for (unsigned long int i=0; i<n; i++) if (M[i][j]==M[i][j]) { avg[j] += M[i][j]; count[j]++; }
  for (unsigned long int j=0; j<m; j++) avg[j] /= count[j];
  FREE1D(count);
  return avg;
}
template float *ColumnAvg<float>(float **, unsigned long int, unsigned long int);
template double *ColumnAvg<double>(double **, unsigned long int, unsigned long int);
//template int *ColumnAvg<int>(int **, unsigned long int, unsigned long int);
//template long int *ColumnAvg<long int>(long int **, unsigned long int, unsigned long int);



//---ColumnStd-----
//
template <class type> type *ColumnStd(type **M, unsigned long int n, unsigned long int m)
{
  type *avg = ColumnAvg(M,n,m);
  type *std;
  unsigned long int *count;
  ALLOCATE1D_INIT(std,m,type,0);
  ALLOCATE1D_INIT(count,m,unsigned long int,0);
  for (unsigned long int j=0; j<m; j++) 
    for (unsigned long int i=0; i<n; i++) if (M[i][j]==M[i][j]) { std[j] += M[i][j]*M[i][j]; count[j]++; }
  for (unsigned long int j=0; j<m; j++) std[j] = (type)sqrt((double)(std[j]/count[j] - avg[j]*avg[j]));
  FREE1D(avg);
  FREE1D(count);
  return std;
}
template float *ColumnStd<float>(float **, unsigned long int, unsigned long int);
template double *ColumnStd<double>(double **, unsigned long int, unsigned long int);
template int *ColumnStd<int>(int **, unsigned long int, unsigned long int);
template long int *ColumnStd<long int>(long int **, unsigned long int, unsigned long int);
template long long int *ColumnStd<long long int>(long long int **, unsigned long int, unsigned long int);




//---ColumnMax-----
//
template <class type> type *ColumnMax(type **M, unsigned long int n, unsigned long int m)
{
  type *max;
  ALLOCATE1D(max,m,type);
  for (unsigned long int j=0; j<m; j++) {
    max[j] = M[0][j];
    for (unsigned long int i=0; i<n; i++) if ((M[i][j]==M[i][j])&&(M[i][j]>max[j])) max[j] = M[i][j];
  }
  return max;
}
template float *ColumnMax<float>(float **, unsigned long int, unsigned long int);
template double *ColumnMax<double>(double **, unsigned long int, unsigned long int);
template int *ColumnMax<int>(int **, unsigned long int, unsigned long int);
template long int *ColumnMax<long int>(long int **, unsigned long int, unsigned long int);
template long long int *ColumnMax<long long int>(long long int **, unsigned long int, unsigned long int);



//---ColumnMin-----
//
template <class type> type *ColumnMin(type **M, unsigned long int n, unsigned long int m)
{
  type *min;
  ALLOCATE1D(min,m,type);
  for (unsigned long int j=0; j<m; j++) {
    min[j] = M[0][j];
    for (unsigned long int i=1; i<n; i++) if ((M[i][j]==M[i][j])&&(M[i][j]<min[j])) min[j] = M[i][j];
  }
  return min;
}
template float *ColumnMin<float>(float **, unsigned long int, unsigned long int);
template double *ColumnMin<double>(double **, unsigned long int, unsigned long int);
template int *ColumnMin<int>(int **, unsigned long int, unsigned long int);
template long int *ColumnMin<long int>(long int **, unsigned long int, unsigned long int);
template long long int *ColumnMin<long long int>(long long int **, unsigned long int, unsigned long int);




//---ColumnCount-----
//
int *ColumnCount(float **M, unsigned long int n, unsigned long int m)
{
  int *count;
  ALLOCATE1D_INIT(count,m,int,0);
  for (unsigned long int j=0; j<m; j++) 
    for (unsigned long int i=0; i<n; i++) if (IsNaN(M[i][j])==false) count[j]++;
  return count;
}









//------------------------------------------------------------------------------------------------//
// MATRIX OPERATIONS                                                                              //
//------------------------------------------------------------------------------------------------//
// ColumnNorm      | Normalize column vectors using mean and standard deviation                   //
// ColumnNormMean  | Normalize column vectors to a zero mean                                      //
//------------------------------------------------------------------------------------------------//


//___ColumnNorm______________
void ColumnNorm(float **M, int n, int m)
{
  float *avg = ColumnAvg(M,n,m);
  float *std = ColumnStd(M,n,m);
  for (int i=0; i<n; i++)
    for (int j=0; j<m; j++) if (M[i][j]==M[i][j]) M[i][j] = (M[i][j]-avg[j])/std[j];
  FREE1D(avg);
  FREE1D(std);
}


//___ColumnNormMean______________
void ColumnNormMean(float **M, int n, int m)
{
  float *avg = ColumnAvg(M,n,m);
  for (int i=0; i<n; i++)
    for (int j=0; j<m; j++)
      M[i][j] -= avg[j];
  FREE1D(avg);
}












//------------------------------------------------------------------------------------------------//
// CLASS Sequences : manages a set of sequences                                                   //
//------------------------------------------------------------------------------------------------//

//----- Constructor ------
//
Sequences::Sequences(char *file, bool uniq)
{
  // initialize
  FileBufferText buffer(file);
  n_sequences = buffer.CountLines();
  n_symbols = 0;
  
  // allocate memory
  ALLOCATE1D(SCORE,n_sequences,float);
  ALLOCATE1D(POS,n_sequences,Sequence);
  ALLOCATE1D(SET,n_sequences,Sequence);

  // read data
  for (long int n=0; n<n_sequences; n++) {
    char *inp = buffer.Next();
    int n_fields = CountTokens(inp,'\t');
    if (n_fields>3) { fprintf(stderr, "<Sequences>: line %ld, no more than 3 fields are permitted!\n", n+1); exit(1); }
    if (n_fields==3) SCORE[n] = atof(GetNextToken(&inp,'\t'));
    else SCORE[n] = 0;
    if (n_fields>=2) { char *s = GetNextToken(&inp,'\t'); POS[n] = Seq(s); }
    else { POS[n] = Seq(1); POS[n][1] = n; }
    SET[n] = Seq(inp);
    if (uniq) SeqUnique(SET[n]);
    for (long int k=1; k<=SET[n][0]; k++) {
      if (SET[n][k]<0) { fprintf(stderr, "<Sequences>: line %ld, negative numbers are not allowed!\n", n+1); exit(1); }
      if (SET[n][k]+1>n_symbols) n_symbols = SET[n][k]+1;
    }
  }
}



//----- Destructor ------
//
Sequences::~Sequences()
{
  FREE1D(SCORE);
  FREE2D(POS,n_sequences);
  FREE2D(SET,n_sequences);
}



//----- GetFreqs ------
//
float *Sequences::GetFreqs()
{
  float *q;
  ALLOCATE1D_INIT(q,n_symbols,float,0.0);
  
  for (long int n=0; n<n_sequences; n++)
     for (long int k=1; k<=SET[n][0]; k++) q[SET[n][k]]++;

  float sumq = 0;
  for (long int s=0; s<n_symbols; s++) sumq += q[s];
  for (long int s=0; s<n_symbols; s++) q[s] /= sumq;
  
  return q;
}



//----- Intersect ------
//
Sequence Sequences::Intersect(Sequence list)
{
  if (list[0]<=0) return Seq(""); 
  
  Sequence Z = SeqCopy(SET[list[1]]);
  if (SeqIsSorted(Z)==false) { fprintf(stderr, "Not sorted!\n"); exit(1); }
  for (long int k=2; k<=list[0]; k++) {
    Sequence Q = SeqIntersect(Z,SET[list[k]]);
    FREE1D(Z);
    Z = Q;  
    if (SeqIsSorted(Z)==false) { fprintf(stderr, "Not sorted!\n"); exit(1); }
  }
  SeqUnique(Z);
  
  return Z;
}

//---------------------------------------------------------------------------------//
// END CLASS Sequences                                                             //
//---------------------------------------------------------------------------------//











//------------------------------------------------------------------------------------------------//
//                                                                                                //
// GPL-dependent code start here                                                                  //
//                                                                                                //
//------------------------------------------------------------------------------------------------//
  
  
  
  

//------------------------------------------------------------------------------------------------//
// PERMUTATIONS AND RESAMPLING                                                                    //
//------------------------------------------------------------------------------------------------//


//------Corr_pvalue------
//
double Corr_pvalue(float *A, float *B, unsigned long int M)
{
  double Ex = 0, Ey = 0;
  double Ex2 = 0, Ey2 = 0;
  double Exy = 0;
  int C = 0; 
  for (unsigned long int i=0; i<M; i++) 
    if ((A[i]==A[i])&&(B[i]==B[i])) {
      Ex += A[i];
      Ex2 += pow((double)A[i],2.0);
      Ey += B[i];
      Ey2 += pow((double)B[i],2.0);
      Exy += A[i]*B[i];
      C++;
    }
  Ex = Ex/C;
  Ey = Ey/C;
  Ex2 = Ex2/C;
  Ey2 = Ey2/C;
  Exy = Exy/C;
  double y = fabs((Exy-Ex*Ey)/sqrt((Ex2-pow(Ex,2.0))*(Ey2-pow(Ey,2.0))));
  if (C<3) return 1.0;
  return gsl_cdf_tdist_Q(y/sqrt((1-y*y)/(C-2)),C-2);
}




//------Resample------
//
unsigned long int *Resample(gsl_rng *rnd_generator, unsigned long int n)
{
  unsigned long int *I;
  ALLOCATE1D(I,n,unsigned long int);
  for (unsigned long int k=0; k<n; k++) I[k] = gsl_rng_uniform_int(rnd_generator,n);
  return I;
}


//------Corr_resample------
//
double Corr_resample(gsl_rng *rnd_generator, float *x, float *y, unsigned long int n, unsigned long int n_tests, double *std)
{
  double mean = 0;
  double var = 0;
  for (unsigned long int t=0; t<n_tests; t++) {
    unsigned long int *I = Resample(rnd_generator,n);
    double c = Corr(x,y,I,n);
    mean += c;
    var += c*c;
    free(I);
  }
  mean = mean/n_tests;
  var = var/n_tests - mean*mean;
  *std = var<=0 ? 0.0 : sqrt(var);
  return mean;
}









//------------------------------------------------------------------------------------------------//
// PERMUTATION TESTS                                                                              //
//------------------------------------------------------------------------------------------------//

//_____PermuteTestCorrelation______________
float *PermuteTestCorrelation(float *X, float *Y, int n, int T)
{
  // normalize
  VectorNorm(X,n);
  VectorNorm(Y,n);

  // make a copy of vector Y
  float *Z;
  ALLOCATE1D(Z,n,float);
  for (int i=0; i<n; i++) Z[i] = Y[i];

  // initialize random generator
  gsl_rng_env_setup();
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

  //for (int i=0; i<n; i++) printf("%f ",Z[i]); printf("\n");

  float *K;
  ALLOCATE1D(K,T,float);

  for (int t=0; t<T; t++) {
    gsl_ran_shuffle(r,Z,n,sizeof(float));
    K[t] = InnerProduct(X,Z,n)/n;
  }

  return K;
}


//____AutoPermute____________
float AutoPermute(float *X, int n, int T)
{
  // make a copy of vector X
  float *Z;
  ALLOCATE1D(Z,n,float);
  for (int i=0; i<n; i++) Z[i] = X[i];

  // initialize random generator
  gsl_rng_env_setup();
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

  // compute p-value
  float pval = 0;
  for (int t=0; t<T; t++) {
    gsl_ran_shuffle(r,Z,n,sizeof(float));
    pval += InnerProduct(X,Z,n)/n;
  }

  return pval/T;
}



//----------InitRandomGenerator------------
//
gsl_rng *InitRandomGenerator(unsigned long int seed)
{
  gsl_rng_env_setup();
  gsl_rng *RANDOM_GENERATOR = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(RANDOM_GENERATOR,seed);
  return RANDOM_GENERATOR;
}












//---------------------------------------------------------------------------------//
// CLASS Matrix                                                                    //
//---------------------------------------------------------------------------------//



//-------Constructor--------
//
Matrix::Matrix(char *file, bool integers, bool verbose)
{
  // Load file or standard input
  FileBuffer *buffer;
  long int n_lines;
  if (file==NULL) buffer = new FileBufferText(LoadStdIn(&n_lines));
  else { buffer = new FileBufferText(file); n_lines = buffer->CountLines(); }
  if (n_lines==0) { fprintf(stderr, "Empty matrix!\n"); exit(1); }
  this->integers = integers;

  // check/read column labels
  char *inp = buffer->Next();
  if (strchr(inp,'\t')!=NULL) {
    has_row_labels = has_col_labels = true;
    n_rows = n_lines - 1;
    GetNextToken(&inp,'\t');
    n_cols = CountTokens(inp,' ');
    col_labels = new string[n_cols];
    row_labels = new string[n_rows];
    for (int c=0; c<n_cols; c++) col_labels[c] = GetNextToken(&inp,' ');
    inp = buffer->Next();      
  }
  else {
    has_row_labels = has_col_labels = false;
    n_rows = n_lines;
    n_cols = CountTokens(inp,' ');
    col_labels = row_labels = NULL;
  }

  // allocate matrix memory
  if (verbose) fprintf(stderr, "* Found %d rows and %d columns; row labels = %s; column labels = %s.\n", n_rows, n_cols, has_row_labels?"ON":"OFF", has_col_labels?"ON":"OFF");
  if (integers==true) { val = NULL; ALLOCATE2D(ival,n_rows,n_cols,int); }
  else { ival = NULL; ALLOCATE2D(val,n_rows,n_cols,float); }

  // read contents
  Progress PRG("Reading matrix entries...",n_rows);
  for (int r=0; r<n_rows; r++,inp=buffer->Next()) {
    if (has_row_labels) row_labels[r] = GetNextToken(&inp,'\t');
    long int n_elem;
    if (integers==true) ival[r] = ReadIntVec(inp,&n_elem);
    else val[r] = ReadFloatVec(inp,&n_elem);
    if (n_elem!=n_cols) { fprintf(stderr, "Line %d: vector length should be %d!\n", r+1+has_col_labels, n_cols); exit(1); }
    PRG.Check();
  }
  PRG.Done();   

  // cleanup
  delete buffer;
}



//-------Destructor--------
//
Matrix::~Matrix()
{
  if (row_labels!=NULL) delete [] row_labels;
  if (col_labels!=NULL) delete [] col_labels;
  FREE2D(val,n_rows);
  FREE2D(ival,n_rows);
}



//-------Print--------
//
void Matrix::Print(char *fmt)
{
  if (has_col_labels==true) {
    printf("\t");
    for (int c=0; c<n_cols; c++) cout << col_labels[c] << ' ';
    printf("\n");
  }
  Progress PRG("Printing matrix...",n_rows);
  for (int r=0; r<n_rows; r++) {
    if (has_row_labels) cout << row_labels[r] << '\t';
    PrintRow(r,fmt);
    printf("\n");
    PRG.Check();
  }
  PRG.Done();   
}



//-------PrintSparse--------
//
void Matrix::PrintSparse(char *fmt, int offset)
{
  Progress PRG("Printing in sparse format...",n_rows);
  for (int r=0; r<n_rows; r++) {
    if (has_row_labels) cout << row_labels[r] << '\t';
    for (int c=0; c<n_cols; c++) 
      if (integers==true) { if ((ival[r][c]==ival[r][c])&&(ival[r][c]!=0)) { printf("%d:", offset+c); printf(fmt, ival[r][c]); printf(" "); } }
      else { if ((val[r][c]==val[r][c])&&(val[r][c]!=0)) { printf("%d:", offset+c); printf(fmt, val[r][c]); printf(" "); } }
    printf("\n");
    PRG.Check();
  }
  PRG.Done();
}



//-------PrintEncoded--------
//
void Matrix::PrintEncoded(int n_bins)
{
  Progress PRG("Printing in encoded format...",n_rows);
  for (int r=0; r<n_rows; r++) {
    if (has_row_labels) cout << row_labels[r] << '\t';
    for (int c=0; c<n_cols; c++) if (ival[r][c]==ival[r][c]) printf("%ld ", (long int)c*n_bins+(long int)ival[r][c]);
    printf("\n");
    PRG.Check();
  }
  PRG.Done();
}



//-------PrintTranspose--------
//
void Matrix::PrintTranspose(char *fmt)
{
  if (has_row_labels==true) {
    printf("\t");
    for (int r=0; r<n_rows; r++) cout << row_labels[r] << ' ';
    printf("\n");
  }
  Progress PRG("Printing matrix...",n_cols);
  for (int c=0; c<n_cols; c++) {
    if (has_col_labels) cout << col_labels[c] << '\t';
    for (int r=0; r<n_rows; r++) {
      if (integers==true) { if (ival[r][c]==ival[r][c]) printf(fmt, ival[r][c]); else printf("NaN"); }
      else { if (val[r][c]==val[r][c]) printf(fmt, val[r][c]); else printf("NaN"); }
      printf(" ");
    }
    printf("\n");
    PRG.Check();
  }
  PRG.Done();   
}



//-------PrintColStats--------
//
void Matrix::PrintColStats(char *fmt)
{
  // compute statistics
  bool *b= ColumnConst(val,n_rows,n_cols);
  float *sum = ColumnSum(val,n_rows,n_cols);
  float *min = ColumnMin(val,n_rows,n_cols);
  float *max = ColumnMax(val,n_rows,n_cols);
  float *avg = ColumnAvg(val,n_rows,n_cols);
  float *std = ColumnStd(val,n_rows,n_cols);

  // print statistics
  if (has_col_labels==true) {
    printf("\t");
    for (int c=0; c<n_cols; c++) cout << col_labels[c] << ' ';
    printf("\n");
  }
  printf("BMP\t"); for (long int c=0; c<n_cols; c++) printf("%i ", b[c]); printf("\n");
  printf("SUM\t"); for (long int c=0; c<n_cols; c++) { printf(fmt, sum[c]); printf(" "); }; printf("\n");
  printf("MIN\t"); for (long int c=0; c<n_cols; c++) { printf(fmt, min[c]); printf(" "); }; printf("\n");
  printf("MAX\t"); for (long int c=0; c<n_cols; c++) { printf(fmt, max[c]); printf(" "); }; printf("\n");
  printf("AVG\t"); for (long int c=0; c<n_cols; c++) { printf(fmt, avg[c]); printf(" "); }; printf("\n");
  printf("STD\t"); for (long int c=0; c<n_cols; c++) { printf(fmt, std[c]); printf(" "); }; printf("\n");

  // free memory
  FREE1D(b);
  FREE1D(sum);
  FREE1D(avg);
  FREE1D(std);
  FREE1D(min);
  FREE1D(max);
}



//-------PrintRowStats--------
//
void Matrix::PrintRowStats(char *fmt)
{
  // print statistics
  printf("\tSUM\tMIN\tMAX\tAVG\tSTD\n");
  Progress PRG("Printing row statistics...",n_rows);
  for (int r=0; r<n_rows; r++) {
    if (has_row_labels) cout << row_labels[r] << '\t';
    //printf("%i", VectorConst(val[r],n_cols)); printf(" ");
    printf(fmt, VectorSum(val[r],n_cols)); printf("\t");
    printf(fmt, VectorMin(val[r],n_cols)); printf("\t");
    printf(fmt, VectorMax(val[r],n_cols)); printf("\t");
    printf(fmt, VectorAvg(val[r],n_cols)); printf("\t");
    printf(fmt, VectorStd(val[r],n_cols));
    printf("\n");
    PRG.Check();
  }
  PRG.Done();   
}



//-------PrintRowPairs--------
//
void Matrix::PrintRowPairs(char *fmt)
{
  Progress PRG("Printing row pairs...",n_rows);
  for (int r1=0; r1<n_rows; r1++) {
    for (int r2=r1+1; r2<n_rows; r2++) {
      if (has_row_labels) cout << row_labels[r1] << '\t' << row_labels[r2] << '\t';
      PrintRow(r1,fmt);
      cout << '\t';
      PrintRow(r2,fmt);
      printf("\n");
    }
    PRG.Check();
  }
  PRG.Done();   
}



//-------PrintRow--------
//
void Matrix::PrintRow(int r, char *fmt)
{
  for (int c=0; c<n_cols; c++) {
    if (integers==true) { if (ival[r][c]==ival[r][c]) printf(fmt, ival[r][c]); else printf("NaN"); }
    else { if (val[r][c]==val[r][c]) printf(fmt, val[r][c]); else printf("NaN"); }
    printf(" ");
  }
}



//-------PrintColLabels--------
//
void Matrix::PrintColLabels()
{
  if (has_row_labels==true) cout << '\t';
  for (int c=0; c<n_cols; c++) cout << col_labels[c] << ' ';
  cout << '\n';
}



//-------ApplyCutoff(float)--------
//
void Matrix::ApplyCutoff(float cutoff, bool upper)
{
  Progress PRG("Applying cutoff...",n_rows);
  for (int r=0; r<n_rows; r++) {
    for (int c=0; c<n_cols; c++)
      if (val[r][c]==val[r][c]) {
        if (upper) val[r][c] = val[r][c]>cutoff ? cutoff:val[r][c];
        else val[r][c] = val[r][c]<cutoff ? cutoff:val[r][c];
      }
    PRG.Check();
  }
  PRG.Done();
}



//-------ApplyCutoff(int)--------
//
void Matrix::ApplyCutoff(int cutoff, bool upper)
{
  Progress PRG("Applying cutoff...",n_rows);
  for (int r=0; r<n_rows; r++) {
    for (int c=0; c<n_cols; c++)
      if (ival[r][c]==ival[r][c]) {
        if (upper) ival[r][c] = ival[r][c]>cutoff ? cutoff:ival[r][c];
        else ival[r][c] = ival[r][c]<cutoff ? cutoff:ival[r][c];
      }
    PRG.Check();
  }
  PRG.Done();
}



//-------Test(float)--------
//
void Matrix::Test(float cutoff, bool equal, bool greater)
{
  Progress PRG("Applying test...",n_rows);
  for (int r=0; r<n_rows; r++) {
    for (int c=0; c<n_cols; c++) {
      if (val[r][c]==val[r][c]) {
        if ((equal==true)&&(val[r][c]==cutoff)) { val[r][c] = 1; continue; }
        if ((greater==true)&&(val[r][c]>cutoff)) { val[r][c] = 1; continue; }
        if ((greater==false)&&(val[r][c]<cutoff)) { val[r][c] = 1; continue; }
        val[r][c] = 0;
      }
    }
    PRG.Check();
  }
  PRG.Done();   
}



//-------Test(int)--------
//
void Matrix::Test(int cutoff, bool equal, bool greater)
{
  Progress PRG("Applying test...",n_rows);
  for (int r=0; r<n_rows; r++) {
    for (int c=0; c<n_cols; c++) {
      if (ival[r][c]==ival[r][c]) {
        if ((equal==true)&&(ival[r][c]==cutoff)) { ival[r][c] = 1; continue; }
        if ((greater==true)&&(ival[r][c]>cutoff)) { ival[r][c] = 1; continue; }
        if ((greater==false)&&(ival[r][c]<cutoff)) { ival[r][c] = 1; continue; }
        ival[r][c] = 0;
      }
    }
    PRG.Check();
  }
  PRG.Done();   
}



//-------NormRows--------
//
void Matrix::NormRows(bool range)
{
  Progress PRG("Normalizing rows...",n_rows);
  for (int r=0; r<n_rows; r++) {
    if (range) VectorNormRange(val[r],n_cols);
    else VectorNorm(val[r],n_cols);
    PRG.Check();
  }
  PRG.Done();
}


//-------NormCols--------
//
void Matrix::NormCols()
{
  float *avg = ColumnAvg(val,n_rows,n_cols);
  float *std = ColumnStd(val,n_rows,n_cols);
  Progress PRG("Normalizing columns...",n_rows);
  for (long int r=0; r<n_rows; r++) {
    for (long int c=0; c<n_cols; c++) 
      if (val[r][c]==val[r][c]) val[r][c] = (val[r][c]-avg[c])/std[c];
    PRG.Check();
  }
  PRG.Done();
  FREE1D(avg);
  FREE1D(std);
}


//-------Shuffle--------
//
void Matrix::Shuffle()
{
  gsl_rng *RANDOM_GENERATOR = InitRandomGenerator(time(NULL));
  Progress PRG("Shuffling rows...",n_rows);
  for (long int r=0; r<n_rows; r++) {
    if (integers==true) gsl_ran_shuffle(RANDOM_GENERATOR,ival[r],n_cols,sizeof(int));
    else gsl_ran_shuffle(RANDOM_GENERATOR,val[r],n_cols,sizeof(float));
    PRG.Check();
  }
  PRG.Done();
  gsl_rng_free(RANDOM_GENERATOR);
}


//-------Multiply--------
//
void Matrix::Multiply(float coeff)
{
  Progress PRG("Multiplying...",n_rows);
  for (int r=0; r<n_rows; r++) {
    for (int c=0; c<n_cols; c++) if (val[r][c]==val[r][c]) val[r][c] *= coeff;
    PRG.Check();
  }
  PRG.Done();
}


//---------------------------------------------------------------------------------//
// END CLASS Matrix                                                                //
//---------------------------------------------------------------------------------//







