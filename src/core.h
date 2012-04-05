//
// Copyright (c) 2011 IBM Corporation. 
// All rights reserved. This program and the accompanying materials are made available under the terms of the GNU General Public License v2.0 
// which accompanies this distribution, and is available at http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
//

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <list>
#include "gzstream.h"
#include "sam.h"
using namespace std;





//------------------------------------------------------------------------------------------------//
// Maximum allowed values                                                                         //
//------------------------------------------------------------------------------------------------//

#define MAX_BUFFER_SIZE 100000000 
#define MAX_FUNCTIONS 50





//------------------------------------------------------------------------------------------------//
// External variables                                                                             //
//------------------------------------------------------------------------------------------------//

extern int _MESSAGES_;






//------------------------------------------------------------------------------------------------//
// Predefined types                                                                               //
//------------------------------------------------------------------------------------------------//

typedef void FuncType(char *);
typedef long int *Sequence;
typedef float *Vector;








//------------------------------------------------------------------------------------------------//
// Macros for memory allocation                                                                   //
//------------------------------------------------------------------------------------------------//

#define ALLOCATE(N,TYPE)  (TYPE *) malloc((N)*sizeof(TYPE));
#define CHECK(PTR)        if ((PTR)==NULL) { fprintf(stderr, "Out of memory!\n"); exit(1); }


#define ALLOCATE1D(PTR,N,TYPE)			\
{						\
  PTR = (TYPE *) malloc((N)*sizeof(TYPE));	\
  CHECK(PTR);					\
}


#define ALLOCATE1D_INIT(PTR,N,TYPE,VAL)							\
{											\
  ALLOCATE1D(PTR,N,TYPE);								\
  for (unsigned long int k=0; k<(unsigned long int)(N); k++) PTR[k] = VAL;		\
}


#define ALLOCATE2D(PTR,N,M,TYPE)					\
{ 									\
  PTR = (TYPE **) malloc((N)*sizeof(TYPE *));				\
  CHECK(PTR);								\
  for (unsigned long int i=0; i<(unsigned long int)(N); i++) {		\
    PTR[i] = (TYPE *) malloc((M)*sizeof(TYPE));				\
    CHECK(PTR[i]);							\
  }									\
}


#define FREE1D(PTR)				\
  if (PTR!=NULL) free(PTR);


#define FREE2D(PTR,N)							\
{									\
  if ((PTR)!=NULL) {							\
    for (unsigned long int i=0; i<(unsigned long int)(N); i++) 		\
      if ((PTR)[i]!=NULL) free((PTR)[i]);				\
    free(PTR);                                                          \
  }									\
}



#define CHECK_BOUNDS(X,DESCR,LBOUND,UBOUND)			\
  { if (((X)<(LBOUND)) || ((X)>=UBOUND)) { fprintf(stderr, "Error: %s out of bounds [%d]!\n", DESCR, X); } }



#define MESSAGES(STATE) { _MESSAGES_ = STATE; }












//------------------------------------------------------------------------------------------------//
// CLASS FuncTable                                                                                //
//------------------------------------------------------------------------------------------------//

class FuncTable
{
 public:
  FuncTable() { nfunc = 0; }
  ~FuncTable() { }
  void Add(char *name, FuncType *f, char *description);
  void Execute(char *name, char *arg);
  void Help();
  FuncType *Lookup(char *name);
 private:
  int nfunc;
  char name[MAX_FUNCTIONS][40];
  char description[MAX_FUNCTIONS][200];
  FuncType *f[MAX_FUNCTIONS];
};

void Interpreter(FuncTable *FTABLE);

//------------------------------------------------------------------------------------------------//
// CLASS FuncTable                                                                                //
//------------------------------------------------------------------------------------------------//










//------------------------------------------------------------------------------------------------//
// CLASS Fields                                                                                   //
//------------------------------------------------------------------------------------------------//
// ResetPointer    | Reset the buffer pointer                                                     //
// GetNextToken    | Return next token and advance the buffer pointer                             //
// GetMatrix       | Convert a buffer to a 2-dimensional array of char[]                          //
// GetFloatMatrix  | Convert a buffer to a 2-dimensional array of float                           //
// GetIntMatrix    | Convert a buffer to a 2-dimensional array of int                             //
// GetVector       | Convert a buffer to a vector of char[]                                       //
// GetFloatVector  | Convert a buffer to a vector of float                                        //
// GetIntVector    | Convert a buffer to a vector of int                                          //
// Clean           | Free memory occupied by the buffer                                           //
//------------------------------------------------------------------------------------------------//

class Fields 
{
 public:
  Fields(char *buffer, const char *delimiters);
  ~Fields();

  void ResetPointer();
  char *GetNextToken(); 
  void Clean();

  // get pointers
  char ***GetPointers();

  // convert to matrix
  char ***GetMatrix(int *n_rows, int *n_columns);
  float **GetFloatMatrix(int *n_rows, int *n_columns);
  int **GetIntMatrix(int *n_rows, int *n_columns);

  // convert to vector
  char **GetVector(int *n_items);
  float *GetFloatVector(int *n_items);
  int *GetIntVector(int *n_items);

  // find a token
  bool Find(char *s);

 public:
  int n_rows, n_columns;
  long int buffer_size;
  char *buffer;
  char *pointer;
};


//------------------------------------------------------------------------------------------------//
// CLASS Fields                                                                                   //
//------------------------------------------------------------------------------------------------//









//------------------------------------------------------------------------------------------------//
// CLASS IntList                                                                                  //
//------------------------------------------------------------------------------------------------//
// Print | Print the list elements                                                                //
// Add   | Add a new element                                                                      //
//------------------------------------------------------------------------------------------------//

struct IntNode
{
  int val;
  IntNode *next;
};

class IntList
{
 public:
  IntList();
  ~IntList();
  void Print();
  void Print(FILE *F);
  void Add(int val);
  void SortAdd(int val);
  bool Match(int *X, int n);

  int n;
  IntNode *list;
  IntNode *last;

 private:
  void FreeList(IntNode *L);

};


//------------------------------------------------------------------------------------------------//
// CLASS IntList                                                                                  //
//------------------------------------------------------------------------------------------------//













//------------------------------------------------------------------------------------------------//
// CLASS Progress                                                                                 //
//------------------------------------------------------------------------------------------------//
//!  This class is used to report progress of loop computations.
//------------------------------------------------------------------------------------------------//
class Progress 
{

 public:
  //! Class constructor.
  /*!
    \param msg 		message to be displayed in stderr.
    \param max_count	if greater than 1, the progress is shown as percentage
  */
  Progress(const char *msg, long int max_count);

  //! Class constructor.
  Progress();

  ~Progress();

  //! Prints message and initializes counter
  /*!
    \param msg 		message to be displayed in stderr.
    \param max_count	if greater than 1, the progress is shown as percentage
  */
  void Init(const char *msg, long int max_count);

  //! Updates counter and reports progress every 1sec.
  /*!
    Call this method inside the loop you want to monitor at the end of each iteration.
  */
  void Check();

  //! Prints final count.
  /*!
    Call this method after the loop. 
  */
  void Done();

 private:
  char *msg;			//!< message
  long int max_count;		//!< maximum number of iterations in the loop
  long int count;		//!< current iteration
  time_t TIME;			//!< time stamp (updated every second)
};


//------------------------------------------------------------------------------------------------//
// CLASS Progress                                                                                 //
//------------------------------------------------------------------------------------------------//





 
//------------------------------------------------------------------------------------------------//
// CLASS FileBuffer                                                                               //
//------------------------------------------------------------------------------------------------//
//! Abstract class for reading lines from file.
//------------------------------------------------------------------------------------------------//
class FileBuffer
{
 public:
  //! Empty class constructor.
  FileBuffer() { };

  //! Class destructor
  virtual ~FileBuffer();

  // methods
  long int CountLines();			//!< Counts the number of lines in the file
  virtual void Reset() = 0;			//!< Resets the file pointer (obviously this does not work for standard input)
  char *Get();						//!< Returns a pointer to the current line
  virtual char *Next() = 0;			//!< Read the next line

  // data
  bool is_stdin;					//!< true if reading from standard input
  unsigned long int n_line;			//!< keeps track of line number
  char *file_name;					//!< file name
  char *BUFFER;						//!< where input line is stored
  unsigned long int BUFFER_SIZE;	//!< buffer size (automatically adjusted during execution to accommodate any line size)
};
//------------------------------------------------------------------------------------------------//
// END CLASS FileBuffer                                                                           //
//------------------------------------------------------------------------------------------------//








//------------------------------------------------------------------------------------------------//
// CLASS FileBufferText                                                                           //
//------------------------------------------------------------------------------------------------//
//! Class for reading lines from text file or standard input.
//------------------------------------------------------------------------------------------------//
class FileBufferText : public FileBuffer
{
 public:
  //! Class constructor.
  /*!
    \param file 		file name, if 'NULL' standard input is read instead
    \param buffer_size 		buffer size (automatically adjusted during execution to accommodate any line size)
  */
  FileBufferText(const char *file, unsigned long int buffer_size=10000);

  //! Class constructor.
  /*!
    \param file_ptr 		file pointer
    \param buffer_size 		buffer size (automatically adjusted during execution to accommodate any line size)
  */
  FileBufferText(FILE *file_ptr, unsigned long int buffer_size=10000);

  //! Class destructor
  virtual ~FileBufferText();

  // methods
  virtual void Reset();				//!< Resets the file pointer (obviously this does not work for standard input)
  virtual char *Next();				//!< Read the next line

  // data
  FILE *file_ptr;					//!< file pointer
};
//------------------------------------------------------------------------------------------------//
// END CLASS FileBufferText                                                                       //
//------------------------------------------------------------------------------------------------//








//------------------------------------------------------------------------------------------------//
// CLASS FileBufferGZ                                                                             //
//------------------------------------------------------------------------------------------------//
//! Class for reading lines from file in GZ format.
//------------------------------------------------------------------------------------------------//
class FileBufferGZ : public FileBuffer
{
 public:
  //! Class constructor.
  /*!
    \param file 		file name, if 'NULL' standard input is read instead
    \param buffer_size 		buffer size (automatically adjusted during execution to accommodate any line size)
  */
  FileBufferGZ(const char *file, unsigned long int buffer_size=10000);

  //! Class destructor
  virtual ~FileBufferGZ();

  virtual void Reset();				//!< Resets the file pointer (obviously this does not work for standard input)
  virtual char *Next();				//!< Read the next line
  bool Read(char *buffer, unsigned long int buffer_size); 		//!< Auxiliary function for reading into buffer

  igzstream *file_stream;				//!< input stream (used for gz files)
};
//------------------------------------------------------------------------------------------------//
// END CLASS FileBufferGZ                                                                         //
//------------------------------------------------------------------------------------------------//








//------------------------------------------------------------------------------------------------//
// CLASS FileBufferBAM                                                                            //
//------------------------------------------------------------------------------------------------//
//! Class for reading lines from file in BAM format.
//------------------------------------------------------------------------------------------------//
class FileBufferBAM : public FileBuffer
{
 public:
  //! Class constructor.
  /*!
    \param file 		file name, if 'NULL' standard input is read instead
    \param buffer_size 		buffer size (automatically adjusted during execution to accommodate any line size)
  */
  FileBufferBAM(const char *file, unsigned long int buffer_size=10000);
  
  //! destructor
  virtual ~FileBufferBAM();

  virtual void Reset();				//!< Resets the file pointer (obviously this does not work for standard input)
  virtual char *Next();				//!< Read the next line

  samfile_t *samfile_ptr;			//!< pointer to SAM/BAM file
  bam1_t *bam_ptr;					//!< pointer to BAM structure
};
//------------------------------------------------------------------------------------------------//
// END CLASS FileBufferBAM                                                                        //
//------------------------------------------------------------------------------------------------//




FileBuffer *CreateFileBuffer(const char *file, unsigned long int buffer_size=10000);










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
char *TokenPreprocess(char *inp, char *inp_delims, char out_delim);
int CountTokens(char *s, char TOKEN_DELIM=' ');
char *GetNextToken(char **pBUF, const char *DELIMS);
char *GetNextToken(char **pBUF, char DELIM);
char *GetNextToken(char **inp);
char *SkipBlank(char *inp);
void Read(char *s, float *v);
void Read(char *s, double *v);
void Read(char *s, int *v);
void Read(char *s, long int *v);
void Read(char *s, unsigned long int *v);
void Read(char *s, long long int *v);
template <class T> void ReadVec(char *s, T **vec, long int *n, char sep=' ');
template <class T> void ReadSparseVec(char *s, T **vec, long int *n);
float *ReadFloatVec(char *s, long int *n, char sep=' ');
double *ReadDoubleVec(char *s, long int *n, char sep=' ');
int *ReadIntVec(char *s, long int *n, char sep=' ');











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
// VectorNonZero | Checks if vector <> 0                                                          //
// VectorCount   | Count nonempty values                                                          //
// VectorVariab  | Compute vector variability measure                                             //
// VectorEntropy | Compute vector entropy                                                         //
// Hist          | Histogram                                                                      //
// GetThreshold  | Determine a threshold for an exponentially decreasing distribution             //
//------------------------------------------------------------------------------------------------//
template <class type> type VectorAvg(type *v, int n);
template <class type> type VectorStd(type *v, int n);
template <class type> type VectorSum(type *v, int n);
template <class type> type VectorSumSq(type *v, int n);
template <class type> type VectorDiff(type *v, int n);
template <class type> type VectorProd(type *v, int n);
template <class type> type VectorMax(type *v, int n);
template <class type> type VectorMin(type *v, int n);
template <class type> int VectorTest(type *v, int n, bool (*func)(type));
template <class type> type VectorRange(type *v, int n);
template <class type> bool VectorNonZero(type *X, int n);
template <class type> int VectorCount(type *X, int n);
template <class type> float VectorL2(type *X, int n);
float VectorVariab(float *X, int n);
float VectorEntropy(float *X, int n, int nbins);
float *Hist(float *X, int N, int n_bins, float min=0, float max=0);
float GetThreshold(float *X, int n, int delta);









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
float *VectorSmooth(float *x, int n, int d);
void VectorSort(float *V, int n);
void VectorSort(int *V, int n);
void VectorSort(long int *V, int n);
void VectorSort(unsigned long int *V, int n);
void VectorSort(double *V, long int n);
int *VectorRank(float *V, int n);
int *VectorRank(double *V, int n);
int *VectorRank(long int *V, int n);
void VectorChop(float *x, int n, float cutoff);
template <class type> type* VectorCopy(type *x, unsigned long int n);
template <class type> void VectorNorm(type *x, unsigned long int n);
template <class type> void VectorNormSum(type *x, unsigned long int n);
template <class type> void VectorNormRange(type *x, int n);
void VectorNormMean(float *x, int n);
void VectorNormSum2(float *x, int n);
float *VectorCut(float *X, int n, bool *vert_cut);
float *VectorCut(float *X, int n, long int mask, int *nn);
float *Apply(float **A, int N, int M, float (*func)(float *,int));
float **Pairwise(float **A, int N, int M, float (*func)(float *,float *,unsigned long int));
float **Pairwise(float **X, int Nx, float **Y, int Ny, int M, float (*func)(float *,float *,unsigned long int));
float *HistPairwise(float **A, int N, int M, float (*func)(float *,float *,unsigned long int), int n_bins);
float *HistPairwise(float **A, float **B, int Na, int Nb, int M, float (*func)(float *,float *,unsigned long int), int n_bins);
float *Int2Float(int *x, int n);
template <class type> void VectorPrint(type *X, int n, const char *format, bool *vert_cut=NULL, FILE *output=stdout);








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
template <class type> type InnerProduct(type *A, type *B, unsigned long int n);
template <class type> type Covariance(type *A, type *B, unsigned long int n);
template <class type> type Euclidean(type *A, type *B, unsigned long int n);
template <class type> double VectorCorr(type *A, type *B, unsigned long int n);
template <class type> type RelativeEntropy(type *p, type *q, unsigned long int n);
template <class type> type Chi2(type *x, type *avg, unsigned long int n);
template <class type> type Chi2std(type *x, type *avg, type *std, unsigned long int n);
float Jaccard(int *A, int *B, int M);
float Corr(float *A, float *B, unsigned long int M);
float Corr(float *A, float *B, unsigned long int *I, unsigned long int M);
int Hamming1(int *x, int *y, int n);











//------------------------------------------------------------------------------------------------//
// MISCELLANEOUS                                                                                  //
//------------------------------------------------------------------------------------------------//
// Abs_                | Absolute value                                                           //
// Min_                | Minimum of two numbers                                                   //
// Max_                | Maximum of two numbers                                                   //
// IsNaN               | Check for NaN values                                                     //
//------------------------------------------------------------------------------------------------//
float Abs_(float x); 
double Abs_(double x); 
bool IsNaN(float x);
int Min_(int x, int y);
int Max_(int x, int y);
long int Min_(long int x, long int y);
long int Max_(long int x, long int y);
size_t Min_(size_t x, size_t y);
size_t Max_(size_t x, size_t y);








//------------------------------------------------------------------------------------------------//
// VECTOR I/O                                                                                     //
//------------------------------------------------------------------------------------------------//
// Exists      | Check whether a file exists                                                      //
// LoadStdIn   | Load standard input into a buffer                                                //
// LoadFile    | Load a file into a buffer                                                        //
// LoadVectors | Load an array of vectors from a file into memory                                 //
// SaveVectors | Save an array of vectors from memory into a file                                 //
// LoadMatrix  | load a matrix from a file with error checking                                    //
//------------------------------------------------------------------------------------------------//
bool Exists(char *file);
char *LoadStdIn();
FILE *LoadStdIn(long int *n_lines, long int buffer_size=MAX_BUFFER_SIZE);
char *LoadFile(FILE *F);
char *LoadFile(char *file);
float **LoadVectors(char *file, int *N, int *M);
void SaveVectors(char *file, char *format, float **data, int *FILTER, int N, int M);
void SaveVectors(char *file, char *format, float **data, int N, int M);
void SaveVectors(char *file, int **data, int N, int M);
float **LoadMatrix(char *file, long int *n_rows, long int *n_cols, char sep=' ');
double **LoadDoubleMatrix(char *file, long int *n_rows, long int *n_cols);
int **LoadIntMatrix(char *file, long int *n_rows, long int *n_cols);


//! Determines file type: 0=text; 1=gz; 2=bam
int GetFileType(const char *file);





//------------------------------------------------------------------------------------------------//
// SYSTEM COMMANDS                                                                                //
//------------------------------------------------------------------------------------------------//
// CmdExecute     | Execute unix shell commands (safe version)                                    //
// CmdCountLines  | Count the number of lines                                                     //
// IsDirectory    | Checks if given file is a directory                                           //
//------------------------------------------------------------------------------------------------//
char *CmdExecute(char *command);
long int CmdCountLines(char *command);
unsigned long int CountLines(char *file, unsigned long int buffer_size);
int IsDirectory(char *file);






//------------------------------------------------------------------------------------------------//
// STRINGS                                                                                        //
//------------------------------------------------------------------------------------------------//
// StrCopy        | Make a new copy of a string                                                   //
// StrChop        | Removes trailing characters                                                   //
// AddSuffix      | Add a suffix to a string                                                      //
// UpperCase      | Convert to uppercase                                                          //
// LowerCase      | Convert to lowercase                                                          //
// IsSuffix       | Test if suffix matches                                                        //
// StrCountLines  | Count newline characters                                                      //
//------------------------------------------------------------------------------------------------//
char *StrCopy(const char *source);
void StrChop(char *source, char c='\n');
char *AddSuffix(char *source, char *suffix);
void UpperCase(char *s);
void LowerCase(char *s);
bool IsSuffix(char *str, char *sfx);
long int StrCountLines(char *s);
unsigned long int MyRand(unsigned long int N);






char complement(char c);
string ReverseComplement(string &x);
void ReverseGene(char *s, unsigned long int len);




//---------------------------------------------------------------------------------------------------//
// CLASS CmdOption                                                                                   //
//---------------------------------------------------------------------------------------------------//


class CmdOption 
{
  public:
   CmdOption(const char *opt, const char *description);
   virtual ~CmdOption();

   // functions 
   virtual void Init();
   virtual void Print();
   virtual void Read(int *argc, char ***argv);

   // variables   
   char *opt;
   char *description;
};


template <class T> class CmdOptionTemplate : public CmdOption
{
  public:
   CmdOptionTemplate(const char *opt, T *ptr, T val, const char *description);
   ~CmdOptionTemplate();

   // functions 
   virtual void Init();
   virtual void Print();
   virtual void Read(int *argc, char ***argv);

   // variables   
   T *ptr;
   T val;
};


class BoolCmdOption : public CmdOption
{
  public:
   BoolCmdOption(const char *opt, bool *ptr, bool val, const char *description);
   ~BoolCmdOption();

   // functions 
   virtual void Init();
   virtual void Print();
   virtual void Read(int *argc, char ***argv);

   // variables   
   bool *ptr;
   bool val;
};


class StrCmdOption : public CmdOption
{
  public:
   StrCmdOption(const char *opt, char **ptr, const char *val, const char *description);
   ~StrCmdOption();

   // functions 
   virtual void Init();
   virtual void Print();
   virtual void Read(int *argc, char ***argv);

   // variables   
   char **ptr;
   char *val;
};








//---------------------------------------------------------------------------------------------------//
// CLASS CmdLine                                                                                     //
//---------------------------------------------------------------------------------------------------//

class CmdLine
{
  public:
    typedef map<string,CmdOption*> OptionMap;
    typedef list<CmdOption*> OptionList;

    // constructor & destructor    
    CmdLine();
    ~CmdLine();

    // read & print
    void SetProgramName(string program_name, string version="");
    void Init();
    int Read(char **argv, int argc);
    void Print();

    // options
    void AddOption(CmdOption *option);
    void AddOption(const char *opt, char *ptr, char val, const char *description);
    void AddOption(const char *opt, int *ptr, int val, const char *description);
    void AddOption(const char *opt, unsigned long int *ptr, unsigned long int val, const char *description);
    void AddOption(const char *opt, long int *ptr, long int val, const char *description);
    void AddOption(const char *opt, bool *ptr, bool val, const char *description);
    void AddOption(const char *opt, float *ptr, float val, const char *description);
    void AddOption(const char *opt, double *ptr, double val, const char *description);
    void AddOption(const char *opt, char **ptr, const char *val, const char *description);

    // print info
    void Usage(const char *text);
    void Usage(char *prog_name, const char *usage);

    // data
    string program_name;
    string version; 
    OptionMap cmd_options;
    OptionList cmd_option_list;
};

//---------------------------------------------------------------------------------------------------//
// END CLASS CmdLine                                                                                 //
//---------------------------------------------------------------------------------------------------//











//---------------------------------------------------------------------------------------------------//
// CLASS CmdLineWithOperations                                                                       //
//---------------------------------------------------------------------------------------------------//

class CmdLineWithOperations : public CmdLine
{
  public:
    struct cmd_info { 
      string usage, description, details, examples; 
      cmd_info(string &usage, string &description, string &details, string &examples):usage(usage),description(description),details(details),examples(examples) {}
    }; 
    typedef map<string,cmd_info*> OperationMap;

    // constructor & destructor    
    CmdLineWithOperations();
    ~CmdLineWithOperations();

    // operations
    void AddOperation(string operation, string usage, string description, string details="", string examples="");
    void SetCurrentOperation(string operation);

    // print info
    void OperationSummary(string usage, string description);
    void OperationUsage();

    // data
    string current_cmd_operation;
    OperationMap cmd_operations;
};

//---------------------------------------------------------------------------------------------------//
// END CLASS CmdLineWithOperations                                                                   //
//---------------------------------------------------------------------------------------------------//











//------------------------------------------------------------------------------------------------//
// INDEX VECTORS                                                                                  //
//------------------------------------------------------------------------------------------------//
// IPartition      | Partition vector into groups, return indices                                 //
// IPrint          | Print index vector                                                           //
// IApply          | Index-based apply                                                            //
//------------------------------------------------------------------------------------------------//
int **IPartition(int *X, int n, int m);
void IPrint(int *X);
float IApply(float **A, int *I, float (*func)(float *,int));
float IApply(float **A, int *X, int N, float (*func)(float *,int));










//------------------------------------------------------------------------------------------------//
// SEQUENCES = unsorted                                                                           //
// MULTISETS = sorted SEQUENCES                                                                   //
//------------------------------------------------------------------------------------------------//
// Seq               | Construct a sequence from a string                                         //
// Seq2Bag           | Convert a sequence to a bag                                                //
// SeqSort           | Sort a sequence, effectively converting it into a multiset                 //
// SeqCopy           | Create a copy of the sequence                                              //
// SeqUnique         | Converts a sequence into the corresponding set                             //
// SeqIntersect      | Returns the intersection of two multisets                                  //
// SeqUnion          | Returns the union of two sequences                                         //
// SeqProduct        | Returns the inner product of two multisets                                 //
// SeqCheckIntersect | Checks if the intersection of two multisets is not empty                   //
// SeqIsSorted       | Checks if sequence is sorted                                               //
// SeqPrint          | Display                                                                    //
//------------------------------------------------------------------------------------------------//
Sequence Seq(const char *s, const char delimiter=' ');
Sequence Seq(long int n);
Sequence Seq2Bag(Sequence S);
Sequence SeqCopy(Sequence X);
void SeqSort(Sequence X);
void SeqUnique(Sequence X);
Sequence SeqIntersect(Sequence X, Sequence Y);
Sequence SeqUnion(Sequence X, Sequence Y);
float SeqProduct(Sequence X, Sequence Y);
bool SeqIsSorted(Sequence X);
bool SeqCheckIntersect(Sequence X, Sequence Y);
void SeqPrint(Sequence X, FILE *stream=stdout);
void BagPrint(Sequence X, FILE *stream=stdout);









//------------------------------------------------------------------------------------------------//
// VECTORS                                                                                        //
//------------------------------------------------------------------------------------------------//
// Vec               | Construct a vector from a string                                           //
// VecCorr           | Vector correlation                                                         //
// VecPrint          | Display                                                                    //
//------------------------------------------------------------------------------------------------//
Vector Vec(char *s, char delimiter=' ');
Vector Vec(long int n);
float VecCorr(Vector V1, Vector V2);
Vector VecDiff(Vector V1, Vector V2);
Vector VecNorm(Vector V);
void VecPrint(Vector V, FILE *stream=stdout);








//------------------------------------------------------------------------------------------------//
// MATRIX OPERATIONS                                                                              //
//------------------------------------------------------------------------------------------------//
// ColumnConst   | Check whether column vectors are constant                                      //
// ColumnSum     | Compute the sum of the values of the columns of matrix M of size n x m         //
// ColumnAvg     | Compute the mean values of the columns of matrix M of size n x m               //
// ColumnStd     | Compute the standard deviation of the columns of matrix M of size n x m        //
// ColumnMax     | Compute the maximum of the columns of matrix M of size n x m                   //
// ColumnMin     | Compute the minimum of the columns of matrix M of size n x m                   //
// ColumnCount   | Count nonempty values of the columns of matrix M of size n x m                 //
// ColumnNorm      | Normalize column vectors using mean and standard deviation                   //
// ColumnNormMean  | Normalize column vectors to a zero mean                                      //
//------------------------------------------------------------------------------------------------//
template <class type> bool *ColumnConst(type **M, unsigned long int n, unsigned long int m);
template <class type> type *ColumnAvg(type **M, unsigned long int n, unsigned long int m);
template <class type> type *ColumnSum(type **M, unsigned long int n, unsigned long int m);
template <class type> type *ColumnStd(type **M, unsigned long int n, unsigned long int m);
template <class type> type *ColumnMax(type **M, unsigned long int n, unsigned long int m);
template <class type> type *ColumnMin(type **M, unsigned long int n, unsigned long int m);
int *ColumnCount(float **M, unsigned long int n, unsigned long int m);
void ColumnNorm(float **M, int n, int m);
void ColumnNormMean(float **M, int n, int m);












//------------------------------------------------------------------------------------------------//
// CLASS Sequences : manages a set of sequences                                                   //
//------------------------------------------------------------------------------------------------//
class Sequences
{
 public:
  Sequences(char *file, bool uniq=false);
  ~Sequences();

  float *GetFreqs();
  Sequence Intersect(Sequence list);

  long int n_sequences;
  long int n_symbols;

  float *SCORE;
  Sequence *POS, *SET;
};


  





  
  
  
  
  
  
  
//------------------------------------------------------------------------------------------------//
//                                                                                                //
// GPL-dependent code start here                                                                  //
//                                                                                                //
//------------------------------------------------------------------------------------------------//
  
  
  
  
  
  



//------------------------------------------------------------------------------------------------//
// PERMUTATIONS AND RESAMPLING                                                                    //
//------------------------------------------------------------------------------------------------//

gsl_rng *InitRandomGenerator(unsigned long int seed);
float *PermuteTestCorrelation(float *X, float *Y, int n, int T);
float AutoPermute(float *X, int n, int T);
unsigned long int *Resample(gsl_rng *rnd_generator, unsigned long int n);
double Corr_pvalue(float *A, float *B, unsigned long int M);
double Corr_resample(gsl_rng *rnd_generator, float *x, float *y, unsigned long int n, unsigned long int n_tests, double *std);












//---------------------------------------------------------------------------------//
// CLASS Matrix                                                                    //
//---------------------------------------------------------------------------------//

class Matrix
{
 public:
  Matrix(char *file, bool integers=false, bool verbose=false);
  ~Matrix();
  
  // print
  void Print(char *fmt="%.2f");
  void PrintSparse(char *fmt="%.2f", int offset=0);
  void PrintEncoded(int n_bins);
  void PrintTranspose(char *fmt="%.2f");
  void PrintRowStats(char *fmt="%.2f");
  void PrintColStats(char *fmt="%.2f");
  void PrintRow(int r, char *fmt="%.2f");
  void PrintRowPairs(char *fmt="%.2f");
  void PrintColLabels();

  // modify
  void ApplyCutoff(float cutoff, bool upper=false);
  void ApplyCutoff(int cutoff, bool upper=false);
  void Test(float cutoff, bool equal=true, bool greater=true);
  void Test(int cutoff, bool equal=true, bool greater=true);
  void NormRows(bool range=false);
  void NormCols();
  void Shuffle();
  void Multiply(float coeff);

  // data
  bool integers;
  int n_rows, n_cols;
  bool has_row_labels, has_col_labels;
  std::string *row_labels, *col_labels;
  float **val;
  int **ival;
};








//---------------------------------------------------------------------------------//
// CLASS MatrixTemplate                                                            //
//---------------------------------------------------------------------------------//

template <class T> class MatrixTemplate
{
 public:
  MatrixTemplate(char *file, bool verbose=false);
  ~MatrixTemplate();
  
  // functions
  void Print(char *fmt="%.2f");
  void Shrink(char *fmt="%.2f");
  void PrintSparse(char *fmt="%.2f", long int offset=0);
  void PrintEncoded(int n_bins);
  void PrintTranspose(char *fmt="%.2f");
  void PrintRowStats(char *fmt="%.2f");
  void PrintColStats(char *fmt="%.2f");
  void PrintColSums(char *fmt="%.2f");
  void PrintRow(int r, char *fmt="%.2f");
  void PrintRowPairs(char *fmt="%.2f");
  void PrintColLabels();
  void ApplyCutoff(T cutoff, bool upper=false);
  void Test(T cutoff, bool equal=true, bool greater=true);
  void Norm(T avg, T std);
  void NormRows(bool range=false);
  void NormCols();
  void Shuffle(gsl_rng *rnd_generator);
  void Multiply(T coeff);
  void Rel(bool first_column=false);
  void PrintRel2(char *fmt="%.2f");
  void Delete(char *fmt, T cutoff, bool equal=true, bool greater=true);

  // data
  long int n_rows, n_cols;
  bool has_row_labels, has_col_labels;
  std::string *row_labels, *col_labels;
  T **val;
};


//-------Constructor--------
//
template <class T> MatrixTemplate<T>::MatrixTemplate(char *file, bool verbose)
{
  // Load file or standard input
  FileBuffer *buffer;
  long int n_lines;
  if (file==NULL) buffer = new FileBufferText(LoadStdIn(&n_lines));
  else { buffer = new FileBufferText(file); n_lines = buffer->CountLines(); }
  if (n_lines==0) { col_labels = row_labels = NULL; val = NULL; n_rows = n_cols = 0; return; }

  // check/read column labels
  char *inp = buffer->Next();
  if (strchr(inp,'\t')!=NULL) {
    has_row_labels = true;
    has_col_labels = inp[0]=='\t';
    n_rows = n_lines - has_col_labels;
    row_labels = new std::string[n_rows];
    if (has_col_labels==true) {
      GetNextToken(&inp,'\t');
      n_cols = CountTokens(inp,' ');
      col_labels = new std::string[n_cols];
      for (long int c=0; c<n_cols; c++) col_labels[c] = GetNextToken(&inp,' ');
      inp = buffer->Next();
    }
    else {
      col_labels = NULL;
      n_cols = CountTokens(inp,' ');
    }
  }
  else {
    has_row_labels = has_col_labels = false;
    n_rows = n_lines;
    n_cols = CountTokens(inp,' ');
    col_labels = row_labels = NULL;
  }

  // allocate matrix memory
  if (verbose) fprintf(stderr, "* Found %ld lines, %ld rows and %ld columns; row labels = %s; column labels = %s.\n", n_lines, n_rows, n_cols, has_row_labels?"ON":"OFF", has_col_labels?"ON":"OFF");
  ALLOCATE1D(val,n_rows,T *);

  // read contents
  Progress PRG("Reading matrix entries...",n_rows);
  for (long int r=0; r<n_rows; r++,inp=buffer->Next()) {
    bool row_label_found = strchr(inp,'\t')!=NULL;
    if (row_label_found!=has_row_labels) { fprintf(stderr, "Line %ld: row_label_found = %s, has_row_labels = %s\n", r+1+has_col_labels, row_label_found?"YES":"NO", has_row_labels?"YES":"NO"); exit(1); } 
    if (has_row_labels) row_labels[r] = GetNextToken(&inp,'\t');
    long int n_elem;
    ReadVec(inp,&val[r],&n_elem);
    if (n_elem!=n_cols) { fprintf(stderr, "Line %ld: vector length should be %ld!\n", r+1+has_col_labels, n_cols); exit(1); }
    PRG.Check();
  }
  PRG.Done();   

  // cleanup
  delete buffer;
}



//-------Destructor--------
//
template <class T> MatrixTemplate<T>::~MatrixTemplate()
{
  if (row_labels!=NULL) delete [] row_labels;
  if (col_labels!=NULL) delete [] col_labels;
  FREE2D(val,n_rows);
}



//-------Print--------
//
template <class T> void MatrixTemplate<T>::Print(char *fmt)
{
  if (has_col_labels==true) {
    if (has_row_labels) printf("\t");
    for (long int c=0; c<n_cols; c++) std::cout << col_labels[c] << ' ';
    printf("\n");
  }
  Progress PRG("Printing matrix...",n_rows);
  for (long int r=0; r<n_rows; r++) {
    if (has_row_labels) std::cout << row_labels[r] << '\t';
    PrintRow(r,fmt);
    printf("\n");
    PRG.Check();
  }
  PRG.Done();   
}



//-------Shrink--------
//
template <class T> void MatrixTemplate<T>::Shrink(char *fmt)
{
  if (n_rows==0) return;
  bool *keep = new bool[n_cols];
  for (long int c=0; c<n_cols; c++) {
    keep[c] = false;
    for (long int r=0; r<n_rows; r++) if ((val[r][c]==val[r][c])&&(val[r][c]!=0)) { keep[c] = true; break; }
  }
  if (has_col_labels==true) {
    if (has_row_labels) printf("\t");
    for (long int c=0; c<n_cols; c++) if (keep[c]==true) std::cout << col_labels[c] << ' ';
    printf("\n");
  }
  Progress PRG("Printing matrix...",n_rows);
  for (long int r=0; r<n_rows; r++) {
    if (has_row_labels) std::cout << row_labels[r] << '\t';
    for (long int c=0; c<n_cols; c++) if (keep[c]==true) {
      if (val[r][c]==val[r][c]) printf(fmt, val[r][c]); else printf("NaN");
      printf(" ");
    }
    printf("\n");
    PRG.Check();
  }
  PRG.Done();
  delete keep;
}



//-------PrintSparse--------
//
template <class T> void MatrixTemplate<T>::PrintSparse(char *fmt, long int offset)
{
  Progress PRG("Printing in sparse format...",n_rows);
  for (long int r=0; r<n_rows; r++) {
    if (has_row_labels) std::cout << row_labels[r] << '\t';
    for (long int c=0; c<n_cols; c++) 
      if ((val[r][c]==val[r][c])&&(val[r][c]!=0)) { printf("%ld:", offset+c); printf(fmt, val[r][c]); printf(" "); } 
    printf("\n");
    PRG.Check();
  }
  PRG.Done();
}



//-------PrintEncoded--------
//
template <class T> void MatrixTemplate<T>::PrintEncoded(int n_bins)
{
  Progress PRG("Printing in encoded format...",n_rows);
  for (long int r=0; r<n_rows; r++) {
    if (has_row_labels) std::cout << row_labels[r] << '\t';
    for (long int c=0; c<n_cols; c++) if (val[r][c]==val[r][c]) printf("%ld ", (long int)c*n_bins+(long int)val[r][c]);
    printf("\n");
    PRG.Check();
  }
  PRG.Done();
}



//-------PrintTranspose--------
//
template <class T> void MatrixTemplate<T>::PrintTranspose(char *fmt)
{
  if (has_row_labels==true) {
    if (has_col_labels) printf("\t");
    for (long int r=0; r<n_rows; r++) std::cout << row_labels[r] << ' ';
    printf("\n");
  }
  Progress PRG("Printing matrix...",n_cols);
  for (long int c=0; c<n_cols; c++) {
    if (has_col_labels) std::cout << col_labels[c] << '\t';
    for (long int r=0; r<n_rows; r++) {
      if (val[r][c]==val[r][c]) printf(fmt, val[r][c]); else printf("NaN");
      printf(" ");
    }
    printf("\n");
    PRG.Check();
  }
  PRG.Done();   
}



//-------PrintColStats--------
//
template <class T> void MatrixTemplate<T>::PrintColStats(char *fmt)
{
  // compute statistics
  T *sum = ColumnSum(val,n_rows,n_cols);
  T *min = ColumnMin(val,n_rows,n_cols);
  T *max = ColumnMax(val,n_rows,n_cols);
  T *avg = ColumnAvg(val,n_rows,n_cols);
  T *std = ColumnStd(val,n_rows,n_cols);

  // print statistics
  if (has_col_labels==true) {
    printf("\t");
    for (long int c=0; c<n_cols; c++) std::cout << col_labels[c] << ' ';
    printf("\n");
  }
  printf("SUM\t"); for (long int c=0; c<n_cols; c++) { printf(fmt, sum[c]); printf(" "); }; printf("\n");
  printf("MIN\t"); for (long int c=0; c<n_cols; c++) { printf(fmt, min[c]); printf(" "); }; printf("\n");
  printf("MAX\t"); for (long int c=0; c<n_cols; c++) { printf(fmt, max[c]); printf(" "); }; printf("\n");
  printf("AVG\t"); for (long int c=0; c<n_cols; c++) { printf(fmt, avg[c]); printf(" "); }; printf("\n");
  printf("STD\t"); for (long int c=0; c<n_cols; c++) { printf(fmt, std[c]); printf(" "); }; printf("\n");

  // free memory
  FREE1D(sum);
  FREE1D(avg);
  FREE1D(std);
  FREE1D(min);
  FREE1D(max);
}



//-------PrintColSums--------
//
template <class T> void MatrixTemplate<T>::PrintColSums(char *fmt)
{
  // return 0, if matrix is empty
  if (n_rows==0) { printf(fmt, 0); printf("\n"); return; }

  // compute statistics
  T *sum = ColumnSum(val,n_rows,n_cols);

  // print statistics
  if (has_col_labels==true) {
    printf("\t");
    for (long int c=0; c<n_cols; c++) std::cout << col_labels[c] << ' ';
    printf("\n");
  }
  for (long int c=0; c<n_cols; c++) { printf(fmt, sum[c]); printf(" "); }; printf("\n");

  // free memory
  FREE1D(sum);
}



//-------PrintRowStats--------
//
template <class T> void MatrixTemplate<T>::PrintRowStats(char *fmt)
{
  // print statistics
  printf("\tSUM\tMIN\tMAX\tAVG\tSTD\tSTD/AVG\tCOUNT\n");
  Progress PRG("Printing row statistics...",n_rows);
  for (long int r=0; r<n_rows; r++) {
    if (has_row_labels) std::cout << row_labels[r] << '\t';
    printf(fmt, VectorSum(val[r],n_cols)); printf("\t");
    printf(fmt, VectorMin(val[r],n_cols)); printf("\t");
    printf(fmt, VectorMax(val[r],n_cols)); printf("\t");
    printf(fmt, VectorAvg(val[r],n_cols)); printf("\t");
    printf(fmt, VectorStd(val[r],n_cols)); printf("\t");
    printf(fmt, VectorStd(val[r],n_cols)/VectorAvg(val[r],n_cols)); printf("\t");
    printf("%ld", n_cols);
    printf("\n");
    PRG.Check();
  }
  PRG.Done();   
}



//-------PrintRowPairs--------
//
template <class T> void MatrixTemplate<T>::PrintRowPairs(char *fmt)
{
  Progress PRG("Printing row pairs...",n_rows);
  for (long int r1=0; r1<n_rows; r1++) {
    for (long int r2=r1+1; r2<n_rows; r2++) {
      if (has_row_labels) std::cout << row_labels[r1] << '\t' << row_labels[r2] << '\t';
      PrintRow(r1,fmt);
      std::cout << '\t';
      PrintRow(r2,fmt);
      printf("\n");
    }
    PRG.Check();
  }
  PRG.Done();   
}



//-------PrintRow--------
//
template <class T> void MatrixTemplate<T>::PrintRow(int r, char *fmt)
{
  for (long int c=0; c<n_cols; c++) {
    if (val[r][c]==val[r][c]) printf(fmt, val[r][c]); else printf("NaN");
    printf(" ");
  }
}



//-------PrintColLabels--------
//
template <class T> void MatrixTemplate<T>::PrintColLabels()
{
  if (has_row_labels==true) std::cout << '\t';
  for (long int c=0; c<n_cols; c++) std::cout << col_labels[c] << ' ';
  std::cout << '\n';
}



//-------ApplyCutoff--------
//
template <class T> void MatrixTemplate<T>::ApplyCutoff(T cutoff, bool upper)
{
  Progress PRG("Applying cutoff...",n_rows);
  for (long int r=0; r<n_rows; r++) {
    for (long int c=0; c<n_cols; c++)
      if (val[r][c]==val[r][c]) {
        if (upper) val[r][c] = val[r][c]>cutoff ? cutoff:val[r][c];
        else val[r][c] = val[r][c]<cutoff ? cutoff:val[r][c];
      }
    PRG.Check();
  }
  PRG.Done();
}



//-------Test--------
//
template <class T> void MatrixTemplate<T>::Test(T cutoff, bool equal, bool greater)
{
  Progress PRG("Applying test...",n_rows);
  for (long int r=0; r<n_rows; r++) {
    for (long int c=0; c<n_cols; c++) {
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



//-------Delete--------
//
template <class T> void MatrixTemplate<T>::Delete(char *fmt, T cutoff, bool equal, bool greater)
{
  Progress PRG("Deleting rows...",n_rows);
  long int *selected = new long int[n_rows];
  long int n_selected = 1;
  selected[0] = 0;
  if (has_col_labels==true) {
    if (has_row_labels) printf("\t");
    for (long int c=0; c<n_cols; c++) std::cout << col_labels[c] << ' ';
    printf("\n");
  }
  if (has_row_labels) std::cout << row_labels[0] << '\t';
  PrintRow(0,fmt);
  printf("\n");
  for (long int r=1; r<n_rows; r++) {
    bool add = true;
    for (long int z=0; z<n_selected; z++) {
      double corr = VectorCorr(val[r],val[selected[z]],n_cols);
      if ((equal==true)&&(corr==cutoff)) { add = false; break; }
      if ((greater==true)&&(corr>cutoff)) { add = false; break; }
      if ((greater==false)&&(corr<cutoff)) { add = false; break; }
    }
    if (add==true) {
      selected[n_selected++] = r;
      if (has_row_labels) std::cout << row_labels[r] << '\t';
      PrintRow(r,fmt);
      printf("\n");
    }
    PRG.Check();
  }
  PRG.Done();   
  if (_MESSAGES_) fprintf(stderr, "* selected %ld rows.\n", n_selected);
  delete selected;
}



//-------Norm--------
//
template <class T> void MatrixTemplate<T>::Norm(T avg, T std)
{
  Progress PRG("Normalizing rows...",n_rows);
  for (long int r=0; r<n_rows; r++) {
    for (long int c=0; c<n_cols; c++) val[r][c] = (val[r][c]-avg)/std;
    PRG.Check();
  }
  PRG.Done();
}


//-------NormRows--------
//
template <class T> void MatrixTemplate<T>::NormRows(bool range)
{
  Progress PRG("Normalizing rows...",n_rows);
  for (long int r=0; r<n_rows; r++) {
    if (range) VectorNormRange(val[r],n_cols);
    else VectorNorm(val[r],n_cols);
    PRG.Check();
  }
  PRG.Done();
}


//-------NormCols--------
//
template <class T> void MatrixTemplate<T>::NormCols()
{
  T *avg = ColumnAvg(val,n_rows,n_cols);
  T *std = ColumnStd(val,n_rows,n_cols);
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
template <class T> void MatrixTemplate<T>::Shuffle(gsl_rng *rnd_generator)
{
  Progress PRG("Shuffling rows...",n_rows);
  for (long int r=0; r<n_rows; r++) {
    gsl_ran_shuffle(rnd_generator,val[r],n_cols,sizeof(T));
    PRG.Check();
  }
  PRG.Done();
}


//-------Multiply--------
//
template <class T> void MatrixTemplate<T>::Multiply(T coeff)
{
  Progress PRG("Multiplying...",n_rows);
  for (long int r=0; r<n_rows; r++) {
    for (long int c=0; c<n_cols; c++) if (val[r][c]==val[r][c]) val[r][c] *= coeff;
    PRG.Check();
  }
  PRG.Done();
}


//-------Rel--------
//
template <class T> void MatrixTemplate<T>::Rel(bool first_column)
{
  Progress PRG("Fold enrichment...",n_rows);
  for (long int r=0; r<n_rows; r++) {
    for (long int c=n_cols-1; c>=1; c--) if (val[r][c]==val[r][c]) val[r][c] /= first_column==true?val[r][0]:val[r][c-1];
    PRG.Check();
  }
  PRG.Done();
}


//-------PrintRel2--------
//
template <class T> void MatrixTemplate<T>::PrintRel2(char *fmt)
{
  Progress PRG("Fold enrichment...",n_rows);
  if (has_col_labels==true) {
    printf("\t");
    for (long int c0=0; c0<n_cols-1; c0++) 
      for (long int c=c0+1; c<n_cols; c++) std::cout << col_labels[c] << '/' << col_labels[c0] << ' ';
    printf("\n");
  }
  for (long int r=0; r<n_rows; r++) {
    if (has_row_labels) std::cout << row_labels[r] << '\t';
    for (long int c0=0; c0<n_cols-1; c0++) 
      for (long int c=c0+1; c<n_cols; c++) { printf(fmt, val[r][c]/val[r][c0]); printf(" "); }
    printf("\n");
    PRG.Check();
  }
  PRG.Done();
}


//---------------------------------------------------------------------------------//
// END CLASS MatrixTemplate                                                        //
//---------------------------------------------------------------------------------//





