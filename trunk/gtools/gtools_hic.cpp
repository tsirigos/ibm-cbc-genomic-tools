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
#include <sstream>
#include "core.h"
#include "genomic_intervals.h"
using namespace std;



// TODO
// ---------------
// * filter: more refined classification (as I had initially)
// * compare align/filter to imakaev
// * [LATER] fragment reg file: check if non-overlapping (don't care about sorted)
// * [LATER] bin: --matrix-output
// * [LATER] filter at the fragment level: too big, too small (--frag-min-size, --frag-max-size) [this needs to be done before creating the index]


//---------------------------------------------------------------------------------//
// Global variables & constants                                                    //
//---------------------------------------------------------------------------------//

const string PROGRAM = "gtools_hic";
const long int BUFFER_SIZE = 10000;
gsl_rng *RANDOM_GENERATOR;
const char *BIN_BITS = "17,20,23,26";


const int n_pair_types = 10;
enum PairType { UNPAIRED, NO_FRAG, UNMAPPED, MULTIHIT, UNMAPPED_MULTIHIT, SINGLE, FRAG_SELF, DIST_SMALL, FRAG_PROXIMAL, FRAG_DISTAL };
const char *PairTypeStr[] = { "unpaired", "no-fragment", "unmapped", "multihit", "unmapped-multihit", "single", "same-fragment", "small-distance", "proximal", "distal" };

typedef map<long int,bool> MapOfRead1;
typedef map<long int,MapOfRead1*> MapOfRead2;
typedef map<string,MapOfRead2*> MapOfRead;
MapOfRead dupMask;



//---------------------------------------------------------------------------------//
// COMMAND-LINE OPTIONS                                                            //
//---------------------------------------------------------------------------------//

bool VERBOSE;
bool HELP;
char *ALIGN_WORK_DIR;
int ALIGN_MIN_LEN;
int ALIGN_LEN_DIFF;
int BOWTIE_THREADS;
char *BOWTIE_PATH;
char *BOWTIE_INDEX;
char *ALIGNED_REG_FILE;
char *ENZYME_REG_FILE;
bool FILTER_SELF_FRAG;
long int FILTER_MIN_DISTANCE;
long int FILTER_MAX_OFFSET;
bool FILTER_DUPLICATES;
char *FILTER_CLASS_FILE;
long int BIN_SIZE;



//-------InitCmdLine-----------
//
CmdLineWithOperations *InitCmdLine(int argc, char *argv[], int *next_arg)
{
  // initialize
  CmdLineWithOperations *cmd_line = new CmdLineWithOperations(); 
  cmd_line->SetProgramName(PROGRAM,VERSION);

  // set operations
  cmd_line->AddOperation("align", "[OPTIONS] READ1-FASTQ READ2-FASTQ", \
  "Iteratively aligns HiC-seq read pairs to reference genome using bowtie2.", \
  "* Input: FASTQ files\n\
  * Output: aligned reads in SAM format (same order as in fastq files)"\
  );
  
  cmd_line->AddOperation("filter", "[OPTIONS] <ALIGNED-READS>", \
  "Filters HiC-seq read pairs for common experimental artifacts.", \
  "* Input: aligned reads in SAM format (grouped by read-id)\n\
  * Output: filtered read pairs in REG format"\
  );

  cmd_line->AddOperation("bin", "[OPTIONS] <FILTERED-READ-PAIRS>", \
  "Bins filtered read pairs to genomic bins of desired resolution.", \
  "* Input: filtered read pairs in REG format\n\
  * Output: binned read pairs"\
  );

  map<string,int> min_files;
  min_files["align"] = 2;
  min_files["filter"] = 0;
  min_files["bin"] = 0;
  
  // print summary of operations
  if (argc<2) { 
    cmd_line->OperationSummary("OPERATION [OPTIONS] <REGION-SET>", \
    "Pipeline for HiC-seq data analysis. For detailed description and list of options choose an operation and use the --help option.", true); 
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
    cmd_line->AddOption("--work-dir", &ALIGN_WORK_DIR, "", "working directory (required)");
    cmd_line->AddOption("--min-len", &ALIGN_MIN_LEN, 30, "minimum truncated read length");
    cmd_line->AddOption("--len-diff", &ALIGN_LEN_DIFF, 10, "read truncation step");
    cmd_line->AddOption("-p", &BOWTIE_THREADS, 1, "number of threads for bowtie2 run");
    cmd_line->AddOption("--bowtie-path", &BOWTIE_PATH, "bowtie2", "full bowtie2 path (version>=2.1.0)");
    cmd_line->AddOption("--bowtie-index", &BOWTIE_INDEX, "genome/bowtie2.index/genome", "full bowtie2 index prefix path");
  }
  else if (op=="filter") {
    cmd_line->AddOption("-E", &ENZYME_REG_FILE, "", "enzyme fragments (BED/GFF/SAM/REG)");
    cmd_line->AddOption("--filter-self", &FILTER_SELF_FRAG, false, "filter read pairs mapping to same fragment");
    cmd_line->AddOption("--min-dist", &FILTER_MIN_DISTANCE, 500, "miminum allowed distance between 5's of reads in read pair");
    cmd_line->AddOption("--max-offset", &FILTER_MAX_OFFSET, 500, "maximum allowed offset of 5's of reads from fragment ends");
    cmd_line->AddOption("--filter-dups", &FILTER_DUPLICATES, false, "filter duplicate read pairs as PCR artifacts");
    cmd_line->AddOption("--class-file", &FILTER_CLASS_FILE, "", "output file for read pair classification details");
  }
  else if (op=="bin") {
    cmd_line->AddOption("--bin-size", &BIN_SIZE, 1000000, "genomic bin size");
  }
  else {
    cerr << "Unknown operation '" << op << "'!\n";
    delete cmd_line;
    exit(1);
  }

  // process command line
  *next_arg = cmd_line->Read(argv+1,argc-1) + 1;
  if (HELP||(argc-*next_arg<min_files[op])) { cmd_line->OperationUsage(); delete cmd_line; exit(1); }
  
  return cmd_line;
}




//---MakeDirectory------
//
bool MakeDirectory(char *dir_name)
{
  struct stat st;
  if (stat(dir_name,&st)!=0) {
    if (VERBOSE) fprintf(stderr, "* Making directory '%s'\n", dir_name);
    mkdir(dir_name,0777);
    return true;
  }
  return false;
}



//---CLASS RegionPair------
//
class RegionPair
{
 public:
  RegionPair(FileBuffer *buffer, GenomicRegionSetIndex *refIndex);
  ~RegionPair();

  GenomicRegionSAM *r[2];
  GenomicRegion *r_matched[2];
  long int start_offset[2], stop_offset[2], min_offset[2];
  int n_frags[2];
  PairType p_type; 
};

//---constructor-------
//
RegionPair::RegionPair(FileBuffer *buffer, GenomicRegionSetIndex *refIndex)
{
  // read #1
  r[0] = new GenomicRegionSAM(buffer);
  
  // check if second read-id is present
  bool found;
  if (buffer->Get()!=NULL) {  
    char *p = StrCopy(buffer->Get());
    char *pp = p;
    char *read_id = GetNextToken(&pp,'\t');
    found = strcmp(r[0]->LABEL,read_id)==0;
    free(p);
  }
  else found = false; 
  if (found==false) {
    p_type = UNPAIRED;
    r[1] = NULL;
    return;
  }

  // read #2
  r[1] = new GenomicRegionSAM(buffer);
  
  // classify according to mappability
  if ((strcmp(r[0]->I[0]->CHROMOSOME,"*")==0)&&(strcmp(r[1]->I[0]->CHROMOSOME,"*")==0)) p_type = UNMAPPED;
  else if ((strstr(r[0]->OPTIONAL,"XS:i:")!=NULL)&&(strstr(r[1]->OPTIONAL,"XS:i:")!=NULL)) p_type = MULTIHIT;
  else if ((strcmp(r[0]->I[0]->CHROMOSOME,"*")!=0)&&(strcmp(r[1]->I[0]->CHROMOSOME,"*")!=0)&&(strstr(r[0]->OPTIONAL,"XS:i:")==NULL)&&(strstr(r[1]->OPTIONAL,"XS:i:")==NULL)) {
    // assign read pair to corresponding enzyme fragments
    for (int q=0; q<=1; q++) {
      n_frags[q] = 0;
      start_offset[q] = stop_offset[q] = min_offset[q] = -1;
      if (r[q]->I.size()!=1) r[q]->PrintError("alignment region is not contiguous!\n"); 
      GenomicInterval read5p = r[q]->I.front();
      read5p.ModifyPos("5p",0);       // assigned fragment will be the one containing the 5p of the read; therefore it will be unique
      r_matched[q] = refIndex->GetMatch(&read5p);
      if (r_matched[q]!=NULL) {
        n_frags[q] = 1; 
        r_matched[q]->I.front()->GetOffsetFrom(&read5p,"5p",true,&start_offset[q],&stop_offset[q]); 
        min_offset[q] = min(abs(start_offset[q]),abs(stop_offset[q]));
      }
    }
    // classify according to fragment assignment
    if ((n_frags[0]==0)||(n_frags[1]==0)) p_type = NO_FRAG; 
    else {
      long int read_dist = abs(r[0]->I.front()->GetCoordinate("5p")-r[1]->I.front()->GetCoordinate("5p"));
      if (r_matched[0]==r_matched[1]) p_type = FRAG_SELF;
      else if ((r[0]->I.front()->STRAND!=r[1]->I.front()->STRAND)&&(read_dist<FILTER_MIN_DISTANCE)) p_type = DIST_SMALL;
      else if ((min_offset[0]<=FILTER_MAX_OFFSET)&&(min_offset[1]<=FILTER_MAX_OFFSET)) p_type = FRAG_PROXIMAL;
      else if ((min_offset[0]>FILTER_MAX_OFFSET)||(min_offset[1]>FILTER_MAX_OFFSET)) p_type = FRAG_DISTAL;
    }
  }
  else if ((strcmp(r[0]->I[0]->CHROMOSOME,"*")!=0)||(strcmp(r[1]->I[0]->CHROMOSOME,"*")!=0)) p_type = SINGLE;
  else p_type = UNMAPPED_MULTIHIT;
}

//---destructor-------
RegionPair::~RegionPair()
{
  if (r[0]) delete r[0];
  if (r[1]) delete r[1];
}

//---ENDCLASS RegionPair-------



//------IsDuplicate------
//
bool IsDuplicate(GenomicInterval *r1, GenomicInterval *r2)
{
  string chr_pair = (string)r1->CHROMOSOME + r1->STRAND + (string)r2->CHROMOSOME + r2->STRAND;
  MapOfRead::iterator i_chr = dupMask.find(chr_pair);
  if (i_chr==dupMask.end()) {
    MapOfRead2 *m2 = dupMask[chr_pair] = new MapOfRead2();
    MapOfRead1 *m1 = (*m2)[r2->START] = new MapOfRead1();
    (*m1)[r1->START] = true;
    return false;
  }
  MapOfRead2 *m2 = i_chr->second;
  MapOfRead2::iterator im2 = m2->find(r2->START);
  if (im2==m2->end()) {
    MapOfRead1 *m1 = (*m2)[r2->START] = new MapOfRead1();
    (*m1)[r1->START] = true;
    return false;
  }
  MapOfRead1 *m1 = im2->second;
  MapOfRead1::iterator im1 = m1->find(r1->START);
  if (im1==m1->end()) {
    (*m1)[r1->START] = true;
    return false;
  }
  return true;
}




//------RunFilter------
//
void RunFilter(char **args, int argn, char *ref_reg_file)
{
  // process args
  char *query_reg_file = argn==0?NULL:args[0];

  // check
  if (strlen(ref_reg_file)==0) { fprintf(stderr, "Error: need to supply enzyme fragment region file!\n"); exit(1); }
  
  // open region sets
  FileBuffer *QueryRegBuffer = CreateFileBuffer(query_reg_file,BUFFER_SIZE);
  GenomicRegionSet RefRegSet(ref_reg_file,BUFFER_SIZE,VERBOSE,true,true);
  GenomicRegionSetIndex RefIndex(&RefRegSet,BIN_BITS);

  // process read pairs
  Progress PRG("Processing aligned reads...",1);
  unsigned long int n_reads = 0;
  unsigned long int n_class[n_pair_types];
  FILE *F = strlen(FILTER_CLASS_FILE)>0?fopen(FILTER_CLASS_FILE,"w"):NULL;
  for (int i=0; i<n_pair_types; i++) n_class[i] = 0;
  for (char *qstr=QueryRegBuffer->Next(); QueryRegBuffer->Get()!=NULL; ) {
    if (qstr[0]=='@') { QueryRegBuffer->Next(); continue; }    // skip header lines
    RegionPair p(QueryRegBuffer,&RefIndex);
    if (p.p_type!=UNPAIRED) {
      if (p.p_type==FRAG_PROXIMAL) {
        GenomicInterval read1_5p = p.r[0]->I.front();
        read1_5p.ModifyPos("5p",0);
        GenomicInterval read2_5p = p.r[1]->I.front();
        read2_5p.ModifyPos("5p",0);
        if ((FILTER_DUPLICATES==false)||(IsDuplicate(&read1_5p,&read2_5p)==false)) { 
          printf("%s\t", p.r[0]->LABEL);
          read1_5p.PrintInterval();
          printf(" ");
          read2_5p.PrintInterval();
          printf("\n");
        }
      }
      // details to be saved in file
      if (F!=NULL) {
        fprintf(F, "%s\t%s\t%s\t", p.r[0]->LABEL, p.r[0]->SEQ, p.r[1]->SEQ);
        p.r[0]->I.front()->PrintInterval(F);
        fprintf(F, "\t");
        p.r[1]->I.front()->PrintInterval(F);
        for (int q=0; q<=1; q++) {
          if (p.r_matched[q]!=NULL) fprintf(F, "\t%s:%ld,%ld,%ld", p.r_matched[q]->LABEL, p.start_offset[q], p.stop_offset[q], p.min_offset[q]);
          else fprintf(F, "\t*******");
        }
        fprintf(F, "\t%s\n", PairTypeStr[p.p_type]);
      }
    }
    n_reads++;
    n_class[p.p_type]++;
    PRG.Check();
  }
  if (F!=NULL) fclose(F);
  PRG.Done();

  // report
  fprintf(stderr, "reads = %lu\n", n_reads);
  for (int i=0; i<n_pair_types; i++) {
    fprintf(stderr, "%s =", PairTypeStr[i]);
    fprintf(stderr, " %lu (%.2f%%)", n_class[i], (double)100*n_class[i]/n_reads);
    fprintf(stderr, "\n");
  }
  
  // cleanup
  delete QueryRegBuffer;
  for (MapOfRead::iterator i_chr=dupMask.begin(); i_chr!=dupMask.end(); i_chr++) {
    MapOfRead2 *m2 = i_chr->second;
    for (MapOfRead2::iterator i_m2=m2->begin(); i_m2!=m2->end(); i_m2++) delete i_m2->second;
    delete m2;
  }
}




//------RunAlign------
//
void RunAlign(char **args, int argn, char *work_dir, int min_len, int len_diff, int n_threads, char *bowtie_path, char *bowtie_index)
{
  // process args
  char *fastqFile1 = args[0];
  char *fastqFile2 = args[1];

  // determine read length
  FileBuffer *buffer = CreateFileBuffer(fastqFile1,BUFFER_SIZE); 
  buffer->Next();  // skip one line
  int read_len = strlen(buffer->Next());
  if (VERBOSE) fprintf(stderr, "* read length is %dnt\n", read_len);
  delete buffer;
  
  // initialize
  if (strlen(work_dir)==0) { fprintf(stderr, "Error: working directory must be specified!\n"); exit(1); }
  if (MakeDirectory(work_dir)==false) { fprintf(stderr, "Error: directory '%s' already exists!\n", work_dir); exit(1); }
  string bowtie_params = "--very-sensitive --reorder --sam-nohead -M 1";
  string fastq_file[2];
  string sam_file[2];
  fastq_file[0] = fastqFile1;
  fastq_file[1] = fastqFile2;
  for (int len=read_len; len>=min_len; len-=len_diff) {
    if (VERBOSE) fprintf(stderr, "Aligning reads (length=%d)...\n", len);
    for (int k=0; k<=1; k++) {
      // align
      stringstream ss_sam;
      ss_sam << work_dir << "/read" << k+1 << ".len=" << len << ".sam";
      sam_file[k] = ss_sam.str();
      stringstream cmd;
      cmd << bowtie_path << " -p " << n_threads << " " << bowtie_params << " -x " << bowtie_index << " -U " << fastq_file[k] << " -S " << sam_file[k];
      cerr << cmd.str() << '\n';
      system(cmd.str().c_str());
      // create new fastq files
      if (len-len_diff<min_len) continue;
      stringstream ss_q;
      ss_q << work_dir << "/read" << k+1 << ".len=" << len-len_diff << ".fastq";
      fastq_file[k] = ss_q.str();
      FileBuffer *buffer = CreateFileBuffer(sam_file[k].c_str(),BUFFER_SIZE); 
      FILE *Fout = fopen(fastq_file[k].c_str(),"w"); 
      Progress PRG("Processing SAM file...",1);
      for (buffer->Next(); buffer->Get()!=NULL; ) {
        GenomicRegionSAM r(buffer);
        GenomicInterval *i = r.I.front();
        if (strcmp(i->CHROMOSOME,"*")==0) fprintf(Fout, "@%s\n%.*s\n+%s\n%.*s\n", r.LABEL, len-len_diff, r.SEQ, r.LABEL, len-len_diff, r.QUAL);
        PRG.Check();
      }
      delete buffer;
      fclose(Fout);
      PRG.Done();
    }
  }

  // merge output sam files and print
  int n = (read_len-min_len)/len_diff+1;
  FileBuffer *samBuffer[2][n];
  for (int k=0; k<=1; k++) {
    for (int len=read_len,i=0; len>=min_len; len-=len_diff,i++) {
      stringstream ss_sam;
      ss_sam << work_dir << "/read" << k+1 << ".len=" << len << ".sam";
      samBuffer[k][i] = CreateFileBuffer(ss_sam.str().c_str(),BUFFER_SIZE); 
      samBuffer[k][i]->Next();
    }
  }
  Progress PRG2("Merging SAM files...",1);
  bool done = false;
  while (done==false) {
    for (int k=0; k<2; k++) {
      for (int i=0; i<n; i++) {
        char *s = samBuffer[k][i]->Get();
        if (s==NULL) { done = true; break; }
        GenomicRegionSAM r(s);
        samBuffer[k][i]->Next(); 
        if ((i==n-1)||(strcmp(r.I.front()->CHROMOSOME,"*")!=0)) { r.Print(); break; }
      }
    }
    PRG2.Check();
  }
  for (int k=0; k<=1; k++) for (int len=read_len,i=0; len>=min_len; len-=len_diff,i++) delete samBuffer[k][i];
  PRG2.Done();
}



//------RunBin------
//
void RunBin(char **args, int argn, long int bin_size)
{
  // process args
  char *filteredRegFile = argn==0?NULL:args[0];

  // open region set
  FileBuffer *filteredRegBuffer = CreateFileBuffer(filteredRegFile,BUFFER_SIZE);

  // process read pairs
  Progress PRG("Binning filtered reads...",1);
  while (filteredRegBuffer->Next()!=NULL) {
    GenomicRegion r(filteredRegBuffer);
    if (r.I.size()!=2) { r.PrintError("each line should be a pair of aligned reads intervals"); exit(1); }
    printf("%s:%ld\t%s:%ld\n", r.I[0]->CHROMOSOME, (r.I[0]->START/bin_size)+1, r.I[1]->CHROMOSOME, (r.I[1]->START/bin_size)+1);
    PRG.Check();
  }
  PRG.Done();


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

  // run 
  if (cmd_line->current_cmd_operation=="align") RunAlign(&argv[next_arg],argc-next_arg,ALIGN_WORK_DIR,ALIGN_MIN_LEN,ALIGN_LEN_DIFF,BOWTIE_THREADS,BOWTIE_PATH,BOWTIE_INDEX);
  else if (cmd_line->current_cmd_operation=="filter") RunFilter(&argv[next_arg],argc-next_arg,ENZYME_REG_FILE);
  else if (cmd_line->current_cmd_operation=="bin") RunBin(&argv[next_arg],argc-next_arg,BIN_SIZE);
  else { cerr << "Unknown method '" << cmd_line->current_cmd_operation << "'!\n"; delete cmd_line; exit(1); }

  // clean up
  delete cmd_line;

  return 0;
}



