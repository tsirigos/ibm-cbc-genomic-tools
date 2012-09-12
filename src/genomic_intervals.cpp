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
#include <iostream>
#include <sstream>
#include <map>
#include <string>
#include <limits>
#include <list>
#include <vector>
#include <iostream>
#include <algorithm>
#include "core.h"
#include "genomic_intervals.h"


using namespace std;
typedef map<string,string> StringMap;
typedef map<string,long int> StringLIntMap;





//---------------------------------------------------------------------------------------------//
// CLASS: Chromosomes                                                                          //
//---------------------------------------------------------------------------------------------//

//---------Constructor-----------
//
Chromosomes::Chromosomes(char *chrom_map_dir, char *chrom_map_name)
{
  // read sequence map file
  this->verbose = false;
  this->load_in_memory = false;
  this->chrom_map_dir = chrom_map_dir;
  this->chrom_map_name = chrom_map_name;
  string chrom_map_file = (string)chrom_map_dir + "/" + (string)chrom_map_name;
  FileBufferText buffer(chrom_map_file.c_str());
  Progress PRG("Reading chromosome map file...",1);
  int c = 0;
  for (; buffer.Next(); c++) {
    char *inp = buffer.Get();
    string chr = GetNextToken(&inp," \t");
    string file = GetNextToken(&inp," \t\n");
    chrom_map[chr] = file;
    PRG.Check();
  }
  PRG.Done();
  current_chromosome_name = "";
  current_chromosome_size = 0;
  current_chromosome_seq = NULL;
}



//---------Constructor-----------
//
Chromosomes::Chromosomes(char *fasta_file, bool verbose)
{
  // read sequence map file
  this->verbose = verbose;
  this->load_in_memory = true;
  FileBufferText buffer(fasta_file);
  Progress PRG("Reading chromosome sequences from FASTA file...",1);
  long int n_line = 0;
  for (char *inp=buffer.Next(); inp!=NULL; inp=buffer.Next()) {
    n_line++;
    if (inp[0]!='>') { fprintf(stderr, "[Chromosomes] Error: Line %ld: chromosome file '%s' is not in FASTA format!\n", n_line, fasta_file); exit(1); }
    ++inp;
    string chromosome_name = GetNextToken(&inp," \t");
    if (chromosome_seq_data.find(chromosome_name)!=chromosome_seq_data.end()) { fprintf(stderr, "[Chromosomes] Error: Line %ld: duplicate chromosome name in '%s'!\n", n_line, fasta_file); exit(1); }
    inp = buffer.Next();
    n_line++;
    if ((inp==NULL)||(inp[0]=='>')) { fprintf(stderr, "[Chromosomes] Error: Line %ld: sequence missing in chromosome file '%s'!\n", n_line, fasta_file); exit(1); }
    chromosome_seq_data[chromosome_name] = pair<char*,size_t>(inp,strlen(inp));
    buffer.BUFFER = new char[buffer.BUFFER_SIZE];									// using this trick to avoid duplicating inp
    //if (verbose) cerr << "* Read chromosome " << chromosome_name << " of size " << chromosome_seq_data[chromosome_name].second << ".\n";
    PRG.Check();
  }
  PRG.Done();
  current_chromosome_name = "";
  current_chromosome_size = 0;
  current_chromosome_seq = NULL;
}



//---------Destructor-----------
//
Chromosomes::~Chromosomes()
{
  if (load_in_memory==false) { if (current_chromosome_seq!=NULL) delete current_chromosome_seq; }
  else for (ChromosomeSeqData::iterator it=chromosome_seq_data.begin(); it!=chromosome_seq_data.end(); it++) delete it->second.first;
}



//------FindChromosome-----
//
bool Chromosomes::FindChromosome(char *chromosome_name)
{
  if (load_in_memory==true) return chromosome_seq_data.find(chromosome_name)!=chromosome_seq_data.end();
  else return chrom_map.find(chromosome_name)!=chrom_map.end();
}



//------LoadChromosome-----
//
void Chromosomes::LoadChromosome(string chromosome_name) 
{
  if (chromosome_name==current_chromosome_name) return;
  
  else if (load_in_memory==true) {
    if (chromosome_seq_data.find(chromosome_name)==chromosome_seq_data.end()) { cerr << "[Chromosomes] Error: chromosome '" << chromosome_name << "' not found in FASTA file!\n"; exit(1); }
    current_chromosome_name = chromosome_name;
    current_chromosome_seq = chromosome_seq_data[chromosome_name].first;
    current_chromosome_size = chromosome_seq_data[chromosome_name].second;
    return;
  }

  else if (chrom_map.find(chromosome_name)!=chrom_map.end()) {
    current_chromosome_name = chromosome_name;
    if (current_chromosome_seq!=NULL) delete current_chromosome_seq;
    string chromosome_file = chrom_map_dir + "/" + chrom_map[chromosome_name];

    FILE *F = fopen(chromosome_file.c_str(),"r");
    if (F==0) { cerr << "[Chromosomes] Error: can't open file '" << chromosome_file << "'!\n"; exit(1); }
    fseek(F,0,SEEK_END);
    long int file_size = ftell(F);
    if (file_size<0) { cerr << "[Chromosomes] Error reading file '" << chromosome_file << "'!\n"; exit(1); }
    fclose(F);
	current_chromosome_seq = new char[file_size+1];
	current_chromosome_seq[0] = 0;
    FileBufferText buffer(chromosome_file.c_str(),10000);
	char *inp = buffer.Next();
    if (inp==NULL) { cerr << "[Chromosomes] Error: chromosome file '" << chromosome_file << "' is empty!\n"; exit(1); }
	if (inp[0]=='>') inp = buffer.Next();
	char *s = current_chromosome_seq;
	while (inp!=NULL) {
	  strcpy(s,inp);
	  s += strlen(inp);
	  inp = buffer.Next();
	}
    current_chromosome_size = strlen(current_chromosome_seq);
  }

  else { cerr << "[Chromosomes] Error: chromosome '" << chromosome_name << "' not found in map file!\n"; exit(1); }

}




//---------GetSeq-----------
//
char *Chromosomes::GetSeq(GenomicInterval *I, bool replace)
{
  LoadChromosome(I->CHROMOSOME);  
  size_t n = I->GetSize();
  char *q = new char[n+1];
  q[n] = 0;
  if (I->STRAND=='+') {
    char *p = current_chromosome_seq + I->START - 1;
    for (size_t k=0; (p[0]!=0)&&(k<n); k++,p++) q[k] = replace&&(p[0]=='N') ? 'a' : p[0];
  }
  else if (I->STRAND=='-')  {
    char *p = current_chromosome_seq + I->STOP - 1;
    for (size_t k=0; (p[0]!=0)&&(k<n); k++,p--) q[k] = replace&&(p[0]=='N') ? 'a' : complement(p[0]);
  }
  else q[0] = 0;
  return q;
}



//---------PrintSeq-----------
//
void Chromosomes::PrintSeq(GenomicInterval *I, bool replace)
{
  LoadChromosome(I->CHROMOSOME);
  size_t n = I->GetSize();
  if ((size_t)I->STOP>current_chromosome_size) { cerr << "[Chromosomes] Error: interval [" << I->CHROMOSOME << ' ' <<  I->STRAND << ' ' << I->START << ' ' << I->STOP << "] exceeds chromosome bounds!\n"; exit(1); }
  if (I->STRAND=='+') {
    char *p = current_chromosome_seq + I->START - 1;
    for (size_t k=0; (p[0]!=0)&&(k<n); k++,p++) printf("%c", replace&&(p[0]=='N') ? 'a' : p[0]);
  }
  else if (I->STRAND=='-') {
    char *p = current_chromosome_seq + I->STOP - 1;
    for (size_t k=0; (p[0]!=0)&&(k<n); k++,p--) printf("%c", replace&&(p[0]=='N') ? 'a' : complement(p[0]));
  }
}



//---------GetSeqLength-----------
//
size_t Chromosomes::GetSeqLength(GenomicInterval *I)
{
  LoadChromosome(I->CHROMOSOME);
  size_t n = I->GetSize();
  size_t len = 0;
  char *p = current_chromosome_seq + I->START - 1;
  for (size_t k=0; (p[0]!=0)&&(k<n); k++,p++) len += p[0]!='N';
  return len;
}




//---------------------------------------------------------------------------------------------//
// END CLASS: Chromosomes                                                                      //
//---------------------------------------------------------------------------------------------//
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//







//---------------------------------------------------------------------------------------------//
// CLASS IMPLEMENTATION: GenomicInterval                                                       //
//---------------------------------------------------------------------------------------------//
//                                                                                             //
// Input format: CHROMOSOME STRAND START STOP                                                  //
//                                                                                             //
//---------------------------------------------------------------------------------------------//





//---------Constructor-----------
//
GenomicInterval::GenomicInterval(GenomicInterval *I)
{
  n_line = I->n_line;
  CHROMOSOME = StrCopy(I->CHROMOSOME);
  STRAND = I->STRAND;
  START = I->START;
  STOP = I->STOP;
  //CheckValid();
}



//---------Constructor-----------
//
GenomicInterval::GenomicInterval(char *chromosome, char strand, long int start, long int stop, unsigned long int n_line)
{
  this->n_line = n_line;
  CHROMOSOME = StrCopy(chromosome);
  STRAND = strand;
  START = start;
  STOP = stop;
  //CheckValid();
}



//---------Constructor-----------
//
GenomicInterval::GenomicInterval(string chromosome, char strand, long int start, long int stop, unsigned long int n_line)
{
  this->n_line = n_line;
  CHROMOSOME = StrCopy(chromosome.c_str());
  STRAND = strand;
  START = start;
  STOP = stop;
  //CheckValid();
}



//---------Constructor-----------
//
GenomicInterval::GenomicInterval(char *inp, char delimiter, unsigned long int n_line)
{
  this->n_line = n_line;
  if (CountTokens(inp,delimiter)!=4) { fprintf(stderr, "Line %ld: expected 4 space-separated tokens!\n", n_line); exit(1); }
  CHROMOSOME = StrCopy(GetNextToken(&inp,delimiter));
  STRAND = ProcessStrand(GetNextToken(&inp,delimiter));
  START = atol(GetNextToken(&inp,delimiter));
  STOP = atol(GetNextToken(&inp,delimiter));
  //CheckValid();
}



//---------Destructor-----------
//
GenomicInterval::~GenomicInterval()
{
  if (CHROMOSOME!=NULL) delete CHROMOSOME;
}



//---------SetCoordinates-----------
//
void GenomicInterval::SetCoordinates(long int start, long int stop)
{
  START = start;
  STOP = stop;
}

  

//---------IsValid-----------
//
bool GenomicInterval::IsValid()
{
  return (START<1)||(START>STOP);
}



//---------CheckValid-----------
//
bool GenomicInterval::CheckValid(bool quiet)
{
  if (START<1) { if (quiet==false) fprintf(stderr, "Line %ld: start/stop coordinates are 1-based\n", n_line); return false; }
  if (START>STOP) { if (quiet==false) fprintf(stderr, "Line %ld: start coordinate cannot exceed stop coordinate\n", n_line); return false; }
  return true;
}



//---------PrintInterval-----------
//
void GenomicInterval::PrintInterval()
{
  printf("%s %c %ld %ld", CHROMOSOME, STRAND, START, STOP);
}



//---------PrintInterval-----------
//
void GenomicInterval::PrintInterval(FILE *file_ptr)
{
  fprintf(file_ptr, "%s %c %ld %ld", CHROMOSOME, STRAND, START, STOP);
}



//---------PrintInterval-----------
//
void GenomicInterval::PrintInterval(long int start, long int stop)
{
  printf("%s %c %ld %ld", CHROMOSOME, STRAND, max(START,start), min(STOP,stop));
}



//---------Print-----------
//
void GenomicInterval::Print(bool point)
{
  if (point) for (long int p=START; p<=STOP; p++) printf("%s %c %ld %ld", CHROMOSOME, STRAND, p, p);
  else PrintInterval();
  printf("\n");
}



//---------PrintWithLabel-----------
//
void GenomicInterval::PrintWithLabel(char *label, bool point)
{
  if (point) for (long int p=START; p<=STOP; p++) printf("%s\t%s %c %ld %ld", label, CHROMOSOME, STRAND, p, p);
  else {
    printf("%s\t", label);
    PrintInterval();
  }
  printf("\n"); 
}



//---------IsBefore-----------
//
bool GenomicInterval::IsBefore(GenomicInterval *I, bool sorted_by_strand)
{
  if (strcmp(CHROMOSOME,I->CHROMOSOME)!=0) return strcmp(CHROMOSOME,I->CHROMOSOME)<0;
  if (sorted_by_strand&&(STRAND!=I->STRAND)) return STRAND<I->STRAND;
  return START<I->START;
}



//---------IsPosBefore-----------
//
bool GenomicInterval::IsPosBefore(GenomicInterval *I)
{
  return START<I->START;
}



//---------IsCompatibleWith-----------
//
bool GenomicInterval::IsCompatibleWith(GenomicInterval *I, bool ignore_strand)
{
  if (strcmp(CHROMOSOME,I->CHROMOSOME)!=0) return false;
  if ((ignore_strand==false)&&(STRAND!=I->STRAND)) return false;
  return true;
}



//---------CalcOverlap-----------
//
long int GenomicInterval::CalcOverlap(GenomicInterval *I, bool ignore_strand)
{
  if (IsCompatibleWith(I,ignore_strand)==false) return 0;
  long int y = min(STOP,I->STOP) - max(START,I->START) + 1;
  return y>0?y:0;
}



//---------CalcDistanceFrom----------
//
long int GenomicInterval::CalcDistanceFrom(GenomicInterval *I, const char *op, const char *I_op)
{
  return GetCoordinate(op) - I->GetCoordinate(I_op);
}




//---------CalcDirection-----------
//
int GenomicInterval::CalcDirection(GenomicInterval *i, bool sorted_by_strand)
{
  int cmp1 = strcmp(CHROMOSOME,i->CHROMOSOME);
  if (cmp1!=0) return cmp1;
  if (sorted_by_strand) {
    int cmp2 = STRAND - i->STRAND;
    if (cmp2!=0) return cmp2;
  }
  if (i->STOP<START) return 1;
  if (STOP<i->START) return -1;
  return 0;
}



//---------GetCoordinate-----------
//
long int GenomicInterval::GetCoordinate(const char *op)
{
  if (strcmp(op,"1")==0) return START;
  else if (strcmp(op,"2")==0) return STOP;
  else if (strcmp(op,"5p")==0) return STRAND=='+'?START:STOP;
  else if (strcmp(op,"3p")==0) return STRAND=='+'?STOP:START;
  else { cerr << "Error: unknown offset reference point operation!\n"; exit(1); }
}



//---------ReverseStrand-----------
//
void GenomicInterval::ReverseStrand()
{
  if (STRAND=='+') STRAND = '-';
  else if (STRAND=='-') STRAND = '+';
}



//---------GetSize-----------
//
size_t GenomicInterval::GetSize()
{
  return START>STOP ? 0 : STOP-START+1;
}


//---------GetSeq-----------
//
char *GenomicInterval::GetSeq(Chromosomes *C, bool replace)
{
  return C->GetSeq(this,replace);
}


//---------PrintSeq-----------
//
void GenomicInterval::PrintSeq(Chromosomes *C, bool replace)
{
  printf(">");
  PrintInterval();
  printf("\n");
  C->PrintSeq(this,replace);
  printf("\n");
}


//---------GetSeqLength--------
//
size_t GenomicInterval::GetSeqLength(Chromosomes *C)
{
  return C->GetSeqLength(this);
}


//---------ShiftPos-----------
//
void GenomicInterval::ShiftPos(long int start_shift, long int stop_shift, bool strand_aware)
{
  if (strand_aware==false) { START += start_shift; STOP += stop_shift; }
  else {
    if (STRAND=='-') { START -= stop_shift; STOP -= start_shift; }
    else { START += start_shift; STOP += stop_shift; }
  }
}



//---------ModifyPos-----------
//
void GenomicInterval::ModifyPos(const char *position_op, long int position_shift)
{
  if (strcmp(position_op,"c")==0) {
    long int center = START + (STOP-START+1)/2;
    START = center - position_shift;
    STOP = center + position_shift;
  }
  else {
    bool choose_start = (strcmp(position_op,"1")==0) || ((STRAND=='+')&&(strcmp(position_op,"5p")==0)) || ((STRAND=='-')&&(strcmp(position_op,"3p")==0));
    long int shift = choose_start ? position_shift : -position_shift;
    if (choose_start) STOP = START + shift;
    else START = STOP + shift;
  }
}



//---------ApplyBounds-----------
//
bool GenomicInterval::ApplyBounds(StringLIntMap *bounds)
{
  if ((bounds==NULL)||(bounds->find(CHROMOSOME)==bounds->end())) return false;
  size_t seq_size = (*bounds)[CHROMOSOME];
  if ((START>(long int)seq_size)||(STOP<1)||(STOP<START)) { 
    START = STOP = 0;
    return true;
  }
  else if ((START<1)||(STOP>(long int)seq_size)) { 
    START = max(START,1L);
    STOP = min(STOP,(*bounds)[CHROMOSOME]);
    return true;
  }
  return false;
}



//---------PrintBEDFormat-----------
//
void GenomicInterval::PrintBEDFormat(char *color, bool convert_chromosome)
{
  long int SCORE = 1000;
  PrintChromosome(CHROMOSOME,convert_chromosome);
  printf("%c%ld%c%ld%c_%c%ld%c%c", BED_SEPARATOR, START-1, BED_SEPARATOR, STOP, BED_SEPARATOR, BED_SEPARATOR, SCORE, BED_SEPARATOR, STRAND);
  printf("%c%ld%c%ld%c%s%c1%c%ld%c0\n", BED_SEPARATOR, START-1, BED_SEPARATOR, STOP, BED_SEPARATOR, color, BED_SEPARATOR, BED_SEPARATOR, STOP-START+1, BED_SEPARATOR);
}


//---------Randomize-----------
//
void GenomicInterval::Randomize(gsl_rng *random_generator, StringLIntMap *bounds)
{
  if (bounds->find(CHROMOSOME)==bounds->end()) { cerr << "Line " << n_line << ": chromosome " << CHROMOSOME << " not found!\n"; return; }
  long int d = STOP - START;
  unsigned long int N = (unsigned long int)((*bounds)[CHROMOSOME]-d);
  START = 1 + gsl_rng_uniform_int(random_generator,N);
  STOP = START + d;
}


//---------ReversePos-----------
//
void GenomicInterval::ReversePos(StringLIntMap *bounds)
{
  if (STRAND=='-') {
    long int start = (*bounds)[CHROMOSOME]-STOP+1;
    long int stop = (*bounds)[CHROMOSOME]-START+1;
    START = start;
    STOP = stop;
  }
}



//---------IsContained-----------
//
bool GenomicInterval::IsContained(GenomicInterval *I, bool ignore_strand)
{
  if (strcmp(CHROMOSOME,I->CHROMOSOME)!=0) return false;
  if ((ignore_strand==false)&&(STRAND!=I->STRAND)) return false;
  return (START>=I->START)&&(STOP<=I->STOP);
}



//---------OverlapsWith-----------
//
bool GenomicInterval::OverlapsWith(GenomicInterval *I, bool ignore_strand)
{
  if (strcmp(CHROMOSOME,I->CHROMOSOME)!=0) return false;
  if ((ignore_strand==false)&&(STRAND!=I->STRAND)) return false; 
  bool outside = (START>I->STOP) || (STOP<I->START);
  return !outside;
}



//---------OverlapsWith-----------
//
bool GenomicInterval::OverlapsWith(GenomicRegion *r, bool ignore_strand)
{
  for (GenomicIntervalSet::iterator p=r->I.begin(); p!=r->I.end(); p++) if (OverlapsWith(*p,ignore_strand)==true) return true;
  return false;
}



//---------GetOffsetFrom-----------
//
void GenomicInterval::GetOffsetFrom(GenomicInterval *ref_int, const char *op, bool ignore_strand, long int *start_offset, long int *stop_offset)
{
  if (IsCompatibleWith(ref_int,ignore_strand)==false) { cerr << "Error: [GetOffsetFrom] intervals must have the same chromosome/strand for this operation!\n"; exit(1); }
  long int ref = ref_int->GetCoordinate(op);
  char strand = ignore_strand?ref_int->STRAND:STRAND;
  bool test = ((strand=='-')&&(strcmp(op,"5p")==0))||((strand=='+')&&(strcmp(op,"3p")==0));
  if (test==true) { *start_offset = ref-STOP; *stop_offset = ref-START; }
  else { *start_offset = START-ref; *stop_offset = STOP-ref; }
}



//---------GetOffsetFrom-----------
//
void GenomicInterval::GetOffsetFrom(GenomicRegion *ref_reg, const char *op, bool ignore_strand, long int *start_offset, long int *stop_offset)
{
  if (ref_reg->IsCompatibleSortedAndNonoverlapping()==false) { cerr << "Error: [GetOffsetFrom] reference intervals must be compatible, sorted and non-overlapping!\n"; exit(1); }
  char strand = ref_reg->I.front()->STRAND;
  bool is_back = (strcmp(op,"2")==0)||((strand=='-')&&(strcmp(op,"5p")==0))||((strand=='+')&&(strcmp(op,"3p")==0));
  GenomicInterval *ref_int = is_back==false?ref_reg->I.front():ref_reg->I.back();
  GetOffsetFrom(ref_int,op,ignore_strand,start_offset,stop_offset);
}



//---------CreateCenter-----------
//
GenomicInterval *GenomicInterval::CreateCenter(bool stop_equals_start)
{
  long int new_start = START+(STOP-START)/2;
  long int new_stop = stop_equals_start?STOP-(STOP-START)/2:new_start;
  GenomicInterval *i = new GenomicInterval(CHROMOSOME,STRAND,new_start,new_stop);
  return i;
}


//---------------------------------------------------------------------------------------------//
// END CLASS: GenomicInterval                                                                  //
//---------------------------------------------------------------------------------------------//
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//










//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//
//---------------------------------------------------------------------------------------------//
// CLASS: GenomicIntervalSetAsArray                                                            //
//---------------------------------------------------------------------------------------------//


//---------Constructor-----------
//
GenomicIntervalSetAsArray::GenomicIntervalSetAsArray()
{
  max_intervals = 10; 
  I = new GenomicInterval*[max_intervals]; 
  n_intervals = 0;
}



//---------Destructor-----------
//
GenomicIntervalSetAsArray::~GenomicIntervalSetAsArray()
{
  delete [] I;
}


//---------------------------------------------------------------------------------------------//
// CLASS: GenomicIntervalSetAsArray                                                            //
//---------------------------------------------------------------------------------------------//
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//






//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//
//---------------------------------------------------------------------------------------------//
// CLASS: GenomicRegion                                                                        //
//---------------------------------------------------------------------------------------------//
//                                                                                             //
// Input formats:                                                                              //
//  1. LABEL <TAB> CHROMOSOME STRAND START[,START]* STOP[,STOP]*                               //
//  2. LABEL <TAB> [CHROMOSOME STRAND START STOP]+                                             //
//                                                                                             //
//---------------------------------------------------------------------------------------------//


//---------Constructor-----------
//
GenomicRegion::GenomicRegion()
{

}



//---------Constructor-----------
//
GenomicRegion::GenomicRegion(FileBuffer *B) 
{
  char *next = B->Get();
  n_line = B->n_line;
  Read(next,n_line); 
  B->Next();
}



//---------Constructor-----------
//
GenomicRegion::GenomicRegion(char *label, GenomicInterval *i)
{
  n_line = 0;
  LABEL = StrCopy(label);
  I.push_back(i);
}



//---------Constructor-----------
//
GenomicRegion::GenomicRegion(char *inp, long int n_line)
{
  this->n_line = n_line;
  Read(inp,n_line); 
}




//---------Destructor-----------
//
GenomicRegion::~GenomicRegion()
{
  if (LABEL!=NULL) delete LABEL;
  for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) delete *i;
  I.clear();
}



//------------------------------------------------------//
//   Read & Print                                       //
//------------------------------------------------------//



//---------Read-----------
//
void GenomicRegion::Read(char *inp, long int n_line, char del1, char del2)
{
  LABEL = StrCopy(GetNextToken(&inp,'\t'));
  //if (strchr(LABEL,' ')!=NULL) PrintError("no spaces allowed in region label!"); 

  if (strchr(inp,del2)==NULL) {				// regular REG format
    int n_tokens = CountTokens(inp,del1);
    if ((n_tokens<4)||(n_tokens%4!=0)) PrintError("invalid number of tokens!");
    long int n_intervals = n_tokens/4;
    for (long int k=0; k<n_intervals; k++) {
      char *chromosome = GetNextToken(&inp,del1);
      char strand = ProcessStrand(GetNextToken(&inp,del1));
      long int start = atol(GetNextToken(&inp,del1));
      long int stop = atol(GetNextToken(&inp,del1));
      I.push_back(new GenomicInterval(chromosome,strand,start,stop,n_line)); 
    }
  }
  else {        							// compact REG format
    int n_tokens = CountTokens(inp,del1);
    if (n_tokens!=4) PrintError("invalid number of tokens in compact format!");
    char *chromosome = GetNextToken(&inp,del1);
    char strand = ProcessStrand(GetNextToken(&inp,del1));
    char *starts = GetNextToken(&inp,del1);
    char *stops = GetNextToken(&inp,del1);
    long int n_intervals = CountTokens(starts,del2);
    if (CountTokens(stops,del2)!=n_intervals) PrintError("number of starts/stops should be equal");
    for (long int k=0; k<n_intervals; k++) {
      long int start = atol(GetNextToken(&starts,del2));
      long int stop = atol(GetNextToken(&stops,del2));
      I.push_back(new GenomicInterval(chromosome,strand,start,stop,n_line));
    }
  }
  // NOTE: need to create a warnings list for invalid intervals
}


//---------Print-----------
//
void GenomicRegion::Print(FILE *file_ptr)
{
  if (I.size()==0) return;
  fprintf(file_ptr, "%s\t", LABEL);
  for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) { (*i)->PrintInterval(file_ptr); if ((*i)!=I.back()) fprintf(file_ptr, " "); }
  fprintf(file_ptr, "\n");
}



//---------PrintIntersection-----------
//
void GenomicRegion::PrintIntersection(GenomicRegion *r, bool ignore_strand, bool merge_labels)
{
  // REQUIREMENTS: both regions must be compatible, sorted and non-overlapping
  bool first = true;
  for (GenomicIntervalSet::iterator it=I.begin(),r_it=r->I.begin(); (it!=I.end())&&(r_it!=r->I.end()); ) { 
    int d = (*it)->CalcDirection(*r_it,!ignore_strand);
    if (d<0) it++;
    else if (d>0) r_it++;
    else {
      long int new_start = max((*it)->START,(*r_it)->START);
      long int new_stop = min((*it)->STOP,(*r_it)->STOP);
      if (first) { 
        if (merge_labels) printf("%s:%s\t", LABEL, r->LABEL);
        else printf("%s\t", LABEL); 
        first = false;
      }
      else printf(" ");
      printf("%s %c %ld %ld", (*it)->CHROMOSOME, (*it)->STRAND, new_start, new_stop);
      if ((*it)->STOP==(*r_it)->STOP) { it++; r_it++; }
      else if ((*it)->STOP<(*r_it)->STOP) it++;
      else if ((*r_it)->STOP<(*it)->STOP) r_it++;
    }	
  }
  if (first==false) printf("\n");
}



//---------PrintConstrained-----------
//
void GenomicRegion::PrintConstrained(GenomicRegion *r, bool merge_labels)
{
  if (IsCompatible(false)==false) PrintError("[GenomicRegion::PrintConstrained] this method is defined for compatible regions only!");
  bool first = true;
  for (GenomicIntervalSet::iterator it=I.begin(); it!=I.end(); it++) { 
    long int new_start = max((*it)->START,r->I.front()->START);
    long int new_stop = min((*it)->STOP,r->I.back()->STOP);
    if (new_start<=new_stop) {
      if (first) { 
        if (merge_labels) printf("%s:%s\t", LABEL, r->LABEL);
        else printf("%s\t", LABEL); 
        first = false;
      }
      else printf(" ");
      printf("%s %c %ld %ld", (*it)->CHROMOSOME, (*it)->STRAND, new_start, new_stop);
    }
  }
  if (first==false) printf("\n");
}



//---------PrintModified-----------
//
void GenomicRegion::PrintModified(char *label, long int start, long int stop)
{
  if (IsCompatible(false)==false) PrintError("intervals should have the same chromosome/strand for this operation!");
  printf("%s\t%s %c %ld %ld\n", label, I.front()->CHROMOSOME, I.front()->STRAND, start, stop); 
}



//---------PrintQ-----------
//
/*
void GenomicRegion::PrintQ(bool reverse)
{
  for (size_t z=0; z<Q.size(); z++) {
    if (L[z].size()>0) cout << L[z] << ':';
    if (reverse==true) cout << ReverseComplement(Q[z]);
    else cout << Q[z]; 
    if (z<Q.size()-1) cout << '|';
  }
  cout << '\n';
}
*/


//---------PrintOffsetFormat-----------
//
void GenomicRegion::PrintOffsetFormat(char *op, bool fraction)
{
  bool ignore_strand = false; 
  printf("%s\t", LABEL);
  if (IsCompatible(ignore_strand)==false) PrintError("region intervals must have the same chromosome/strand for this operation!");
  I.front()->PrintInterval();
  size_t size = I.front()->GetSize();
  printf(" @%s@", op); 
  long int ref = I.front()->GetCoordinate(op);
  GenomicIntervalSet::iterator i = I.begin();
  i++;
  for ( ; i!=I.end(); i++) {
    long int start = (*i)->START - ref;
    long int stop = (*i)->STOP - ref;
    if (I.front()->STRAND=='-') { long int s = start; start = -stop; stop = -s; }
    if (fraction) printf(" %f %f", (float)start/size, (float)stop/size);
    else printf(" %ld %ld", start, stop);
    if (start>stop) PrintError("[PrintOffsetFormat] this must be a bug!");
  }
  printf("\n");
}



//---------PrintREG-----------
//
void GenomicRegion::PrintREG(bool compact)
{
  if (compact) {
    if (IsCompatible(false)==false) PrintError("region intervals must have the same chromosome/strand to perform conversion into compact REG format!");
    printf("%s\t%s %c ", LABEL, I.front()->CHROMOSOME, I.front()->STRAND);
    for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) printf("%ld%c", (*i)->START, (*i)==I.back()?' ':',');
    for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) printf("%ld%s", (*i)->STOP, (*i)==I.back()?"":",");
    printf("\n");
  }
  else GenomicRegion::Print();
  
  
  // NOTE: this code is for the SEQ format
  /*
  if (point) {
    if (fasta==true) {
      for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) {
        for (long int s=0,p=(*i)->START; p<=(*i)->STOP; p++,s++) {
          if ((*i)->STRAND=='-') PrintError("This operation is not tested for negative strand yet!");        // SOS: this can be fixed!!
          cout << Q[0][s] << '\t' << (*i)->CHROMOSOME << ' ' << (*i)->STRAND << ' ' << p << ' ' << p << '\n';
    }
  }
  else if (one) {
    if (fasta==true) PrintError("FASTA format for this operation is not implemented yet!");        // SOS: this can be fixed!!
    for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) { printf("%s\t", LABEL); (*i)->PrintInterval(); printf("\n"); }
  }
  else {
    printf(">");
    printf("%s\t", LABEL);
    PrintIntervals(compact);
    printf("\n");
    PrintQ();
  }
  */
}



//---------PrintError-----------
//
void GenomicRegion::PrintError(string error_msg)
{
  fprintf(stderr,"\n");
  fprintf(stderr, "Error: Line %ld: %s\n", n_line, error_msg.c_str());
  exit(1);
}




//------------------------------------------------------//
//   Get & Set                                          //
//------------------------------------------------------//


//---------SetLabel-----------
//
void GenomicRegion::SetLabel(char *label)
{
  if (LABEL!=NULL) delete LABEL;
  LABEL = StrCopy(label);
}



//---------SetStrand-----------
//
void GenomicRegion::SetStrand(char strand)
{
  if (IsCompatible(false)==false) PrintError("region intervals must have the same chromosome/strand for this operation!");
  for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) (*i)->STRAND = strand;
}



//---------GetChromosome-----------
//
char *GenomicRegion::GetChromosome()
{
  return I.front()->CHROMOSOME;
}



//---------GetSize-----------
//
size_t GenomicRegion::GetSize(bool skip_gaps)
{
  if (skip_gaps==true) {
	size_t size = 0;
    for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) size += (*i)->GetSize();
    return size;
  }
  else return I.back()->STOP - I.front()->START + 1;
}



//---------GetSize-----------
//
size_t GenomicRegion::GetSeqLength(Chromosomes *C)
{
  size_t len = 0;
  for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) len += (*i)->GetSeqLength(C);
  return len;
}



//---------DeleteIntervals-----------
//
void GenomicRegion::DeleteIntervals(GenomicIntervalSet::iterator i, GenomicIntervalSet::iterator j)
{
  for (GenomicIntervalSet::iterator k=i; k!=j; k++) delete (*k);
  I.erase(i,j);
}










//------------------------------------------------------//
//   Check                                              //
//------------------------------------------------------//



//---------IsCompatible-----------
//
bool GenomicRegion::IsCompatible(bool ignore_strand)
{
  if (I.size()<=1) return true;
  for (GenomicIntervalSet::iterator i0=I.begin(),i=i0+1; i!=I.end(); i++) 
    if ((*i0)->IsCompatibleWith(*i,ignore_strand)==false) return false;
  return true;
}



//---------IsCompatibleWith-----------
//
bool GenomicRegion::IsCompatibleWith(GenomicRegion *r, bool ignore_strand)
{
  if (IsCompatible(ignore_strand)==false) return false;
  if (r->IsCompatible(ignore_strand)==false) return false;
  return I.front()->IsCompatibleWith(r->I.front(),ignore_strand);
}



//---------IsCompatibleSorted-----------
//
bool GenomicRegion::IsCompatibleSorted(bool ignore_strand)
{
  if (IsCompatible(ignore_strand)==false) return false;
  GenomicIntervalSet::iterator i = I.begin();
  i++;
  GenomicIntervalSet::iterator j = I.begin();
  for ( ; i!=I.end(); i++,j++) if ((*i)->START<(*j)->START) return false;
  return true;
}



//---------IsCompatibleSortedAndNonoverlapping-----------
//
bool GenomicRegion::IsCompatibleSortedAndNonoverlapping()
{
  if (IsCompatibleSorted(false)==false) return false;
  GenomicIntervalSet::iterator i = I.begin();
  i++;
  GenomicIntervalSet::iterator j = I.begin();
  for ( ; i!=I.end(); i++,j++) if ((*i)->START<=(*j)->STOP) return false;
  return true;
}



//---------OverlapsWith-----------
//
bool GenomicRegion::OverlapsWith(GenomicRegion *r, bool ignore_strand)
{
  for (GenomicIntervalSet::iterator i1=I.begin(); i1!=I.end(); i1++)
    for (GenomicIntervalSet::iterator i2=r->I.begin(); i2!=r->I.end(); i2++) if ((*i1)->OverlapsWith(*i2,ignore_strand)) return true;
  return false;
}



//---------IsBefore-----------
//
bool GenomicRegion::IsBefore(GenomicRegion *r, bool sorted_by_strand)
{
  return I.front()->IsBefore(r->I.front(),sorted_by_strand);
}



//---------IsPosBefore-----------
//
bool GenomicRegion::IsPosBefore(GenomicRegion *r)
{
  return I.front()->IsPosBefore(r->I.front());
}



//---------CalcOverlap-----------
//
long int GenomicRegion::CalcOverlap(GenomicRegion *r, bool ignore_strand)
{
  long int y = 0;
  for (GenomicIntervalSet::iterator i1=I.begin(); i1!=I.end(); i1++)
    for (GenomicIntervalSet::iterator i2=r->I.begin(); i2!=r->I.end(); i2++) y += (*i1)->CalcOverlap(*i2,ignore_strand);
  return y;
}



//---------CalcDirection-----------
//
int GenomicRegion::CalcDirection(GenomicInterval *i, bool sorted_by_strand)
{
  int cmp1 = strcmp(I.front()->CHROMOSOME,i->CHROMOSOME);
  if (cmp1!=0) return cmp1;
  if (sorted_by_strand) {
    int cmp2 = I.front()->STRAND - i->STRAND;
    if (cmp2!=0) return cmp2;
  }
  if (i->STOP<I.front()->START) return 1;
  if (I.back()->STOP<i->START) return -1;
  return 0;
}



//---------CalcDirection-----------
//
int GenomicRegion::CalcDirection(GenomicRegion *r, bool sorted_by_strand)
{
  int cmp1 = strcmp(I.front()->CHROMOSOME,r->I.front()->CHROMOSOME);
  if (cmp1!=0) return cmp1;
  if (sorted_by_strand) {
    int cmp2 = I.front()->STRAND - r->I.front()->STRAND;
    if (cmp2!=0) return cmp2;
  }
  if (r->I.back()->STOP<I.front()->START) return 1;
  if (I.back()->STOP<r->I.front()->START) return -1;
  return 0;
}



//------------------------------------------------------//
//   Line-based (horizontal) operations                 //
//------------------------------------------------------//


//---------RunAlign-----------
//
void GenomicRegion::RunAlign(Chromosomes *C)
{
  PrintError("input sequences are required for this operation, e.g. use SAM format!");
}


  
//---------PrintBEDFormat-----------
//
void GenomicRegion::PrintBEDFormat(char *color, bool convert_chromosome)
{
  if (IsCompatibleSortedAndNonoverlapping()==false) PrintError("region intervals must be compatible, sorted and non-overlapping for this operation!");
  long int SCORE = 1000;
  PrintChromosome(I.front()->CHROMOSOME,convert_chromosome);
  printf("%c%ld%c%ld%c%s%c%ld%c%c", BED_SEPARATOR, I.front()->START-1, BED_SEPARATOR, I.back()->STOP, BED_SEPARATOR, LABEL, BED_SEPARATOR, SCORE, BED_SEPARATOR, I.front()->STRAND);
  if ((strlen(color)>0)||(I.size()>1)) printf("%c%ld%c%ld%c%s", BED_SEPARATOR, I.front()->START-1, BED_SEPARATOR, I.back()->STOP, BED_SEPARATOR, strlen(color)>0?color:"0");
  if (I.size()>1) {
    printf("%c%lu%c", BED_SEPARATOR, (unsigned long int)I.size(), BED_SEPARATOR);
    for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) printf("%ld%s", (*i)->STOP-(*i)->START+1, (*i)==I.back()?"":",");
    printf("%c0", BED_SEPARATOR);
    GenomicIntervalSet::iterator i=I.begin(), j=I.begin();
    for (i++; i!=I.end(); i++,j++) {
      if ((*i)->START<=(*j)->STOP) { fprintf(stderr, "Line %lu: exons cannot overlap!\n", n_line); exit(1); }
      printf(",%ld", (*i)->START-I.front()->START);
    }
  }
  printf("\n");
}



//---------ApplyBounds-----------
//
bool GenomicRegion::ApplyBounds(StringLIntMap *bounds)
{
  bool invalid = false; 
  for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) 
    if ((*i)->ApplyBounds(bounds)==true) invalid = true;
  Fix();
  return invalid;
}



//---------Center-----------
//
void GenomicRegion::Center()
{
  for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) {
    long int start = (*i)->START+((*i)->STOP-(*i)->START)/2;
    long int stop = (*i)->STOP-((*i)->STOP-(*i)->START)/2;
    (*i)->START = start;
    (*i)->STOP = stop;
  }
}



//---------ModifyChromosomeNames-----------
//
void GenomicRegion::ModifyChromosomeNames(map<string,string> &chromosome_conversion_map)
{
  for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) {
    map<string,string>::iterator x = chromosome_conversion_map.find((*i)->CHROMOSOME);
	if (x!=chromosome_conversion_map.end()) (*i)->CHROMOSOME = StrCopy(x->second.c_str());
  }
}



//---------Connect-----------
//
void GenomicRegion::Connect()
{
  if (I.size()<=1) return;
  if (IsCompatible(false)==false) PrintError("region intervals must have the same chromosome/strand for this operation!"); 
  long int start = I.front()->START;
  long int stop = I.front()->STOP;
  GenomicIntervalSet::iterator i = I.begin();
  i++;
  for ( ; i!=I.end(); i++) {
    if ((*i)->START<start) start = (*i)->START;
    if ((*i)->STOP>stop) stop = (*i)->STOP;
  }
  I.front()->START = start;
  I.front()->STOP = stop;
  i = I.begin();
  i++;
  for ( ; i!=I.end(); i++) delete (*i);
  i = I.begin();
  i++;
  I.erase(i,I.end());
}



//---------Diff-----------
//
void GenomicRegion::Diff()
{
  if (I.size()==0) return;
  if (I.size()==1) {
    for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) delete (*i);
    I.clear();
    return;
  }
  if (IsCompatibleSortedAndNonoverlapping()==false) PrintError("region intervals must be compatible, sorted and non-overlapping for this operation!");
  GenomicIntervalSet::iterator i = I.begin();
  GenomicIntervalSet::iterator j = i;
  GenomicIntervalSet::iterator k = i;
  for (i++; i!=I.end(); i++,k++) 
    if ((*i)->START-1>=(*k)->STOP+1) { (*j)->START = (*k)->STOP+1; (*j)->STOP = (*i)->START-1; j++; }
  for (GenomicIntervalSet::iterator i=j; i!=I.end(); i++) delete (*i);
  I.erase(j,I.end());
}



//---------RunCalcDistances----------
//
void GenomicRegion::RunCalcDistances(char *op1, char *op2)
{
  if (I.size()<=1) return;     		// NOTE: maybe we want to do this: printf("%s\tNaN\n", LABEL);
  else {
    printf("%s\t", LABEL);
    GenomicIntervalSet::iterator i = I.begin();
    GenomicIntervalSet::iterator j = i;
    for (i++; i!=I.end(); i++,j++) 
	  if (strcmp((*i)->CHROMOSOME,(*j)->CHROMOSOME)==0) printf("%ld%c", (*i)->CalcDistanceFrom(*j,op2,op1), *i==I.back()?'\n':' ');
	  else printf("Inf%c", *i==I.back()?'\n':' ');
  }
}



//---------Divide-----------
//
void GenomicRegion::Divide()
{
  for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) {
    if ((*i)->START>=(*i)->STOP) continue;
    long int middle = (*i)->START + (*i)->GetSize()/2 - 1;
    i = I.insert(i,new GenomicInterval((*i)->CHROMOSOME,(*i)->STRAND,(*i)->START,middle,n_line));
    i++;
    (*i)->START = middle+1;
  }
}



//---------RunDivide-----------
//
void GenomicRegion::RunDivide()
{
  printf("%s\t", LABEL);
  for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) {
    if ((*i)->START>=(*i)->STOP) printf("%s %c %ld %ld", (*i)->CHROMOSOME, (*i)->STRAND, (*i)->START, (*i)->STOP); 
    else {
      long int middle = (*i)->START + (*i)->GetSize()/2 - 1;
      printf("%s %c %ld %ld", (*i)->CHROMOSOME, (*i)->STRAND, (*i)->START, middle);
      printf(" ");
      printf("%s %c %ld %ld", (*i)->CHROMOSOME, (*i)->STRAND, middle+1, (*i)->STOP);
    }
    if (*i!=I.back()) printf(" ");
  }
  printf("\n");
}



//---------Fix-----------
//
void GenomicRegion::Fix()
{
  for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); ) {
    if ((*i)->CheckValid(true)==false) { delete *i; i = I.erase(i); if (i==I.end()) break; }
    else i++;
  }
}



//---------Intersect-----------
//
void GenomicRegion::Intersect()
{
  if (IsCompatible(false)==false) PrintError("region intervals must have the same chromosome/strand for this operation!");
  if (I.size()<=1) return;
  long int start = I.front()->START;
  long int stop = I.front()->STOP;
  GenomicIntervalSet::iterator i = I.begin();
  for (i++; i!=I.end(); i++) {
    start = max(start,(*i)->START);
    stop = min(stop,(*i)->STOP);
  }
  i = I.begin();
  if (stop>=start) { (*i)->START = start; (*i)->STOP = stop; i++; }
  DeleteIntervals(i,I.end());
}



//---------RunSize-----------
//
void GenomicRegion::RunSize()
{
  printf("%s\t%lu\n", LABEL, (long unsigned int)GetSize(true));
}



//---------ReversePos-----------
//
void GenomicRegion::ReversePos(StringLIntMap *bounds)
{
  for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) (*i)->ReversePos(bounds);
}



//---------ModifyPos-----------
//
void GenomicRegion::ModifyPos(const char *position_op, long int position_shift)
{
  for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) (*i)->ModifyPos(position_op,position_shift);
}



//---------Randomize-----------
//
void GenomicRegion::Randomize(gsl_rng *random_generator, StringLIntMap *bounds)
{
  if (IsCompatibleSortedAndNonoverlapping()==false) PrintError("region intervals must be compatible, sorted and non-overlapping for this operation!");
  char *CHROMOSOME = I.front()->CHROMOSOME; 
  if (bounds->find(CHROMOSOME)==bounds->end()) PrintError("chromosome not found in genome region file!");		// NOTE: if warning, then use this code:  DeleteIntervals(I.begin(),I.end()); return;
  long int N = (*bounds)[CHROMOSOME]-(I.back()->STOP-I.front()->START);
  if (N<=0) PrintError("out of chromosomal bounds!");
  long int shift = gsl_rng_uniform_int(random_generator,(unsigned long int)N)-I.front()->START+1; 
  for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) { (*i)->START += shift; (*i)->STOP += shift; }
}



//---------Select-----------
//
void GenomicRegion::Select(bool first, bool last, bool from5p, bool from3p)
{
  if (I.size()==0) return; 
  if (IsCompatible(false)==false) PrintError("region intervals must have the same chromosome/strand for this operation!");
  char strand = I.front()->STRAND;
  bool pos_first = first||(from5p&&(strand=='+'))||(from3p&&(strand=='-'));
  bool pos_last = last||(from5p&&(strand=='-'))||(from3p&&(strand=='+'));
  if ((pos_first==false)&&(pos_last==false)) DeleteIntervals(I.begin(),I.end());
  else if (I.size()==1) return;
  else if ((pos_first==true)&&(pos_last==true)) {
    if (I.size()==2) return;
    Sort();
    GenomicIntervalSet::iterator i = I.begin()+1;
    (*i)->START = I.back()->START;
    (*i)->STOP = I.back()->STOP;
    DeleteIntervals(I.begin()+2,I.end());	
  }
  else if ((pos_first==true)&&(pos_last==false)) {
    Sort();
    DeleteIntervals(I.begin()+1,I.end());	
  }
  else if ((pos_first==false)&&(pos_last==true)) {
    Sort();
    GenomicIntervalSet::iterator i = I.begin();
    (*i)->START = I.back()->START;
    (*i)->STOP = I.back()->STOP;
    DeleteIntervals(I.begin()+1,I.end());	
  }
}



//---------ShiftPos-----------
//
void GenomicRegion::ShiftPos(long int start_shift, long int stop_shift, bool strand_aware)
{
  for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) (*i)->ShiftPos(start_shift,stop_shift,strand_aware);
}



//---------RunShuffle-----------
//
void GenomicRegion::RunShuffle(gsl_rng *random_generator, GenomicRegionSet *refReg, StringLIntMap *index, StringVecLIntMap *loc)
{
  if (I.size()!=1) PrintError("single-interval regions expected for this operation!\n");
  string v = (string)I.front()->CHROMOSOME + (I.front()->STRAND=='+'?"+":"-");
  long int l = I.front()->STOP - I.front()->START + 1;
  if (loc->find(v)!=loc->end()) {
    unsigned long int N = (unsigned long int)((*loc)[v].back()-l);
    long int d = gsl_rng_uniform_int(random_generator,N);
    vector<long int>::iterator x = lower_bound((*loc)[v].begin()+1,(*loc)[v].end(),d) - 1;
    long int rr = (*index)[v]+x-(*loc)[v].begin();
    if ((rr<0)||(rr>=refReg->n_regions)) { fprintf(stderr, "Error [GenomicRegion::RunShuffle]: this must be a bug!\n"); exit(1); }
    long int start = refReg->R[rr]->I.front()->START + d - *x;
    long int stop = start + l - 1;
    long int new_stop = min(stop,refReg->R[rr]->I.front()->STOP);
    if (new_stop>=start) PrintModified(LABEL,start,new_stop);
    while (stop>refReg->R[rr]->I.front()->STOP) {
      stop = refReg->R[rr+1]->I.front()->START + stop - refReg->R[rr]->I.front()->STOP - 1;
      rr++;
      long int new_stop = min(stop,refReg->R[rr]->I.front()->STOP);
      if (new_stop>=refReg->R[rr]->I.front()->START) PrintModified(LABEL,refReg->R[rr]->I.front()->START,new_stop);
    }
  }
  //else AddToWarnings("chromosome/strand combination not found in reference regions!");   // NOTE: warnings are not implemented yet
}



//---------Sort-----------
//
void GenomicRegion::Sort()
{
  if (I.size()<=1) return;
  SortGenomicIntervals(&I);
}



//---------RunSplit-----------
//
void GenomicRegion::RunSplit()
{
  for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) { printf("%s\t", LABEL); (*i)->PrintInterval(); printf("\n"); }
}



//---------ReverseStrand-----------
//
void GenomicRegion::ReverseStrand()
{
  for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) (*i)->ReverseStrand();
}



//---------ModifyStrand-----------
//
void GenomicRegion::ModifyStrand(char *strand_op)
{
  if (IsCompatible(false)==false) PrintError("region intervals must have the same chromosome/strand for this operation!");
  if (strcmp(strand_op,"r")==0) {
    for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++)
    	if ((*i)->STRAND=='+') (*i)->STRAND = '-';
    	else if ((*i)->STRAND=='-') (*i)->STRAND = '+';
  }
  else if (strcmp(strand_op,"+")==0) {
    for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) (*i)->STRAND = '+';
  }
  else if (strcmp(strand_op,"-")==0) {
    for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) (*i)->STRAND = '-';
  }

  
  /* NOTE: this code below is for GenomicRegionSEQ
  if (strcmp(strand_op,"r")==0) {
    for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++)
    	if ((*i)->STRAND=='+') (*i)->STRAND = '-';
    	else if ((*i)->STRAND=='-') (*i)->STRAND = '+';
    Print();
    if (fasta==true) PrintQ(I.front()->STRAND!='+');
  }
  else if (strcmp(strand_op,"+")==0) {
    for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) (*i)->STRAND = '+';
    Print();
    if (fasta==true) PrintQ(I.front()->STRAND!='+');
  }
  else if (strcmp(strand_op,"-")==0) {
    for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) (*i)->STRAND = '-';
    Print();
    if (fasta==true) PrintQ(I.front()->STRAND!='-');
  }
*/  

}



//---------Union-----------
//
void GenomicRegion::Union()
{
  if (I.size()<=1) return; 
  if (IsCompatible(false)==false) PrintError("region intervals must have the same chromosome/strand for this operation!");
  Sort();
  long int start = I.front()->START;
  long int stop = I.front()->STOP;
  GenomicIntervalSet::iterator i = I.begin(), j = i;
  for (i++; i!=I.end(); i++) {
    if (stop<(*i)->START) {
      (*j)->START = start;
      (*j)->STOP = stop; 
      j++;
      start = (*i)->START;
      stop = (*i)->STOP;
    }
    else stop = max(stop,(*i)->STOP);
  }
  (*j)->START = start;
  (*j)->STOP = stop; 
  j++;
  DeleteIntervals(j,I.end());
}



//---------RunUnion-----------
//
void GenomicRegion::RunUnion()
{
  if (IsCompatible(false)==false) PrintError("region intervals must have the same chromosome/strand for this operation!");
  Sort();
  printf("%s\t", LABEL);
  long int start = I.front()->START;
  long int stop = I.front()->STOP;
  GenomicIntervalSet::iterator i=I.begin();
  for (i++; i!=I.end(); i++) {
    if (stop<(*i)->START) {
      printf("%s %c %ld %ld ", I.front()->CHROMOSOME, I.front()->STRAND, start, stop);
      start = (*i)->START;
      stop = (*i)->STOP;
    }
    else stop = max(stop,(*i)->STOP);
  }
  printf("%s %c %ld %ld\n", I.front()->CHROMOSOME, I.front()->STRAND, start, stop);
}



//---------PrintWindows-----------
//
void GenomicRegion::PrintWindows(long int win_step, long int win_size)
{
  if (I.size()!=1) PrintError("this operation requires single-interval regions!\n");
  for (long int z=I.front()->START,k=1; z+win_size-1<=I.front()->STOP; z+=win_step,k++) {
    printf("%s#%ld\t%s %c %ld %ld\n", LABEL, k, I.front()->CHROMOSOME, I.front()->STRAND, z, z+win_size-1);
  }
}



//---------PrintRemoveN--------
//
void GenomicRegion::PrintRemoveN(Chromosomes *C)
{
  if (IsCompatible(false)==false) PrintError("region intervals must have the same chromosome/strand for this operation!");
  C->LoadChromosome(I.front()->CHROMOSOME);
  for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) {
    char *p = C->current_chromosome_seq + (*i)->START - 1;
    long int start = (*i)->START;
    while (start<=(*i)->STOP) {
      while ((start<=(*i)->STOP)&&(p[0]=='N')) { start++; p++; }
      if (start>(*i)->STOP) break;
      long int stop = start;
      while ((stop<=(*i)->STOP)&&(p[0]!='N')) { stop++; p++; }
      printf("%s\t%s %c %ld %ld\n", LABEL, (*i)->CHROMOSOME, (*i)->STRAND, start, stop-1);
      start = stop;
    }
  }
}



//---------PrintSeqLength-----------
//
void GenomicRegion::PrintSeqLength(Chromosomes *C)
{
  printf("%s\t%lu\n", LABEL, (long unsigned int)GetSeqLength(C));
}



//---------GetSeq-----------
//
char *GenomicRegion::GetSeq(Chromosomes *C, bool replace)
{
  if (C->FindChromosome(I.front()->CHROMOSOME)==false) return NULL;
  size_t size = this->GetSize(true);
  char *seq = new char[size+1];
  seq[size] = 0;
  char *q = seq;
  if (I.front()->STRAND=='+') {							// NOTE: this assumes that intervals are compatible
    for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) { 
      char *s = (*i)->GetSeq(C,replace);
      size_t n = (*i)->GetSize(); 
      strncpy(q,s,n);
      q += n;
      delete s;
    }
  }
  else if (I.front()->STRAND=='-') {
    for (GenomicIntervalSet::reverse_iterator i=I.rbegin(); i!=I.rend(); i++) { 
      char *s = (*i)->GetSeq(C,replace);
      size_t n = (*i)->GetSize(); 
      strncpy(q,s,n);
      q += n;
      delete s;
    }
  }
  else seq[0] = 0;
  return seq;
}



//---------PrintSeq-----------
//
void GenomicRegion::PrintSeq(Chromosomes *C, bool replace)
{
  if (IsCompatible(false)==false) PrintError("region intervals must have the same chromosome/strand for this operation!");
  printf(">");
  Print();
  if (I.front()->STRAND=='+') {	
    for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) C->PrintSeq((*i),replace);
  }
  else if (I.front()->STRAND=='-') {
    for (GenomicIntervalSet::reverse_iterator i=I.rbegin(); i!=I.rend(); i++) C->PrintSeq((*i),replace);
  }
  printf("\n");
}



//---------PrintVerifySeq--------
//
void GenomicRegion::PrintVerifySeq(Chromosomes *C, bool ignore)
{
  char *seq = GetSeq(C,false);
  if (strcmp(LABEL,seq)==0) { printf("+\t"); Print(); }
  else {
    ReverseGene(seq,strlen(seq));
    if (strcmp(LABEL,seq)==0) { printf("-\t"); ReverseStrand(); Print(); }
    else if (ignore==false) { printf("0\t"); ReverseGene(seq,strlen(seq)); Print(); }
  }	
  free(seq);
}



//---------PrintSearch-----------
//
void GenomicRegion::PrintSearch(char *pattern, bool header, bool summary)
{
  fprintf(stderr, "Line %ld: this operation requires SEQ format as input!\n", n_line); exit(1);
}



/*
//---------PrintNoGaps-----------
//
void GenomicRegion::PrintNoGaps()
{
  printf(">");
  Print();
  for (size_t z=0; z<Q.size(); z++) {
    if (L[z].size()>0) cout << L[z] << ':';
    for (size_t a=0; a<Q[z].size(); a++) if (Q[z][a]!='-') cout << Q[z][a];
    if (z<Q.size()-1) cout << '|';
  }
  printf("\n");
}
*/






//------------------------------------------------------//
//   Operations returning new regions                   //
//------------------------------------------------------//

  
//---------Constrain-----------
//
GenomicRegion *GenomicRegion::Constrain(GenomicRegion *r, char *label)
{
  GenomicRegion *new_r = NULL;
  for (GenomicIntervalSet::iterator it=I.begin(); it!=I.end(); it++) {
    long int new_start = max((*it)->START,r->I.front()->START);
    long int new_stop = min((*it)->STOP,r->I.back()->STOP);
    if (new_start<=new_stop) { 
      GenomicInterval *i = new GenomicInterval(I.front()->CHROMOSOME,I.front()->STRAND,new_start,new_stop);
      if (new_r==NULL) new_r = new GenomicRegion(label==NULL?LABEL:label,i);
      else new_r->I.push_back(i);
    }
  }
  return new_r;
}



//---------Diff-----------
//
GenomicRegion *GenomicRegion::Diff(GenomicRegion *r)
{
  GenomicInterval *i = new GenomicInterval(I.front()->CHROMOSOME,I.front()->STRAND,r->I.front()->STOP+1,I.front()->START-1);
  GenomicRegion *rr = new GenomicRegion((char*)"_",i);
  return rr;
}



//---------------------------------------------------------------------------------------------//
// END CLASS: GenomicRegion                                                                    //
//---------------------------------------------------------------------------------------------//
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//














//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//
//---------------------------------------------------------------------------------------------//
// CLASS: GenomicRegionSEQ                                                                     //
//---------------------------------------------------------------------------------------------//
//                                                                                             //
// Input formats:                                                                              //
//  1. LABEL <TAB> CHROMOSOME STRAND START[,START]* STOP[,STOP]* <NEW-LINE> SEQUENCE           //
//  2. LABEL <TAB> [CHROMOSOME STRAND START STOP]+  <NEW-LINE> SEQUENCE                        //
//                                                                                             //
//---------------------------------------------------------------------------------------------//


//---------Constructor-----------
//
GenomicRegionSEQ::GenomicRegionSEQ(FileBuffer *B) 
{
  char *next = B->Get();
  n_line = B->n_line;
  if (next[0]!='>') PrintError("sequence header should start with '>'!\n");
  ++next;
  Read(next,n_line); 
  SEQ = StrCopy(B->Next());
  n_line = B->n_line;
  if (SEQ[0]=='>') PrintError("sequence should not start with '>'!\n");
  B->Next();
}


//---------Destructor-----------
//
GenomicRegionSEQ::~GenomicRegionSEQ()
{
  if (LABEL!=NULL) delete LABEL;
  for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) delete *i;
  I.clear();
  delete SEQ;
}


//---------PrintSearch-----------
//
void GenomicRegionSEQ::PrintSearch(char *pattern, bool header, bool summary)
{
  if (I.size()!=1) { fprintf(stderr, "Line %ld: this operation requires single-interval regions!\n", n_line); exit(1); }
  if (strlen(pattern)==0) { fprintf(stderr, "Line %ld: pattern is empty!\n", n_line); exit(1); }
  if (header) { printf(">"); Print(); }
  
  long int n = CountTokens(pattern,'|');
  char **P = new char*[n];
  for (long int k=0; k<n; k++) P[k] = StrCopy(GetNextToken(&pattern,'|'));
  size_t seq_len = strlen(SEQ);
  for (size_t pos=0; pos<seq_len; pos++) {
    bool match = false;
    size_t plen = 0;
    for (long int k=0; (k<n)&&(match==false); k++) {
      plen = strlen(P[k]);
      if (seq_len-pos<plen) continue;
      for (size_t j=0; j<plen; j++) 
        if (SEQ[pos+j]!=P[k][j]) { match = false; break; }
        else match = true;
    }
    if (match) {
      if (I.front()->STRAND=='+') printf("_\t%s %c %ld %ld\n", I.front()->CHROMOSOME, I.front()->STRAND, I.front()->START+pos, I.front()->START+pos+plen-1);
      else printf("_\t%s %c %ld %ld\n", I.front()->CHROMOSOME, I.front()->STRAND, I.front()->STOP-pos-plen+1, I.front()->STOP-pos);
    }
  }

  delete [] P;
}



/*
//---------PrintSearch(FULL OLD VERSION)-----------
//
void GenomicRegionSEQ::PrintSearch(char *pattern, bool header, bool summary)
{
  if (I.size()!=1) { fprintf(stderr, "Line %ld: this operation requires single-interval regions!\n", n_line); exit(1); }
  if (strlen(pattern)==0) { fprintf(stderr, "Line %ld: pattern is empty!\n", n_line); exit(1); }
  if (header) { printf(">"); Print(); }
  
  long int n = CountTokens(pattern,'|');
  char **P = new char*[n];
  for (long int k=0; k<n; k++) P[k] = StrCopy(GetNextToken(&pattern,'|'));
  for (size_t pos=0; pos<Q[0].size(); pos++) {
    bool match = false;
    size_t plen = 0;
    for (long int k=0; (k<n)&&(match==false); k++) {
      plen = strlen(P[k]);
      if (Q[0].size()-pos<plen) continue;
      for (size_t j=0; j<plen; j++) 
        if (Q[0][pos+j]!=P[k][j]) { match = false; break; }
        else match = true;
    }
    if (match) {
      if (I.front()->STRAND=='+') printf("_\t%s %c %ld %ld\n", I.front()->CHROMOSOME, I.front()->STRAND, I.front()->START+pos, I.front()->START+pos+plen-1);
      else printf("_\t%s %c %ld %ld\n", I.front()->CHROMOSOME, I.front()->STRAND, I.front()->STOP-pos-plen+1, I.front()->STOP-pos);
    }
  }

  delete [] P;
}
*/

//---------------------------------------------------------------------------------------------//
// END CLASS: GenomicRegionSEQ                                                                 //
//---------------------------------------------------------------------------------------------//
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//









  
  
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//
//---------------------------------------------------------------------------------------------//
// CLASS: GenomicRegionBEDToREG                                                                //
//---------------------------------------------------------------------------------------------//


//---------Constructor-----------
//
GenomicRegionBEDToREG::GenomicRegionBEDToREG(FileBuffer *B)
{
  char *next = B->Get();
  n_line = B->n_line;
  Read(next,n_line); 
  B->Next();
}



//---------Constructor-----------
//
GenomicRegionBEDToREG::GenomicRegionBEDToREG(char *inp, long int n_line)
{
  this->n_line = n_line;
  Read(inp,n_line); 
}



//---------Destructor-----------
//
GenomicRegionBEDToREG::~GenomicRegionBEDToREG()
{

}


//---------Read-----------
//
void GenomicRegionBEDToREG::Read(char *inp, long int n_line)
{
  char sep = strchr(inp,'\t')==NULL?' ':'\t';
  long int n_tokens = CountTokens(inp,sep);
  if (n_tokens<3) PrintError("number of tokens should be at least 3 for BED format!");
  char *chromosome = GetNextToken(&inp,sep);
  long int start = atol(GetNextToken(&inp,sep))+1;
  long int stop = atol(GetNextToken(&inp,sep));
  char strand = '+';
  LABEL = n_tokens==3?StrCopy("_"):StrCopy(GetNextToken(&inp,sep));
  if (n_tokens>=6) { 
    /*long int score =*/ atol(GetNextToken(&inp,sep));
    strand = ProcessStrand(GetNextToken(&inp,sep)); 
  }
  if (n_tokens!=12) {
    I.push_back(new GenomicInterval(chromosome,strand,start,stop,n_line));
  }
  else {
    /*long int thickStart =*/ atol(GetNextToken(&inp,sep));
    /*long int thickEnd =*/ atol(GetNextToken(&inp,sep));
    /*char *itemRgb =*/ GetNextToken(&inp,sep);
    long int n_intervals = atol(GetNextToken(&inp,sep));
    for (long int k=0; k<n_intervals; k++) I.push_back(new GenomicInterval(chromosome,strand,start,stop,n_line));
    char *blockSizes = GetNextToken(&inp,sep);
    for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) (*i)->STOP = atol(GetNextToken(&blockSizes,','));
    char *blockStarts = GetNextToken(&inp,sep);
    for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) { (*i)->START += atol(GetNextToken(&blockStarts,',')); (*i)->STOP += (*i)->START-1; }
  }
}
  
  



//---------------------------------------------------------------------------------------------//
// END CLASS: GenomicRegionBEDToREG                                                            //
//---------------------------------------------------------------------------------------------//
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//





  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

  
  
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//
//---------------------------------------------------------------------------------------------//
// CLASS: GenomicRegionBED                                                                     //
//---------------------------------------------------------------------------------------------//


//---------Constructor-----------
//
GenomicRegionBED::GenomicRegionBED(FileBuffer *B)
{
  char *next = B->Get();
  n_line = B->n_line;
  Read(next,n_line); 
  B->Next();
}



//---------Constructor-----------
//
GenomicRegionBED::GenomicRegionBED(char *inp, long int n_line)
{
  this->n_line = n_line;
  Read(inp,n_line); 
}



//---------Constructor-----------
//
GenomicRegionBED::GenomicRegionBED(char *label, GenomicInterval *i)
{
  n_tokens = 6;
  LABEL = StrCopy(label);
  score = 1000;
  I.push_back(i);
}



//---------Destructor-----------
//
GenomicRegionBED::~GenomicRegionBED()
{
  if (n_tokens==12) delete itemRgb;
}


//---------Read-----------
//
void GenomicRegionBED::Read(char *inp, long int n_line)
{
  char sep = strchr(inp,'\t')==NULL?' ':'\t';
  n_tokens = CountTokens(inp,sep);
  if (n_tokens<3) PrintError("number of tokens should be at least 3 for BED format!");
  char *chromosome = GetNextToken(&inp,sep);
  long int start = atol(GetNextToken(&inp,sep))+1;
  long int stop = atol(GetNextToken(&inp,sep));
  char strand = '+';
  LABEL = n_tokens==3?StrCopy("_"):StrCopy(GetNextToken(&inp,sep));
  if (n_tokens>=5) score = atol(GetNextToken(&inp,sep));
  if (n_tokens>=6) strand = ProcessStrand(GetNextToken(&inp,sep)); 
  if (n_tokens>=8) { thickStart = atol(GetNextToken(&inp,sep)); thickEnd = atol(GetNextToken(&inp,sep)); }
  if (n_tokens>=9) itemRgb = StrCopy(GetNextToken(&inp,sep));
  if (n_tokens!=12) {
    I.push_back(new GenomicInterval(chromosome,strand,start,stop,n_line));
  }
  else {
    long int n_intervals = atol(GetNextToken(&inp,sep));
    for (long int k=0; k<n_intervals; k++) I.push_back(new GenomicInterval(chromosome,strand,start,stop,n_line));
    char *blockSizes = GetNextToken(&inp,sep);
    for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) (*i)->STOP = atol(GetNextToken(&blockSizes,','));
    char *blockStarts = GetNextToken(&inp,sep);
    for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) { (*i)->START += atol(GetNextToken(&blockStarts,',')); (*i)->STOP += (*i)->START-1; }
  }
}
  
  

//---------Print-----------
//
void GenomicRegionBED::Print(FILE *file_ptr)
{
  if (I.size()==0) return;
  fprintf(file_ptr, "%s%c%ld%c%ld", I.front()->CHROMOSOME, BED_SEPARATOR, I.front()->START-1, BED_SEPARATOR, I.back()->STOP);
  if (n_tokens>=4) {
    fprintf(file_ptr, "%c%s", BED_SEPARATOR, LABEL);
    if (n_tokens>=5) {
      fprintf(file_ptr, "%c%ld", BED_SEPARATOR, score);
      if (n_tokens>=6) {
        fprintf(file_ptr, "%c%c", BED_SEPARATOR, I.front()->STRAND);
        if (n_tokens>=8) {
          fprintf(file_ptr, "%c%ld%c%ld", BED_SEPARATOR, thickStart, BED_SEPARATOR, thickEnd);
          if (n_tokens>=9) {
            fprintf(file_ptr, "%c%s", BED_SEPARATOR, itemRgb);
            if (n_tokens==12) {
              fprintf(file_ptr, "%c%lu%c", BED_SEPARATOR, (unsigned long int)I.size(), BED_SEPARATOR);
              for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) fprintf(file_ptr, "%ld%s", (*i)->STOP-(*i)->START+1, *i==I.back()?"":",");
              fprintf(file_ptr, "%c0", BED_SEPARATOR);
              GenomicIntervalSet::iterator i=I.begin(), j=I.begin();
              for (i++; i!=I.end(); i++,j++) {
                if ((*i)->START<=(*j)->STOP) { fprintf(stderr, "Line %lu: exons cannot overlap!\n", n_line); exit(1); }
                fprintf(file_ptr, ",%ld", (*i)->START-I.front()->START);
              }
            }
          }
        }
      }
    }
  }
  fprintf(file_ptr, "\n");
}



//---------Print-----------
//
void GenomicRegionBED::Print(FILE *file_ptr, char *color, bool convert_chromosome)
{
  if (I.size()==0) return;
  if (convert_chromosome) fprintf(file_ptr, "chr%s%c", strcmp(I.front()->CHROMOSOME,"MT")==0?"M":I.front()->CHROMOSOME, BED_SEPARATOR);
  else fprintf(file_ptr, "%s%c", I.front()->CHROMOSOME, BED_SEPARATOR);
  fprintf(file_ptr, "%ld%c%ld", I.front()->START-1, BED_SEPARATOR, I.back()->STOP);
  if (n_tokens>=4) {
    fprintf(file_ptr, "%c%s", BED_SEPARATOR, LABEL);
    if (n_tokens>=5) {
      fprintf(file_ptr, "%c%ld", BED_SEPARATOR, score);
      if (n_tokens>=6) {
        fprintf(file_ptr, "%c%c", BED_SEPARATOR, I.front()->STRAND);
        if (n_tokens>=8) {
          fprintf(file_ptr, "%c%ld%c%ld", BED_SEPARATOR, thickStart, BED_SEPARATOR, thickEnd);
          if (n_tokens>=9) {
            fprintf(file_ptr, "%c%s", BED_SEPARATOR, strcmp(color,"")==0?itemRgb:color);
            if (n_tokens==12) {
              fprintf(file_ptr, "%c%lu%c", BED_SEPARATOR, (unsigned long int)I.size(), BED_SEPARATOR);
              for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) fprintf(file_ptr, "%ld%s", (*i)->STOP-(*i)->START+1, *i==I.back()?"":",");
              fprintf(file_ptr, "%c0", BED_SEPARATOR);
              GenomicIntervalSet::iterator i=I.begin(), j=I.begin();
              for (i++; i!=I.end(); i++,j++) {
                if ((*i)->START<=(*j)->STOP) { fprintf(stderr, "Line %lu: exons cannot overlap!\n", n_line); exit(1); }
                fprintf(file_ptr, ",%ld", (*i)->START-I.front()->START);
              }
            }
          }
        }
      }
    }
  }
  fprintf(file_ptr, "\n");
}



//---------PrintIntersection-----------
//
void GenomicRegionBED::PrintIntersection(GenomicRegion *r, bool ignore_strand, bool merge_labels)
{
  // REQUIREMENTS: both regions must be compatible, sorted and non-overlapping
  bool first = true;
  for (GenomicIntervalSet::iterator it=I.begin(),r_it=r->I.begin(); (it!=I.end())&&(r_it!=r->I.end()); ) { 
    int d = (*it)->CalcDirection(*r_it,!ignore_strand);
    if (d<0) it++;
    else if (d>0) r_it++;
    else {
      long int new_start = max((*it)->START,(*r_it)->START);
      long int new_stop = min((*it)->STOP,(*r_it)->STOP);
      if (first) { 
        if (merge_labels) printf("%s:%s%c", LABEL, r->LABEL, BED_SEPARATOR);
        else printf("%s%c", LABEL, BED_SEPARATOR);
        first = false;
      }
      else printf("%c", BED_SEPARATOR);
      printf("%s%c%c%c%ld%c%ld", (*it)->CHROMOSOME, BED_SEPARATOR, (*it)->STRAND, BED_SEPARATOR, new_start, BED_SEPARATOR, new_stop);
      if ((*it)->STOP==(*r_it)->STOP) { it++; r_it++; }
      else if ((*it)->STOP<(*r_it)->STOP) it++;
      else if ((*r_it)->STOP<(*it)->STOP) r_it++;
    }	
  }
  if (first==false) printf("\n");
}



//---------PrintConstrained-----------
//
void GenomicRegionBED::PrintConstrained(GenomicRegion *r, bool merge_labels)
{
  GenomicRegion *new_r;
  if (merge_labels) {
    char *new_label = new char[strlen(LABEL)+1+strlen(r->LABEL)+1];
    sprintf(new_label, "%s:%s", LABEL, r->LABEL);
    new_r = Constrain(r,new_label);
    delete new_label;
  }
  else new_r = Constrain(r);
  new_r->Print();
  delete new_r;
}



//---------PrintModified-----------
//
void GenomicRegionBED::PrintModified(char *label, long int start, long int stop)
{
  if (IsCompatible(false)==false) PrintError("intervals should have the same chromosome/strand for this operation!");
  printf("%s%c%ld%c%ld%c%s%c%ld%c%c\n", I.front()->CHROMOSOME, BED_SEPARATOR, start-1, BED_SEPARATOR, stop, BED_SEPARATOR, label, BED_SEPARATOR, score, BED_SEPARATOR, I.front()->STRAND);
}



//------------------------------------------------------//
//   Line-based (horizontal) operations                 //
//------------------------------------------------------//



//---------UpdateThick-----------
//
void GenomicRegionBED::UpdateThick()
{
  if (I.size()>0) {
    thickStart = max(thickStart,I.front()->START-1);
    thickEnd = min(thickEnd,I.back()->STOP);
    if (thickEnd<=thickStart) thickStart = thickEnd = I.front()->START-1;
  }
}



//---------PrintBEDFormat-----------
//
void GenomicRegionBED::PrintBEDFormat(char *color, bool convert_chromosome)
{
  Print(stdout,color,convert_chromosome);
}



//---------ApplyBounds-----------
//
bool GenomicRegionBED::ApplyBounds(StringLIntMap *bounds)
{
  bool invalid = GenomicRegion::ApplyBounds(bounds);
  UpdateThick();
  return invalid;
}



//---------Center-----------
//
void GenomicRegionBED::Center()
{
  GenomicRegion::Center();
  UpdateThick();
}



//---------Connect-----------
//
void GenomicRegionBED::Connect()
{
  GenomicRegion::Connect();
  UpdateThick();
}



//---------Diff-----------
//
void GenomicRegionBED::Diff()
{
  GenomicRegion::Diff();
  UpdateThick();
}



//---------RunCalcDistances-----------
//
void GenomicRegionBED::RunCalcDistances(char *op1, char *op2)
{
  GenomicRegion::RunCalcDistances(op1,op2);
}



//---------Divide-----------
//
void GenomicRegionBED::Divide()
{
  GenomicRegion::Divide();
  UpdateThick();
}



//---------RunDivide-----------
//
void GenomicRegionBED::RunDivide()
{
  Divide();
  Print();
}



//---------Fix-----------
//
void GenomicRegionBED::Fix()
{
  GenomicRegion::Fix();
  UpdateThick();
}



//---------Intersect-----------
//
void GenomicRegionBED::Intersect()
{
  PrintError("the result of the 'intersect' operation is empty or trivial for BED inputs!");
}



//---------ModifyPos-----------
//
void GenomicRegionBED::ModifyPos(const char *position_op, long int position_shift)
{
  GenomicRegion::ModifyPos(position_op,position_shift);
  GenomicRegion::Union();
  //UpdateThick();
}



//---------Randomize-----------
//
void GenomicRegionBED::Randomize(gsl_rng *random_generator, StringLIntMap *bounds)
{
  long int old_start = I.front()->START;
  GenomicRegion::Randomize(random_generator,bounds);
  long int d = I.front()->START - old_start;
  thickStart += d;
  thickEnd += d;  
}



//---------Select-----------
//
void GenomicRegionBED::Select(bool first, bool last, bool from5p, bool from3p)
{
  GenomicRegion::Select(first,last,from5p,from3p);
  UpdateThick();
}




//---------ShiftPos-----------
//
void GenomicRegionBED::ShiftPos(long int start_shift, long int stop_shift, bool strand_aware)
{
  GenomicRegion::ShiftPos(start_shift,stop_shift,strand_aware);
  GenomicRegion::Union();
  UpdateThick();
}




//---------RunShuffle-----------
//
void GenomicRegionBED::RunShuffle(gsl_rng *random_generator, GenomicRegionSet *refReg, StringLIntMap *index, StringVecLIntMap *loc)
{
  GenomicRegion::RunShuffle(random_generator,refReg,index,loc);
}



//---------Sort-----------
//
void GenomicRegionBED::Sort()
{
  return; 		// no need to do anything because BED regions are by definition sorted
}


  
//---------RunSplit-----------
//
void GenomicRegionBED::RunSplit()
{
  for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) PrintModified(LABEL,(*i)->START,(*i)->STOP);
}



//---------Union-----------
//
void GenomicRegionBED::Union()
{
  return; 		// no need to do anything because BED regions are by definition non-overlapping
}



//---------PrintWindows-----------
//
void GenomicRegionBED::PrintWindows(long int win_step, long int win_size)
{
  if (I.size()!=1) PrintError("this operation requires single-interval regions!\n");
  for (long int z=I.front()->START,k=1; z+win_size-1<=I.front()->STOP; z+=win_step,k++) {
    printf("%s%c%ld%c%ld", I.front()->CHROMOSOME, BED_SEPARATOR, z-1, BED_SEPARATOR, z+win_size-1);
    if (n_tokens>=4) {
      printf("%c%s", BED_SEPARATOR, LABEL);
      if (n_tokens>=6) printf("%c%ld%c%c", BED_SEPARATOR, score, BED_SEPARATOR, I.front()->STRAND);
    }
    printf("\n");	
  }
}



//------------------------------------------------------//
//   Operations returning new regions                   //
//------------------------------------------------------//

  
//---------Constrain-----------
//
GenomicRegionBED *GenomicRegionBED::Constrain(GenomicRegion *r, char *label)
{
  GenomicRegionBED *new_r = NULL;
  for (GenomicIntervalSet::iterator it=I.begin(); it!=I.end(); it++) {
    long int new_start = max((*it)->START,r->I.front()->START);
    long int new_stop = min((*it)->STOP,r->I.back()->STOP);
    if (new_start<=new_stop) { 
      GenomicInterval *i = new GenomicInterval(I.front()->CHROMOSOME,I.front()->STRAND,new_start,new_stop);
      if (new_r==NULL) new_r = new GenomicRegionBED(label==NULL?LABEL:label,i);
      else new_r->I.push_back(i);
    }
  }
  new_r->n_tokens = n_tokens;
  new_r->score = n_tokens>=5?score:0;
  new_r->thickStart = n_tokens>=7?max(thickStart,r->I.front()->START-1):r->I.front()->START-1;
  new_r->thickEnd = n_tokens>=8?min(thickEnd,r->I.back()->STOP):r->I.back()->STOP;
  new_r->itemRgb = n_tokens>=9?StrCopy(itemRgb):NULL;
  return new_r;
}



//---------Diff-----------
//
GenomicRegion *GenomicRegionBED::Diff(GenomicRegion *r)
{
  GenomicInterval *i = new GenomicInterval(I.front()->CHROMOSOME,I.front()->STRAND,r->I.front()->STOP+1,I.front()->START-1,n_line);
  GenomicRegion *rr = new GenomicRegionBED((char*)"_",i);
  return rr;
}





//---------------------------------------------------------------------------------------------//
// END CLASS: GenomicRegionBED                                                                 //
//---------------------------------------------------------------------------------------------//
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//














  
  

//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//
//---------------------------------------------------------------------------------------------//
// CLASS: GenomicRegionSAMToREG                                                                //
//---------------------------------------------------------------------------------------------//


//---------Constructor-----------
//
GenomicRegionSAMToREG::GenomicRegionSAMToREG(FileBuffer *B)
{
  char *next = B->Get();
  n_line = B->n_line;
  Read(next,n_line); 
  B->Next();
}



//---------Constructor-----------
//
GenomicRegionSAMToREG::GenomicRegionSAMToREG(char *inp, long int n_line)
{
  this->n_line = n_line;
  Read(inp,n_line); 
}



//---------Destructor-----------
//
GenomicRegionSAMToREG::~GenomicRegionSAMToREG()
{

}



//---------Read-----------
//
void GenomicRegionSAMToREG::Read(char *inp, long int n_line)
{
  char sep = '\t';
  long int n_tokens = CountTokens(inp,sep);
  if (n_tokens<11) PrintError("number of tokens should be at least 11 for SAM format!");
  LABEL = StrCopy(GetNextToken(&inp,sep));
  unsigned long int FLAG = (unsigned long int)atol(GetNextToken(&inp,sep));
  char *RNAME = GetNextToken(&inp,sep);
  long int POS = atol(GetNextToken(&inp,sep));
  /*MAPQ =*/ atol(GetNextToken(&inp,sep));
  char *CIGAR = GetNextToken(&inp,sep);
  /*
  RNEXT = StrCopy(GetNextToken(&inp,sep));
  PNEXT = atol(GetNextToken(&inp,sep));
  TLEN = atol(GetNextToken(&inp,sep));
  SEQ = StrCopy(GetNextToken(&inp,sep));
  QUAL = StrCopy(GetNextToken(&inp,sep));
  OPTIONAL = n_tokens>11?StrCopy(GetNextToken(&inp,'\n')):NULL;
  */
  
  long int start = POS;
  char strand = CalcStrandFromFlag(FLAG);
  char *p = CIGAR;
  char token_type;
  long int reference_len = 0;
  for (long int token_len=GetNextTokenOfCIGAR(&p,&token_type); token_len!=-1; token_len=GetNextTokenOfCIGAR(&p,&token_type)) {
    if (token_type!='N') reference_len += token_len;
    else {
      I.push_back(new GenomicInterval(RNAME,strand,start,start+reference_len-1,n_line));
      start = start + reference_len + token_len;
      reference_len = 0;
    }
  }
  if (reference_len>0) I.push_back(new GenomicInterval(RNAME,strand,start,start+reference_len-1,n_line));
}
  
  

//---------CalcStrandFromFlag-----------
//
char GenomicRegionSAMToREG::CalcStrandFromFlag(unsigned long int flag)
{
  return (flag&0x0010)==0x0010?'-':'+';
}



//---------GetNextTokenOfCIGAR-----------
//
long int GenomicRegionSAMToREG::GetNextTokenOfCIGAR(char **cigar, char *type)
{
  if ((*cigar)[0]!=0) {
    char *p = *cigar;
    while (((*cigar)[0]>='0')&&((*cigar)[0]<='9')) (*cigar)++;
    *type = (*cigar)[0];
    char *token = new char[(*cigar)-p+1];
    size_t i = 0;
    for ( ; p!=(*cigar); i++,p++) token[i] = p[0];
    token[i] = 0;
    long int val = atol(token);
    delete token;
    (*cigar)++;
    return val;
  }
  else { *type = ' '; return -1; }
}




//---------------------------------------------------------------------------------------------//
// END CLASS: GenomicRegionSAMToREG                                                            //
//---------------------------------------------------------------------------------------------//
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//





  
  
  
  
  
   
  
  
  

  
  
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//
//---------------------------------------------------------------------------------------------//
// CLASS: GenomicRegionSAM                                                                     //
//---------------------------------------------------------------------------------------------//


//---------Constructor-----------
//
GenomicRegionSAM::GenomicRegionSAM(FileBuffer *B)
{
  char *next = B->Get();
  n_line = B->n_line;
  Read(next,n_line); 
  B->Next();
}



//---------Constructor-----------
//
GenomicRegionSAM::GenomicRegionSAM(char *inp, long int n_line)
{
  this->n_line = n_line;
  Read(inp,n_line); 
}



//---------Destructor-----------
//
GenomicRegionSAM::~GenomicRegionSAM()
{
  delete CIGAR;
  delete RNEXT;
  delete SEQ;
  delete QUAL;
  if (OPTIONAL!=NULL) delete OPTIONAL;
}


//---------Read-----------
//
void GenomicRegionSAM::Read(char *inp, long int n_line)
{
  char sep = '\t';
  n_tokens = CountTokens(inp,sep);
  if (n_tokens<11) PrintError("number of tokens should be at least 11 for SAM format!");
  LABEL = StrCopy(GetNextToken(&inp,sep));
  FLAG = (unsigned long int)atol(GetNextToken(&inp,sep));
  char *RNAME = GetNextToken(&inp,sep);
  long int POS = atol(GetNextToken(&inp,sep));
  MAPQ = atol(GetNextToken(&inp,sep));
  CIGAR = StrCopy(GetNextToken(&inp,sep));
  RNEXT = StrCopy(GetNextToken(&inp,sep));
  PNEXT = atol(GetNextToken(&inp,sep));
  TLEN = atol(GetNextToken(&inp,sep));
  SEQ = StrCopy(GetNextToken(&inp,sep));
  QUAL = StrCopy(GetNextToken(&inp,sep));
  OPTIONAL = n_tokens>11?StrCopy(GetNextToken(&inp,'\n')):NULL;

  long int start = POS;
  char strand = CalcStrandFromFlag(FLAG);
  if (strcmp(CIGAR,"*")==0) { stringstream ss; ss << strlen(SEQ) << "M"; delete CIGAR; CIGAR = StrCopy(ss.str().c_str()); }   // WARNING: Invalid CIGAR string!!!
  TokenizeCIGAR();
  if ((strcmp(SEQ,"*")!=0)&&((long int)strlen(SEQ)!=CalcFragmentLengthFromCIGAR())) {
    stringstream err_msg;
    err_msg << "length of aligned fragment does not match CIGAR string: \n"	
            << "  LABEL = " << LABEL << '\n' 
            << "  CIGAR = " << CIGAR << '\n' 
            << "  length(SEQ) = " << strlen(SEQ) << '\n';
    PrintError(err_msg.str());
  }
  long int reference_len = 0;
  for (CIGARTokens::iterator i=T.begin(); i!=T.end(); i++) {
    if (i->second!='N') { if (strchr(refCIGARop,i->second)!=NULL) reference_len += i->first; }
    else {
      I.push_back(new GenomicInterval(RNAME,strand,start,start+reference_len-1,n_line));
      start = start + reference_len + i->first;
      reference_len = 0;
    }
  }
  if (reference_len>0) I.push_back(new GenomicInterval(RNAME,strand,start,start+reference_len-1,n_line));


}
  
  

//---------Print-----------
//
void GenomicRegionSAM::Print(FILE *file_ptr)
{
  if (I.size()==0) return;
  fprintf(file_ptr, "%s\t%ld\t%s\t%ld\t%ld\t%s\t%s\t%ld\t%ld\t%s\t%s", LABEL, FLAG, I.front()->CHROMOSOME, I.front()->START, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL); 
  if (n_tokens>11) fprintf(file_ptr, "\t%s", OPTIONAL);
  fprintf(file_ptr, "\n");
}



//---------PrintConstrained-----------
//
void GenomicRegionSAM::PrintConstrained(GenomicRegion *r, bool merge_labels)
{
  long int new_start = max(I.front()->START,r->I.front()->START);
  long int d_start = r->I.front()->START - I.front()->START;
  long int d_stop = I.back()->STOP - r->I.back()->STOP;
  char *new_cigar = TrimCIGAR(d_start,d_stop);
  long int new_TLEN = 0;							// NOTE: this cannot be updated, so set it to 0
  char *new_seq = TrimSequence(SEQ,d_start,d_stop);
  char *new_qual = TrimSequence(QUAL,d_start,d_stop);
  if (merge_labels) printf("%s:%s", LABEL, r->LABEL);
  else printf("%s", LABEL);
  printf("\t%ld\t%s\t%ld\t%ld\t%s\t%s\t%ld\t%ld\t%s\t%s", FLAG, I.front()->CHROMOSOME, new_start, MAPQ, new_cigar, RNEXT, PNEXT, new_TLEN, new_seq, new_qual); 
  if (n_tokens>11) printf("\t%s", OPTIONAL);
  printf("\n");
  delete new_cigar;
  delete new_seq;
  delete new_qual;
}



//---------PrintModified-----------
//
void GenomicRegionSAM::PrintModified(char *label, long int start, long int stop)
{
  if (IsCompatible(false)==false) PrintError("intervals should have the same chromosome/strand for this operation!");
  printf("%s\t%ld\t%s\t%ld\t%d\t%ldM\t%s\t%d\t%d\t%s\t%s", label, FLAG, I.front()->CHROMOSOME, start, 255, stop-start+1, RNEXT, 0, 0, "*", "*"); 
  if (n_tokens>11) printf("\t%s", OPTIONAL);
  printf("\n");
}



//------------------------------------------------------//
//   Auxiliary functions                                //
//------------------------------------------------------//

//---------CalcStrandFromFlag-----------
//
char GenomicRegionSAM::CalcStrandFromFlag(unsigned long int flag)
{
  return (flag&0x0010)==0x0010?'-':'+';
}



//---------UpdateFlagFromStrand-----------
//
void GenomicRegionSAM::UpdateFlagFromStrand(char strand)
{
  FLAG = (FLAG&0xFFEF)|(strand=='-'?0x0010:0x0000); 
}



//---------GetNextTokenOfCIGAR-----------
//
long int GenomicRegionSAM::GetNextTokenOfCIGAR(char **cigar, char *type)
{
  if ((*cigar)[0]!=0) {
    char *p = *cigar;
    while (((*cigar)[0]>='0')&&((*cigar)[0]<='9')) (*cigar)++;
    *type = (*cigar)[0];
    char *token = new char[(*cigar)-p+1];
    size_t i = 0;
    for ( ; p!=(*cigar); i++,p++) token[i] = p[0];
    token[i] = 0;
    long int val = atol(token);
    delete token;
    (*cigar)++;
    return val;
  }
  else { *type = ' '; return -1; }
}



//---------TokenizeCIGAR-----------
//
void GenomicRegionSAM::TokenizeCIGAR()
{
  char S[] = "MIDNSHP-X";
  char *p = CIGAR;
  char token_type;
  for (long int token_len=GetNextTokenOfCIGAR(&p,&token_type); token_len!=-1; token_len=GetNextTokenOfCIGAR(&p,&token_type)) 
    if (strchr(S,token_type)!=NULL) T.push_back(pair<long int,char>(token_len,token_type));
    else PrintError((string)"unknown CIGAR operation type '" + token_type + "'!");
}



//---------PrintGappedSequence-----------
//
void GenomicRegionSAM::PrintGappedSequence(char *seq, char *gap_token_types, bool lowercase_for_sort_clipping)
{
  char *q = seq;
  for (CIGARTokens::iterator tok_it=T.begin(); tok_it!=T.end(); tok_it++) {
    if (tok_it->second=='H') continue;
    if (tok_it->second=='N') printf("--N=%ld--", tok_it->first); 
    else {
      for (long int k=0; k<tok_it->first; k++) 
        if (strchr(gap_token_types,tok_it->second)!=NULL) printf("-"); 
        else { printf("%c", (lowercase_for_sort_clipping&&(tok_it->second=='S'))?tolower(q[0]):q[0]); q++; }
    }
  }
}


  
//---------GetGappedSequence-----------
//
char *GenomicRegionSAM::GetGappedSequence(char *seq, char *gap_token_types, bool lowercase_for_sort_clipping)
{
  size_t n = 0;
  for (CIGARTokens::iterator tok_it=T.begin(); tok_it!=T.end(); tok_it++) if (tok_it->second!='H') n += tok_it->second=='N'?3:tok_it->first;
  char *gapped_seq = new char[n+1];
  char *p = gapped_seq;
  char *q = seq;
  for (CIGARTokens::iterator tok_it=T.begin(); tok_it!=T.end(); tok_it++) {
    if (tok_it->second=='H') continue;
    if (tok_it->second=='N') { p[0] = '.'; p[1] = '.'; p[2] = '.'; p += 3; }			// NOTE: this is non-standard; is there a standard approach to this??
    else {
      for (long int k=0; k<tok_it->first; k++,p++) 
        if (strchr(gap_token_types,tok_it->second)!=NULL) p[0] = ALIGNMENT_GAP;
        else { p[0] = (lowercase_for_sort_clipping&&(tok_it->second=='S'))?tolower(q[0]):q[0]; q++; }
    }
  }
  p[0] = 0;
  return gapped_seq;
}


  
//---------InvalidateSeqData-----------
//
void GenomicRegionSAM::InvalidateSeqData()
{
  MAPQ = 255;
  PNEXT = 0;
  TLEN = 0;
  long int reference_len = CalcReferenceLengthFromCIGAR();
  stringstream cigar_stream;
  cigar_stream << reference_len << "M";
  delete CIGAR;
  CIGAR = StrCopy(cigar_stream.str().c_str());
  T.clear();
  TokenizeCIGAR();
  delete SEQ;
  SEQ = StrCopy("*");
  delete QUAL;
  QUAL = StrCopy("*");
  if (OPTIONAL!=NULL) { delete OPTIONAL; OPTIONAL = NULL; n_tokens = 11; }
}



//---------UpdateCigarFromIntervals-----------
//
void GenomicRegionSAM::UpdateCigarFromIntervals()
{
  stringstream cigar_stream;
  for (GenomicIntervalSet::iterator i=I.begin(),j=i+1; i!=I.end(); i++,j++) {
    long int d = (*i)->GetSize();
    cigar_stream << d << "X";
    if (*i!=I.back()) { 
      long int d = (*j)->START - (*i)->STOP - 1;
      cigar_stream << d << "N";
    } 
  }
  delete CIGAR;
  CIGAR = StrCopy(cigar_stream.str().c_str());
  T.clear();
  TokenizeCIGAR();
}



//---------CalcReferenceLengthFromCIGAR-----------
//
long int GenomicRegionSAM::CalcReferenceLengthFromCIGAR()
{
  long int reference_len = 0;
  for (CIGARTokens::iterator i=T.begin(); i!=T.end(); i++) if (strchr(refCIGARop,i->second)!=NULL) reference_len += i->first;
  return reference_len;
}



//---------CalcFragmentLengthFromCIGAR-----------
//
long int GenomicRegionSAM::CalcFragmentLengthFromCIGAR()
{
  long int fragment_len = 0;
  for (CIGARTokens::iterator i=T.begin(); i!=T.end(); i++) if (strchr(fragCIGARop,i->second)!=NULL) fragment_len += i->first;
  return fragment_len;
}



//---------TrimCIGAR-----------
//
char *GenomicRegionSAM::TrimCIGAR(long int d_ref_start, long int d_ref_stop)
{
  // NOTE: 'N' CIGAR op cannot be present here (this method is valid only for single-interval regions, i.e. no split fragments allowed).
  CIGARTokens Tcopy = T;
  long int k_start = 0;
  if (d_ref_start>0) {
    for ( ; k_start<(long int)Tcopy.size(); k_start++) {
      if (strchr(refCIGARop,Tcopy[k_start].second)==NULL) continue;
      if (d_ref_start>=Tcopy[k_start].first) { d_ref_start -= Tcopy[k_start].first; Tcopy[k_start].first = 0; }
      else { Tcopy[k_start].first -= d_ref_start; break; }
    }
  }
  long int k_stop = Tcopy.size()-1;
  if (d_ref_stop>0) {
    for ( ; k_stop>=0; k_stop--) {
      if (strchr(refCIGARop,Tcopy[k_stop].second)==NULL) continue;
      if (d_ref_stop>=Tcopy[k_stop].first) { d_ref_stop -= Tcopy[k_stop].first; Tcopy[k_stop].first = 0; }
      else { Tcopy[k_stop].first -= d_ref_stop; break; }
    }
  }
  stringstream cigar_stream; 
  for (long int k=k_start; k<=k_stop; k++) cigar_stream << Tcopy[k].first << Tcopy[k].second;
  return StrCopy(cigar_stream.str().c_str());
}



//---------TrimSequence-----------
//
char *GenomicRegionSAM::TrimSequence(char *seq, long int d_ref_start, long int d_ref_stop)
{
  if (strcmp(seq,"*")==0) return StrCopy(seq);
  long int fragment_start = 0;
  if (d_ref_start>0) {
    for (long int k=0; k<(long int)T.size(); k++) {
      if (strchr(fragCIGARop,T[k].second)!=NULL) fragment_start += T[k].second=='I'?T[k].first:min(d_ref_start,T[k].first);
      if (strchr(refCIGARop,T[k].second)!=NULL) { if (d_ref_start<T[k].first) break; else d_ref_start -= T[k].first; }
    }
  }
  long int fragment_stop = strlen(seq)-1;
  if (d_ref_stop>0) {
    for (long int k=T.size()-1; k>=0; k--) {
      if (strchr(fragCIGARop,T[k].second)!=NULL) fragment_stop -= T[k].second=='I'?T[k].first:min(d_ref_stop,T[k].first);
      if (strchr(refCIGARop,T[k].second)!=NULL) { if (d_ref_stop<T[k].first) break; else d_ref_stop -= T[k].first; }
    }
  }
  long int n = fragment_stop-fragment_start+1;
  if (n<=0) return StrCopy("*");
  char *q = new char[n+1];
  strncpy(q,seq+fragment_start,n);
  q[n] = 0;
  return q;
}



//------------------------------------------------------//
//   Line-based (horizontal) operations                 //
//------------------------------------------------------//


//---------RunAlign-----------
//
void GenomicRegionSAM::RunAlign(Chromosomes *C)
{
  // extract reference sequence
  char *ref_seq;
  if (I.front()->STRAND=='+') ref_seq = GetSeq(C); 
  else {
    GenomicRegion::ModifyStrand((char*)"+"); 							// if this is not done, GetSeq will rev-comp the sequence
    ref_seq = GetSeq(C);
    GenomicRegion::ModifyStrand((char*)"-"); 							// restore strands
  }
  if (ref_seq==NULL) return;
  printf(">");
  Print();

  // print (gapped) fragment sequence 
  char *gapped_fragment = GetGappedSequence(SEQ,(char*)"DP");
  printf("%s\n", gapped_fragment);

  // print (gapped) reference sequence
  char *gapped_reference = GetGappedSequence(ref_seq,(char*)"IPS");
  printf("%s\n", gapped_reference);

  // print extended CIGAR string
  for (char *t=gapped_fragment,*r=gapped_reference; t[0]!=0; t++,r++) 
    if (t[0]=='.') printf("N");
    else if ((t[0]==ALIGNMENT_GAP)&&(r[0]==ALIGNMENT_GAP)) printf("P");
    else if (t[0]==ALIGNMENT_GAP) printf("D");
    else if ((r[0]==ALIGNMENT_GAP)&&(t[0]>='a')&&(t[0]<='z')) printf("S");
    else if (r[0]==ALIGNMENT_GAP) printf("I");
    else if (t[0]==r[0]) printf("=");
    else if (t[0]!=r[0]) printf("X");
    else printf("?");
  printf("\n");

  // print (gapped) quality score
  char *gapped_quality = GetGappedSequence(QUAL,(char*)"DP",false);
  printf("%s\n", gapped_quality);
    
  // clean-up
  delete ref_seq;
  delete gapped_fragment;
  delete gapped_reference;
  delete gapped_quality;
}


  
//---------ApplyBounds-----------
//
bool GenomicRegionSAM::ApplyBounds(StringLIntMap *bounds)
{
  bool invalid = GenomicRegion::ApplyBounds(bounds);
  if (invalid==true) {
    InvalidateSeqData();
    UpdateCigarFromIntervals();
  }  
  return invalid;
}



//---------Center-----------
//
void GenomicRegionSAM::Center()
{
  GenomicRegion::Center();
  InvalidateSeqData();
  UpdateCigarFromIntervals();  
}



//---------Connect-----------
//
void GenomicRegionSAM::Connect()
{
  if (I.size()>1) {
    for (char *p=CIGAR; p[0]!=0; p++) if (p[0]=='N') p[0] = 'D';				// convert skipped regions into deletions
    for (CIGARTokens::iterator it=T.begin(); it!=T.end(); it++) if (it->second=='N') it->second = 'D';
    GenomicRegion::Connect();
  }
}



//---------Diff-----------
//
void GenomicRegionSAM::Diff()
{
  GenomicRegion::Diff();		
  InvalidateSeqData();
  UpdateCigarFromIntervals();  
}




//---------RunCalcDistances-----------
//
void GenomicRegionSAM::RunCalcDistances(char *op1, char *op2)
{
  GenomicRegion::RunCalcDistances(op1,op2);
}



//---------Divide-----------
//
void GenomicRegionSAM::Divide()
{
  PrintError("the result of the 'divide' operation cannot be represented in SAM format; convert to REG first and try again!");
}



//---------RunDivide-----------
//
void GenomicRegionSAM::RunDivide()
{
  Divide();
  Print();
}



//---------Fix-----------
//
void GenomicRegionSAM::Fix()
{
  GenomicRegion::Fix();
}



//---------Intersect-----------
//
void GenomicRegionSAM::Intersect()
{
  PrintError("the result of the 'intersect' operation is empty or trivial for SAM inputs!");
}




//---------Randomize-----------
//
void GenomicRegionSAM::Randomize(gsl_rng *random_generator, StringLIntMap *bounds)
{
  GenomicRegion::Randomize(random_generator,bounds);
  InvalidateSeqData();
}



//---------ModifyPos-----------
//
void GenomicRegionSAM::ModifyPos(const char *position_op, long int position_shift)
{
  GenomicRegion::ModifyPos(position_op,position_shift);
  GenomicRegion::Union();
  InvalidateSeqData();
  UpdateCigarFromIntervals();  
}




//---------Select-----------
//
void GenomicRegionSAM::Select(bool first, bool last, bool from5p, bool from3p)
{
  GenomicRegion::Select(first,last,from5p,from3p);
  InvalidateSeqData();
  UpdateCigarFromIntervals();  
}




//---------ShiftPos-----------
//
void GenomicRegionSAM::ShiftPos(long int start_shift, long int stop_shift, bool strand_aware)
{
  GenomicRegion::ShiftPos(start_shift,stop_shift,strand_aware);
  GenomicRegion::Union();
  InvalidateSeqData();
  UpdateCigarFromIntervals();  
}



//---------RunShuffle-----------
//
void GenomicRegionSAM::RunShuffle(gsl_rng *random_generator, GenomicRegionSet *refReg, StringLIntMap *index, StringVecLIntMap *loc)
{
  InvalidateSeqData();
  UpdateCigarFromIntervals();  
  GenomicRegion::RunShuffle(random_generator,refReg,index,loc);
}



//---------Sort-----------
//
void GenomicRegionSAM::Sort()
{
  return; 		// no need to do anything because SAM regions are by definition sorted
}


  
//---------RunSplit-----------
//
void GenomicRegionSAM::RunSplit()
{
  CIGARTokens::iterator tok_it = T.begin();
  char *seq = SEQ;
  char *qual = QUAL;
  for (GenomicIntervalSet::iterator i=I.begin(); i!=I.end(); i++) {
    printf("%s\t%ld\t%s\t%ld\t%ld\t", LABEL, FLAG, (*i)->CHROMOSOME, (*i)->START, MAPQ);
    long int fragment_len = 0;
    for ( ; (tok_it!=T.end())&&(tok_it->second!='N'); tok_it++) {
      printf("%ld%c", tok_it->first, tok_it->second);
      if (strchr(fragCIGARop,tok_it->second)!=NULL) fragment_len += tok_it->first;
    }
    if (tok_it!=T.end()) tok_it++;
    printf("\t%s\t%ld\t%ld\t", RNEXT, PNEXT, TLEN);									// NOTE: TLEN is not modified
	if (strcmp(seq,"*")==0) printf("*");
	else for (long int k=0; k<fragment_len; k++) { printf("%c", seq[0]); seq++; }
    printf("\t");
	if (strcmp(qual,"*")==0) printf("*");
    else for (long int k=0; k<fragment_len; k++) { printf("%c", qual[0]); qual++; }
    if (n_tokens>11) printf("\t%s", OPTIONAL);
    printf("\n");
  }
}



//---------ModifyStrand-----------
//
void GenomicRegionSAM::ModifyStrand(char *strand_op)
{
  char old_strand = I.front()->STRAND;
  GenomicRegion::ModifyStrand(strand_op);
  UpdateFlagFromStrand(I.front()->STRAND);
  if (old_strand!=I.front()->STRAND) ReverseGene(SEQ,strlen(SEQ));
}


  
//---------Union-----------
//
void GenomicRegionSAM::Union()
{
  return; 		// no need to do anything because SAM regions are by definition non-overlapping
}



//---------PrintWindows-----------
//
void GenomicRegionSAM::PrintWindows(long int win_step, long int win_size)
{
  PrintError("the 'windows' operation has not yet been implemented for SAM inputs; convert to REG/BED first and try again!");
}



//---------------------------------------------------------------------------------------------//
// END CLASS: GenomicRegionSAM                                                                 //
//---------------------------------------------------------------------------------------------//
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//











  
  

//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//
//---------------------------------------------------------------------------------------------//
// CLASS: GenomicRegionGFFToREG                                                                //
//---------------------------------------------------------------------------------------------//


//---------Constructor-----------
//
GenomicRegionGFFToREG::GenomicRegionGFFToREG(FileBuffer *B)
{
  char *next = B->Get();
  n_line = B->n_line;
  Read(next,n_line); 
  B->Next();
}



//---------Constructor-----------
//
GenomicRegionGFFToREG::GenomicRegionGFFToREG(char *inp, long int n_line)
{
  this->n_line = n_line;
  Read(inp,n_line); 
}



//---------Destructor-----------
//
GenomicRegionGFFToREG::~GenomicRegionGFFToREG()
{

}


//---------Read-----------
//
void GenomicRegionGFFToREG::Read(char *inp, long int n_line)
{
  char sep = '\t';
  long int n_tokens = CountTokens(inp,sep);
  if ((n_tokens<8)||(n_tokens>10)) PrintError("wrong number of tokens for GFF format!");
  char *SEQNAME = GetNextToken(&inp,sep);
  /*SOURCE =*/ StrCopy(GetNextToken(&inp,sep));
  /*FEATURE =*/ StrCopy(GetNextToken(&inp,sep));
  long int START = atol(GetNextToken(&inp,sep));
  long int END = atol(GetNextToken(&inp,sep));
  /*SCORE =*/ StrCopy(GetNextToken(&inp,sep));
  char STRAND = GetNextToken(&inp,sep)[0];
  /*FRAME =*/ GetNextToken(&inp,sep)[0];
  /*LABEL =*/ n_tokens==8?StrCopy("_"):StrCopy(GetNextToken(&inp,sep));
  /*COMMENT =*/ n_tokens>9?StrCopy(GetNextToken(&inp,sep)):NULL;
  I.push_back(new GenomicInterval(SEQNAME,STRAND,START,END,n_line));
}
  
  
//---------------------------------------------------------------------------------------------//
// END CLASS: GenomicRegionGFFToREG                                                            //
//---------------------------------------------------------------------------------------------//
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//





  
  
  
  
  
   
  
  
  

  
  
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//
//---------------------------------------------------------------------------------------------//
// CLASS: GenomicRegionGFF                                                                     //
//---------------------------------------------------------------------------------------------//


//---------Constructor-----------
//
GenomicRegionGFF::GenomicRegionGFF(FileBuffer *B)
{
  char *next = B->Get();
  n_line = B->n_line;
  Read(next,n_line); 
  B->Next();
}



//---------Constructor-----------
//
GenomicRegionGFF::GenomicRegionGFF(char *inp, long int n_line)
{
  this->n_line = n_line;
  Read(inp,n_line); 
}



//---------Destructor-----------
//
GenomicRegionGFF::~GenomicRegionGFF()
{
  delete SOURCE;
  delete FEATURE;
  delete SCORE;
  if (COMMENT!=NULL) delete COMMENT;
}


//---------Read-----------
//
void GenomicRegionGFF::Read(char *inp, long int n_line)
{
  char sep = '\t';
  n_tokens = CountTokens(inp,sep);
  if ((n_tokens<8)||(n_tokens>10)) PrintError("wrong number of tokens for GFF format!");
  char *SEQNAME = GetNextToken(&inp,sep);
  SOURCE = StrCopy(GetNextToken(&inp,sep));
  FEATURE = StrCopy(GetNextToken(&inp,sep));
  long int START = atol(GetNextToken(&inp,sep));
  long int END = atol(GetNextToken(&inp,sep));
  SCORE = StrCopy(GetNextToken(&inp,sep));
  char STRAND = GetNextToken(&inp,sep)[0];
  FRAME = GetNextToken(&inp,sep)[0];
  LABEL = n_tokens==8?StrCopy("_"):StrCopy(GetNextToken(&inp,sep));
  COMMENT = n_tokens>9?StrCopy(GetNextToken(&inp,sep)):NULL;
  I.push_back(new GenomicInterval(SEQNAME,STRAND,START,END,n_line));
}
  
  

//---------Print-----------
//
void GenomicRegionGFF::Print(FILE *file_ptr)
{
  if (I.size()==0) return;
  fprintf(file_ptr, "%s\t%s\t%s\t%ld\t%ld\t%s\t%c\t%c", I.front()->CHROMOSOME, SOURCE, FEATURE, I.front()->START, I.front()->STOP, SCORE, I.front()->STRAND, FRAME);
  if (n_tokens>8) fprintf(file_ptr, "\t%s", LABEL);
  if (n_tokens>9) fprintf(file_ptr, "\t%s", COMMENT);
  fprintf(file_ptr, "\n");
}



//---------PrintConstrained-----------
//
void GenomicRegionGFF::PrintConstrained(GenomicRegion *r, bool merge_labels)
{
  long int new_start = max(I.front()->START,r->I.front()->START);
  long int new_stop = min(I.front()->STOP,r->I.back()->STOP);
  if (new_start>new_stop) return;
  printf("%s\t%s\t%s\t%ld\t%ld\t%s\t%c\t%c", I.front()->CHROMOSOME, SOURCE, FEATURE, new_start, new_stop, SCORE, I.front()->STRAND, FRAME);
  if (n_tokens>8) {
    if (merge_labels) printf("\t%s:%s", LABEL, r->LABEL);
    else printf("\t%s", LABEL);
  }
  if (n_tokens>9) printf("\t%s", COMMENT);
  printf("\n");
}




//---------PrintModified-----------
//
void GenomicRegionGFF::PrintModified(char *label, long int start, long int stop)
{
  if (IsCompatible(false)==false) PrintError("intervals should have the same chromosome/strand for this operation!");
  printf("%s\t%s\t%s\t%ld\t%ld\t%s\t%c\t%c", I.front()->CHROMOSOME, SOURCE, FEATURE, start, stop, SCORE, I.front()->STRAND, FRAME);
  if (n_tokens>8) printf("\t%s", label);
  if (n_tokens>9) printf("\t%s", COMMENT);
  printf("\n");
}



//---------Diff-----------
//
void GenomicRegionGFF::Diff()
{
  PrintError("the result of the 'difference' operation on GFF inputs is always empty!");
}



//---------RunCalcDistances-----------
//
void GenomicRegionGFF::RunCalcDistances(char *op1, char *op2)
{
  PrintError("the result of the 'distances' operation on GFF inputs is always empty!");
}



//---------Divide-----------
//
void GenomicRegionGFF::Divide()
{
  PrintError("the result of the 'divide' operation cannot be represented in GFF format; convert to REG first and try again!");
}



//---------RunDivide-----------
//
void GenomicRegionGFF::RunDivide()
{
  Divide();
  Print();
}



//---------Intersect-----------
//
void GenomicRegionGFF::Intersect()
{
  PrintError("the result of the 'intersect' operation is empty or trivial for GFF inputs!");
}



//---------Select-----------
//
void GenomicRegionGFF::Select(bool first, bool last, bool from5p, bool from3p)
{
  PrintError("the result of the 'select' operation is empty or trivial for GFF inputs!");
}



//---------RunSplit-----------
//
void GenomicRegionGFF::RunSplit()
{
  Print();
}



//---------PrintWindows-----------
//
void GenomicRegionGFF::PrintWindows(long int win_step, long int win_size)
{
  PrintError("the 'windows' operation has not yet been implemented for GFF inputs; convert to REG/BED first and try again!");
}




//---------------------------------------------------------------------------------------------//
// END CLASS: GenomicRegionGFF                                                                 //
//---------------------------------------------------------------------------------------------//
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//





  
  




  
  

//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//
//---------------------------------------------------------------------------------------------//
// CLASS: GenomicRegionSet                                                                     //
//---------------------------------------------------------------------------------------------//



//---------Constructor--------
//
GenomicRegionSet::GenomicRegionSet(char *file, unsigned long int buffer_size, bool verbose, bool load_in_memory, bool hide_header)
{
  // initialize
  this->file = file==NULL?NULL:StrCopy(file);
  this->file_ptr = NULL;
  this->buffer_size = buffer_size;
  this->verbose = verbose;
  this->from_stdin = file==NULL;
  this->load_in_memory = from_stdin==true?false:load_in_memory;
  this->hide_header = hide_header;
  this->buffer = NULL;
  Init();
}



//---------Constructor--------
//
GenomicRegionSet::GenomicRegionSet(FILE *file_ptr, unsigned long int buffer_size, bool verbose, bool load_in_memory, bool hide_header)
{
  // initialize
  this->file = NULL;
  this->file_ptr = file_ptr;
  this->buffer_size = buffer_size;
  this->verbose = verbose;
  this->from_stdin = file==NULL;
  this->load_in_memory = from_stdin==true?false:load_in_memory;
  this->buffer = NULL;
  Init();
}



//---------Destructor--------
//
GenomicRegionSet::~GenomicRegionSet()
{
  if (file!=NULL) delete file;
  if (buffer!=NULL) delete buffer;
  if (R!=NULL) delete [] R;
}



//---------ProcessFileHeader--------
//
void GenomicRegionSet::ProcessFileHeader(bool hide)
{
  char *next = buffer->Next(); 
  if (next==NULL) format = "EMPTY";
  else if ((strstr(next,"browser ")==next)||(strstr(next,"track ")==next)) {	// UCSC browser track info
    format = "";
    for ( ; (next!=NULL)&&((strstr(next,"browser ")==next)||(strstr(next,"track ")==next)); next=buffer->Next()) if (hide==false) printf("%s\n", next); 
  }
  else if (next[0]=='@') { 											// SAM header
    format = "SAM"; 
    for ( ; (next!=NULL)&&(next[0]=='@'); next=buffer->Next()) if (hide==false) printf("%s\n", next); 
  }
  else if ((next[0]=='#')&&(next[1]=='#')) { 						// GFF header
    format = "GFF"; 
    for ( ; (next!=NULL)&&(next[0]=='#')&&(next[1]=='#'); next=buffer->Next()) if (hide==false) printf("%s\n", next); 
  }
  else format = "";
}



//---------DetectFileFormat--------
//
void GenomicRegionSet::DetectFileFormat()
{
  ProcessFileHeader(hide_header);
  char *input_line = buffer->Get();
  if (input_line==NULL) format = "EMPTY";
  if (format!="") return;
  if (input_line[0]=='>') format = "SEQ";
  else {
    long int n_tokens = CountTokens(input_line,'\t');
    if (n_tokens==1) format = "BED";
    if (n_tokens==2) format = "REG";
    if ((n_tokens>=3)&&(n_tokens<=6)) format = "BED";
    else if (n_tokens>=6) {
      char *pp = StrCopy(input_line);
      char *p = pp;
      char *t;
      for (int i=1; i<=6; i++) t = GetNextToken(&p,'\t');
      if ((strchr(t,'+')!=NULL)||(strchr(t,'-')!=NULL)) format = "BED";
      else if (n_tokens>=11) format = "SAM";
      else if ((n_tokens>=8)&&(n_tokens<=10)) format = "GFF";
      delete pp;
    }
  }
}



//---------CreateGenomicRegion--------
//
GenomicRegion *GenomicRegionSet::CreateGenomicRegion(FileBuffer *buffer)
{
  if (format=="REG") return new GenomicRegion(buffer);
  else if (format=="SEQ") return new GenomicRegionSEQ(buffer);
  else if (format=="BED") return new GenomicRegionBED(buffer);
  else if (format=="SAM") return new GenomicRegionSAM(buffer);
  else if (format=="GFF") return new GenomicRegionGFF(buffer);
  PrintError("unsupported input format!\n");
  return NULL;
}



//---------Init--------
//
void GenomicRegionSet::Init()
{
  // open region set file
  buffer = file_ptr==NULL?CreateFileBuffer(file,buffer_size):new FileBufferText(file_ptr,buffer_size);
  if (load_in_memory) {
    ProcessFileHeader(true);
    n_regions = 0;
    for (char *next=buffer->Get(); next!=NULL; next=buffer->Next()) n_regions++;
    delete buffer; 
    buffer = file_ptr==NULL?CreateFileBuffer(file,buffer_size):new FileBufferText(file_ptr,buffer_size);
    DetectFileFormat();   
    if (format=="SEQ") n_regions /= 2;
  }
  else {
    n_regions = 1;
    DetectFileFormat(); 
  }
  if (format=="EMPTY") n_regions = 0;
  
  // print info
  if (verbose) { 
    cerr << "Reading from '" << (file==NULL?"<standard input>":file) << "'; ";
    cerr << "read-from-stdin = " << (from_stdin?"true":"false") << "; "; 
    cerr << "load-in-memory = " << (load_in_memory?"true":"false") << "; "; 
    cerr << "number of regions = " << n_regions << "; "; 
    cerr << "format = " << format << "\n";
  }  

  // read data
  R = n_regions>0?new GenomicRegion*[n_regions]:NULL;
  if (verbose) progress.Init("Reading regions...",n_regions);
  for (long int n=0; n<n_regions; n++) {
    R[n] = CreateGenomicRegion(buffer);
    progress.Check();
  }
  progress.Done();
  r_index = 0;
}



//---------Reset--------
//
void GenomicRegionSet::Reset()
{
  if (from_stdin) PrintError("stdin cannot be reset!\n"); 

  if (load_in_memory) {

  }
  else {
    if (buffer!=NULL) delete buffer;
    buffer = CreateFileBuffer(file,buffer_size);
    n_regions = 1;
    buffer->Next();
    R[0] = CreateGenomicRegion(buffer);
  }
  
  r_index = 0;
}



//---------Begin--------
//
GenomicRegion *GenomicRegionSet::Begin()
{
  Reset();
  return Get();
}



//---------Get--------
//
GenomicRegion *GenomicRegionSet::Get()
{
  if (n_regions==0) return NULL;
  return load_in_memory==false ? R[0] : (r_index>=n_regions ? NULL : R[r_index]);
}



//---------Next--------
//
GenomicRegion *GenomicRegionSet::Next(bool retain_current)
{
  if (load_in_memory==false) {
    if (buffer->Get()==NULL) return NULL;
    if ((R[0]!=NULL)&&(retain_current==false)) delete R[0];
    R[0] = CreateGenomicRegion(buffer);
    return R[0];
  }
  else {
    ++r_index;
    return r_index>=n_regions?NULL:R[r_index];
  }
}




//---------Next--------
//
GenomicRegion *GenomicRegionSet::Next(bool sorted_by_strand, bool retain_current)
{
  GenomicRegion *r0 = Get();
  if (r0==NULL) return NULL;
  GenomicRegion *r = Next(true);
  if ((r!=NULL)&&(r->I.front()->IsBefore(r0->I.front(),sorted_by_strand))) r->PrintError("input regions are not sorted (sorted-by-strand = " + (string)(sorted_by_strand?"true":"false") + ")!");
  if ((load_in_memory==false)&&(retain_current==false)) delete r0;
  return r;
}




//---------Release--------
//
void GenomicRegionSet::Release(GenomicRegion *r)
{
  if ((load_in_memory==false)&&(r!=NULL)) delete r;
}



//---------GetRetain--------
//
GenomicRegion *GenomicRegionSet::GetRetain()
{
  GenomicRegion *r = Get();
  if (load_in_memory==false) R[0] = NULL; 
  return r;
}



//---------PrintError-----------
//
void GenomicRegionSet::PrintError(string error_msg)
{
  fprintf(stderr,"\n");
  fprintf(stderr, "Error: %s\n", error_msg.c_str());
  exit(1);
}



//---------CountRegions--------
//
long int GenomicRegionSet::CountRegions(bool use_label_counts)
{
  if (from_stdin) PrintError("this operation cannot be executed if the file is being read from stdin!");  

  Progress PRG("Counting...",n_regions);
  long int n_counts = 0;
  for (GenomicRegion *r=Begin(); r!=NULL; r=Next(),PRG.Check()) n_counts += use_label_counts?max(1L,atol(r->LABEL)):1;
  PRG.Done();
  //Reset();

  return n_counts;
}




//---------------------------------------//
//  Line-based (horizontal) operations   //
//---------------------------------------//


//---------RunAlign--------
//
void GenomicRegionSet::RunAlign(Chromosomes *C)
{
  if (n_regions==0) return;
  if (C==NULL) PrintError("chromosome sequences required for this operation!");
  Progress PRG("Aligning input sequences...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) r->RunAlign(C);
  PRG.Done();
}



//---------RunConvertToBED--------
//
void GenomicRegionSet::RunConvertToBED(char *title, char *color, char *position, bool convert_chromosome)
{
  if (n_regions==0) return;
  Progress PRG("Printing in BED format...",n_regions);
  if (strlen(position)>0) printf("browser position %s\n", position);
  if (strlen(title)>0) printf("track name='%s' itemRgb=On visibility=1\n", title);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) r->PrintBEDFormat(color,convert_chromosome);
  PRG.Done();
}



//---------RunBounds--------
//
void GenomicRegionSet::RunBounds(StringLIntMap *bounds)
{
  if (n_regions==0) return;
  Progress PRG("Printing bounded regions...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check())
    if (bounds&&(bounds->find(r->I.front()->CHROMOSOME)!=bounds->end())) { r->ApplyBounds(bounds); r->Print(); }
  PRG.Done();
}



//---------RunCenter--------
//
void GenomicRegionSet::RunCenter()
{
  if (n_regions==0) return;
  Progress PRG("Printing center of intervals...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) { r->Center(); r->Print(); }
  PRG.Done();
}



//---------RunModifyChromosomeNames--------
//
void GenomicRegionSet::RunModifyChromosomeNames(char *chromosome_conversion_file)
{
  if (n_regions==0) return;
  Progress PRG("Modifying chromosome names...",n_regions);
  if (strlen(chromosome_conversion_file)==0) PrintError("chromosome name conversion file is missing!\n"); 
  FileBufferText buffer(chromosome_conversion_file);
  map<string,string> chromosome_conversion_map;
  for (char *inp=buffer.Next(); inp!=NULL; inp=buffer.Next()) {
    char *x = GetNextToken(&inp," \t"); 
    char *y = GetNextToken(&inp," \t"); 
	chromosome_conversion_map[x] = y;
  }
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) { r->ModifyChromosomeNames(chromosome_conversion_map); r->Print(); }
  PRG.Done();
}



//---------RunConnect--------
//
void GenomicRegionSet::RunConnect()
{
  if (n_regions==0) return;
  Progress PRG("Printing connected regions...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) { r->Connect(); r->Print(); }
  PRG.Done();
}



//---------RunDiff--------
//
void GenomicRegionSet::RunDiff()
{
  if (n_regions==0) return;
  Progress PRG("Printing difference...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) { r->Diff(); r->Print(); }
  PRG.Done();
}



//---------RunCalcDistances--------
//
void GenomicRegionSet::RunCalcDistances(char *op1, char *op2)
{
  if (n_regions==0) return;
  Progress PRG("Printing distance...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) r->RunCalcDistances(op1,op2);
  PRG.Done();
}



//---------RunDivide--------
//
void GenomicRegionSet::RunDivide()
{
  if (n_regions==0) return;
  Progress PRG("Printing divided regions...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) r->RunDivide(); 
  PRG.Done();
}



//---------RunFix--------
//
void GenomicRegionSet::RunFix()
{
  if (n_regions==0) return;
  Progress PRG("Printing fixed regions...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) { r->Fix(); r->Print(); }
  PRG.Done();
}



//---------RunIntersection--------
//
void GenomicRegionSet::RunIntersection()
{
  if (n_regions==0) return;
  Progress PRG("Printing intersections...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) { r->Intersect(); r->Print(); }
  PRG.Done();
}



//---------RunSize--------
//
void GenomicRegionSet::RunSize()
{
  if (n_regions==0) return;
  Progress PRG("Printing region size...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) r->RunSize();

  PRG.Done();
}



//---------RunModifyPos--------
//
void GenomicRegionSet::RunModifyPos(const char *position_op, long int position_shift)
{
  if ((strcmp(position_op,"1")!=0)&&(strcmp(position_op,"2")!=0)&&(strcmp(position_op,"5p")!=0)&&(strcmp(position_op,"3p")!=0)&&(strcmp(position_op,"c")!=0)) { cerr << "Error: unknown operation!\n"; exit(1); }
  if (n_regions==0) return;
  Progress PRG("Printing modified regions...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) { r->ModifyPos(position_op,position_shift); r->Print(); }
  PRG.Done();
}



//---------RunConvertToREG--------
//
void GenomicRegionSet::RunConvertToREG(bool compact)
{
  if (n_regions==0) return;
  Progress PRG("Printing regions...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) r->PrintREG(compact);
  PRG.Done();
}



//---------RunRandomize--------
//
void GenomicRegionSet::RunRandomize(gsl_rng *random_generator, StringLIntMap *bounds)
{
  if (n_regions==0) return;
  Progress PRG("Printing randomized regions...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) { r->Randomize(random_generator,bounds); r->Print(); }
  PRG.Done();
}



//---------RunSelect--------
//
void GenomicRegionSet::RunSelect(bool first, bool last, bool from5p, bool from3p)
{
  if (n_regions==0) return;
  Progress PRG("Printing selected regions...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) { r->Select(first,last,from5p,from3p); r->Print(); }
  PRG.Done();
}



//---------RunShiftPos--------
//
void GenomicRegionSet::RunShiftPos(long int start_shift, long int stop_shift, bool strand_aware)
{
  if (n_regions==0) return;
  Progress PRG("Printing shifted regions...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) { r->ShiftPos(start_shift,stop_shift,strand_aware); r->Print(); }
  PRG.Done();
}



//------RunShuffle------
//
void GenomicRegionSet::RunShuffle(gsl_rng *random_generator, char *ref_reg_file)
{
  bool sorted_by_strand = true;

  // step #1
  GenomicRegionSet refReg(ref_reg_file,10000,verbose,true,false);
  Progress PRG1("Processing reference regions...",refReg.n_regions);
  StringLIntMap index;
  StringVecLIntMap loc;
  for (long int rr=0; rr<refReg.n_regions; rr++) {
    if (refReg.R[rr]->I.size()!=1) refReg.R[rr]->PrintError("reference regions must be single-interval regions for this operation!");
    if ((rr>0)&&(refReg.R[rr]->IsBefore(refReg.R[rr-1],sorted_by_strand))) refReg.R[rr]->PrintError("reference regions must be sorted for this operation!"); 
    if ((rr>0)&&(refReg.R[rr]->OverlapsWith(refReg.R[rr-1],!sorted_by_strand))) refReg.R[rr]->PrintError("reference regions must be non-overlapping for this operation!"); 
    string v = (string)refReg.R[rr]->I.front()->CHROMOSOME + (refReg.R[rr]->I.front()->STRAND=='+'?"+":"-");
    long int l = refReg.R[rr]->I.front()->STOP - refReg.R[rr]->I.front()->START + 1;
    if (index.find(v)==index.end()) { index[v] = rr; loc[v].push_back(0); loc[v].push_back(l); }
    else loc[v].push_back(loc[v].back()+l);
    PRG1.Check();
  }
  PRG1.Done();
  //for (map2::iterator x=loc.begin(); x!=loc.end(); x++) { string v = x->first; cerr << v << " (" << index[v] << ") = [";  for (vector<long int>::iterator y=loc[v].begin(); y!=loc[v].end(); y++) cerr << ' ' << *y; cerr << " ]\n"; }

  // step #2
  Progress PRG2("Shuffling regions...",1);
  for (GenomicRegion *r=Get(); r!=NULL; PRG2.Check(),r=Next()) r->RunShuffle(random_generator,&refReg,&index,&loc);
  PRG2.Done();
}


//---------RunSort--------
//
void GenomicRegionSet::RunSort()
{
  if (n_regions==0) return;
  Progress PRG("Printing sorted regions...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) { r->Sort(); r->Print(); }
  PRG.Done();
}


//---------RunSplit--------
//
void GenomicRegionSet::RunSplit()
{
  if (n_regions==0) return;
  Progress PRG("Splitting regions...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) r->RunSplit();
  PRG.Done();
}


/*
//---------PrintModifyLabel--------
//
void GenomicRegionSet::PrintModifyLabel()
{
  if (n_regions==0) return;
  Progress PRG("Printing regions with modified labels...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) r->PrintModifyLabel();
  PRG.Done();
}
*/


//---------PrintSeqLength--------
//
void GenomicRegionSet::PrintSeqLength(Chromosomes *C)
{
  if (n_regions==0) return;
  Progress PRG("Printing sequence length...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) r->PrintSeqLength(C);
  PRG.Done();
}


/*
//---------PrintNoGaps--------
//
void GenomicRegionSet::PrintNoGaps()
{
  if (n_regions==0) return;
  Progress PRG("Printing no-gap regions...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) r->PrintNoGaps();
  PRG.Done();
}
*/


//---------PrintRemoveN--------
//
void GenomicRegionSet::PrintRemoveN(Chromosomes *C)
{
  if (n_regions==0) return;
  Progress PRG("Printing connected regions...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) r->PrintRemoveN(C);
  PRG.Done();
}



//---------PrintOffsetFormat--------
//
void GenomicRegionSet::PrintOffsetFormat(char *op, bool fraction)
{
  if (n_regions==0) return;
  Progress PRG("Printing in offset format...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) r->PrintOffsetFormat(op,fraction);
  PRG.Done();
}



//---------PrintReversePos--------
//
void GenomicRegionSet::PrintReversePos(StringLIntMap *bounds)
{
  if (n_regions==0) return;
  Progress PRG("Printing connected regions...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) { r->ReversePos(bounds); r->Print(); }
  PRG.Done();
}



//---------PrintSearch--------
//
void GenomicRegionSet::PrintSearch(char *pattern, bool header, bool summary)
{
  if (n_regions==0) return;
  Progress PRG("Searching...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) r->PrintSearch(pattern,header,summary);
  PRG.Done();
}



//---------StoreInTempFile-----------
//
GenomicRegionSet *GenomicRegionSet::StoreInTempFile()
{
  bool sorted_by_strand = true;
  GenomicRegion *r = Get();
  if ((r==NULL)||(r->I.front()->STRAND!='+')) return NULL;
  GenomicRegion *r0 = GetRetain();
  FILE *temp_file_ptr = tmpfile();
  r0->Print(temp_file_ptr);
  //long int N = 1;
  //fprintf(stderr, "chr=%s: ", r0->I.front()->CHROMOSOME); 
  r = Next();
  while ((r!=NULL)&&(r->IsCompatibleWith(r0,!sorted_by_strand)==true)) {
    if (r->IsBefore(r0,sorted_by_strand)) r->PrintError("regions must be sorted for this operation!"); 
    delete r0;
    r0 = GetRetain();
    r0->Print(temp_file_ptr);
    r = Next();
    progress.Check();
    //N++;
  }
  //fprintf(stderr, "%ld regions stored\n", N);
  rewind(temp_file_ptr);
  return new GenomicRegionSet(temp_file_ptr,10000,false,false,false); 
}



//---------PrintMergeSort-----------
//
GenomicRegion *GenomicRegionSet::PrintMergeSort(GenomicRegionSet *rtempSet, char strand)
{
  GenomicRegion *rtemp = rtempSet==NULL?NULL:rtempSet->Get();
  GenomicRegion *r = Get();
  if ((rtemp==NULL)&&(r==NULL)) return NULL;
  char *chromosome = StrCopy(rtemp==NULL?r->GetChromosome():rtemp->GetChromosome());
  //fprintf(stderr, "chr=%s merging...\n", chromosome);
  while (rtemp&&r&&(strcmp(r->GetChromosome(),chromosome)==0)) {
    if (rtemp->IsPosBefore(r)) { rtemp->SetStrand(strand); rtemp->Print(); rtemp = rtempSet->Next(); }
    else { r->SetStrand(strand); r->Print(); r = Next(); progress.Check(); }
  }
  //fprintf(stderr, "done merging, printing leftovers...\n");
  while (rtemp) { rtemp->SetStrand(strand); rtemp->Print(); rtemp = rtempSet->Next(); }
  while (r&&(strcmp(r->GetChromosome(),chromosome)==0)) { r->SetStrand(strand); r->Print(); r = Next(); progress.Check(); }
  delete chromosome;
  return r;
}



//---------PrintReverse-----------
//
GenomicRegion *GenomicRegionSet::PrintReverse(GenomicRegionSet *rtempSet)
{
  GenomicRegion *rtemp = rtempSet==NULL?NULL:rtempSet->Get();
  GenomicRegion *r = Get();
  if ((rtemp==NULL)&&(r==NULL)) return NULL;
  char *chromosome = StrCopy(rtemp==NULL?r->GetChromosome():rtemp->GetChromosome());
  while (r&&(strcmp(r->GetChromosome(),chromosome)==0)) { r->ReverseStrand(); r->Print(); r = Next(); progress.Check(); }
  while (rtemp) { rtemp->ReverseStrand(); rtemp->Print(); rtemp = rtempSet->Next(); }
  delete chromosome;
  return r;
}



//---------RunModifyStrand--------
//
void GenomicRegionSet::RunModifyStrand(char *strand_op, bool sorted)
{
  if (sorted) RunModifyStrandSorted(strand_op);
  else {
    if (n_regions==0) return;
    Progress PRG("Printing modified regions...",n_regions);
    for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) { r->ModifyStrand(strand_op); r->Print(); }
    PRG.Done();
  }
}



//---------RunModifyStrandSorted-----------
//
void GenomicRegionSet::RunModifyStrandSorted(char *strand_op)
{
  if (format=="SEQ") PrintError("SEQ format for this operation is not implemented yet!");        // NOTE: this can be fixed!!
  if (n_regions==0) return;
  if (strcmp(strand_op,"r")==0) {
    progress.Init("Printing modified regions...",verbose?0:1);
    while (true) { 
      GenomicRegionSet *rtempSet = StoreInTempFile();
      GenomicRegion *r = PrintReverse(rtempSet);
      delete rtempSet;
      if (r==NULL) break;
    }
    progress.Done();
  }
  else if ((strcmp(strand_op,"+")==0)||(strcmp(strand_op,"-")==0)) {
    progress.Init("Printing modified regions...",verbose?0:1);
    while (true) { 
      GenomicRegionSet *rtempSet = StoreInTempFile();
      GenomicRegion *r = PrintMergeSort(rtempSet,strand_op[0]);
      if (rtempSet!=NULL) { fclose(rtempSet->file_ptr); delete rtempSet; }
      if (r==NULL) break;
    }
    progress.Done();
  }
  else PrintError("unknown strand operation!"); 
}




//---------RunUnion--------
//
void GenomicRegionSet::RunUnion()
{
  if (n_regions==0) return;
  Progress PRG("Printing unions...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) { r->Union(); r->Print(); } 
  PRG.Done();
}



//---------PrintVerifySeq--------
//
void GenomicRegionSet::PrintVerifySeq(Chromosomes *C, bool ignore)
{
  if (n_regions==0) return;
  Progress PRG("Verifying sequences...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) r->PrintVerifySeq(C,ignore);
  PRG.Done();
}



//---------RunConvertToWIG--------
//
void GenomicRegionSet::RunConvertToWIG(char *title, char *color, char *position, char *options, long int span, bool convert_chromosome)
{
  if (n_regions==0) return;
  bool ignore_strand = true;
  Progress PRG("Printing in Wiggle format...",n_regions);
  if (strlen(position)>0) printf("browser position %s\n", position);
  if (strlen(title)>0) printf("track type=wiggle_0 name='%s' color=%s %s\n", title, color, options);
  for (GenomicRegion *r=Get(); r!=NULL; ) {
    if (r->I.size()!=1) r->PrintError("only single-interval regions are accepted for this operation!");
    char *chromosome = StrCopy(r->I.front()->CHROMOSOME);
    printf("variableStep chrom=");
    PrintChromosome(chromosome,convert_chromosome);
    printf(" span=%ld\n", span);
    while (strcmp(r->I.front()->CHROMOSOME,chromosome)==0) {
      printf("%ld %s\n", r->I.front()->START, r->LABEL);
      r = Next(!ignore_strand,false);
      if (r==NULL) break;
      PRG.Check();
    }
    delete chromosome;
  }
  PRG.Done();
}



//---------RunSlidingWindows--------
//
void GenomicRegionSet::RunSlidingWindows(long int win_step, long int win_size)
{
  if (n_regions==0) return;
  Progress PRG("Printing windows...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) r->PrintWindows(win_step,win_size);
  PRG.Done();
}



//---------RunExtractSeq--------
//
void GenomicRegionSet::RunExtractSeq(Chromosomes *C, bool replace)
{
  if (n_regions==0) return;
  if (C==NULL) PrintError("chromosome sequences required for this operation!");
  Progress PRG("Extracting sequences...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next(),PRG.Check()) r->PrintSeq(C,replace);
  PRG.Done();
}





//---------------------------------------//
//  File-based (vertical) operations     //
//---------------------------------------//


//---------RunGlobalAnnotate--------
//
void GenomicRegionSet::RunGlobalAnnotate(StringLIntMap *bounds)
{
  if (load_in_memory==false) PrintError("need to load entire region set into memory for this operation!");
  if (n_regions==0) return;

  Progress PRG("Annotating regions...",n_regions);
  for (long int n=0; n<n_regions; n++,PRG.Check()) {
    if (R[n]->IsCompatibleSortedAndNonoverlapping()==false) R[n]->PrintError("region intervals must be compatible, sorted and non-overlapping for this operation!");
    if (R[n]->IsCompatible(false)==false) R[n]->PrintError("region intervals must have the same chromosome and strand for this operation!");
    if (bounds->find(R[n]->I.front()->CHROMOSOME)==bounds->end()) continue;
    
    long int prev_stop = ((n==0)||(strcmp(R[n]->I.front()->CHROMOSOME,R[n-1]->I.front()->CHROMOSOME)!=0)) ? 0 : R[n-1]->I.front()->STOP;
    long int next_start = ((n==n_regions-1)||(strcmp(R[n]->I.front()->CHROMOSOME,R[n+1]->I.front()->CHROMOSOME)!=0)) ? (*bounds)[R[n]->I.front()->CHROMOSOME]+1 : R[n+1]->I.front()->START;

    if (prev_stop+1<=R[n]->I.front()->START-1) {
      if (R[n]->I.front()->STRAND=='+') cout << "upstream"; else cout << "downstream";
      cout << '|' << R[n]->LABEL << '\t' << R[n]->I.front()->CHROMOSOME << ' ' << R[n]->I.front()->STRAND << ' ' << prev_stop+1 << ' ' << R[n]->I.front()->START-1 << '\n';
    }

    cout << (R[n]->I.front()->STRAND=='+'?"TSS":"TES") << '|' << R[n]->LABEL << '\t';
    R[n]->I.front()->Print(); 
    GenomicIntervalSet::iterator i=R[n]->I.begin(), j=R[n]->I.begin();
    for (i++; i!=R[n]->I.end(); i++,j++) {
      if ((*i)->START-1>=(*j)->STOP+1) {
        cout << "intron" << '|' << R[n]->LABEL << '\t' << R[n]->I.front()->CHROMOSOME << ' ' << R[n]->I.front()->STRAND << ' ' << (*j)->STOP+1 << ' ' << (*i)->START-1 << '\n';
      }
      if (*i!=R[n]->I.back()) {
        cout << "exon" << '|' << R[n]->LABEL << '\t';
        (*i)->Print();
      }
    }
    cout << (R[n]->I.front()->STRAND=='+'?"TES":"TSS") << '|' << R[n]->LABEL << '\t';
    R[n]->I.back()->Print(); 
    
    if (R[n]->I.front()->STOP+1<=next_start-1) {
      if (R[n]->I.front()->STRAND=='+') cout << "downstream"; else cout << "upstream";
      cout << '|' << R[n]->LABEL << '\t' << R[n]->I.front()->CHROMOSOME << ' ' << R[n]->I.front()->STRAND << ' ' << R[n]->I.front()->STOP+1 << ' ' << next_start-1 << '\n';
    }
  }
  PRG.Done();

}



//---------RunGlobalCluster--------
//
void GenomicRegionSet::RunGlobalCluster(bool merge)
{
  if (load_in_memory==false) { fprintf(stderr, "Sorry, need to load entire region set into memory :-(\n"); exit(1); }
  if (n_regions==0) return;

  // initialize clusters
  long int *clusters = new long int[n_regions];
  clusters[0] = 0;
  R[0]->Print();
  long int n_clusters = 1;

  // start clustering
  Progress PRG("Clustering regions...",n_regions);
  for (long int n=1; n<n_regions; n++) {
    bool new_cluster = true;
    for (long int c=0; c<n_clusters; c++) if (R[clusters[c]]->OverlapsWith(R[n],false)==true) { new_cluster = false; break; }
    if (new_cluster) { clusters[n_clusters++] = n; R[n]->Print(); }
    PRG.Check();
  }
  PRG.Done();
}



//---------RunGlobalCalcDistances--------
//
void GenomicRegionSet::RunGlobalCalcDistances(char *op1, char *op2)
{
  if (n_regions==0) return;
  bool sorted_by_strand = true;
  Progress PRG("Computing distances between successive regions...",n_regions);
  for (GenomicRegion *r=Get(),*r0=NULL; r!=NULL; r=Next()) {
    if (r->I.size()!=1) r->PrintError("this operation requires single-interval regions!");
    if (r0!=NULL) {
      if (r->I.front()->IsBefore(r0->I.front(),sorted_by_strand)) r->PrintError("input regions are not sorted (sorted-by-strand = " + (string)(sorted_by_strand?"true":"false") + ")!");
      printf("%s\t%s\t", r0->LABEL, r->LABEL);
      if (r->I.front()->IsCompatibleWith(r0->I.front(),!sorted_by_strand)) printf("%ld", r->I.front()->CalcDistanceFrom(r0->I.front(),op2,op1));
      else printf("NaN");
      printf("\n");
    }
    if (load_in_memory==false) { R[0] = NULL; if (r0!=NULL) delete r0; } 
    r0 = r;
    PRG.Check();
  }
  PRG.Done();
}


//---------RunGlobalSort--------
//
void GenomicRegionSet::RunGlobalSort(bool sorted_by_strand, int bin_bits)
{
	GenomicRegionSetBins *bins = BinGenomicRegions(this,sorted_by_strand,bin_bits);

	Progress PRG("Printing sorted regions...",n_regions);
	for (GenomicRegionSetBins::iterator p=bins->begin(); p!=bins->end(); p++) {
		ChromosomeBins *chrom_bins = p->second;
		GenomicRegionList *b = chrom_bins->b_fwd;
		while (true) {
			GenomicRegionList *bb = b;
			for (unsigned long int w=0; w<chrom_bins->n_bins; w++,bb++) {
				if (bb->size()==0) continue;
				bb->sort(CompareBinnedGenomicRegions);
				for (GenomicRegionList::iterator z=bb->begin(); z!=bb->end(); z++) { (*z)->Print(); PRG.Check(); }
			}
			if (sorted_by_strand==false) break;
			else if (b==chrom_bins->b_fwd) b = chrom_bins->b_rev;
			else break;
		}
	}
	bins->clear();
	delete bins;
	PRG.Done();
}



//-----RunGlobalInvert----------
//
void GenomicRegionSet::RunGlobalInvert(StringLIntMap *bounds)
{
  if (n_regions==0) return;
  bool sorted_by_strand = true;
  Progress PRG("Inverting regions...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; ) {
    if (r->I.size()!=1) r->PrintError("not a single-interval region!");
    GenomicRegion *r0 = r;
    if (load_in_memory==false) R[0] = NULL;
    if (bounds->find(r0->I.front()->CHROMOSOME)==bounds->end()) { fprintf(stderr, "Line %ld: chromosome %s not found!\n", buffer->n_line, r0->I.front()->CHROMOSOME); exit(1); } 
    long int chrom_size = (*bounds)[r0->I.front()->CHROMOSOME];
    if (r0->I.front()->START>1) r0->PrintModified((char*)"_",1L,r0->I.front()->START-1); 
    while (((r=Next())!=NULL)&&(r->I.front()->IsCompatibleWith(r0->I.front(),!sorted_by_strand))) {
      if (r->I.size()!=1) r->PrintError("not a single-interval region!");
      PRG.Check();
      if (r->I.front()->IsBefore(r0->I.front(),sorted_by_strand)) r->PrintError("input regions are not sorted (sorted-by-strand = " + (string)(sorted_by_strand?"true":"false") + ")!");
      if (r->I.front()->START>r0->I.front()->STOP+1) r->PrintModified((char*)"_",r0->I.front()->STOP+1,r->I.front()->START-1);
      if (load_in_memory==false) { R[0] = NULL; delete r0; }
      r0 = r;
    }
    if (r0->I.front()->STOP+1<chrom_size) r0->PrintModified((char*)"_",r0->I.front()->STOP+1,chrom_size);
    if (load_in_memory==false) delete r0;
  }
  PRG.Done();
}




//---------RunGlobalLink--------
//
void GenomicRegionSet::RunGlobalLink(bool sorted_by_strand, long int max_difference)
{
  if (n_regions==0) return;
  Progress PRG("Linking regions...",n_regions);
  for (GenomicRegion *r=Get(); r!=NULL; ) {
    if (r->I.size()!=1) r->PrintError("not a single-interval region!");
    GenomicRegion *r0 = r; 
    long int new_stop = r0->I.front()->STOP; 
    while (true) {
      if ((r=Next(sorted_by_strand,r==r0))!=NULL) {
        if (r->I.size()!=1) r->PrintError("not a single-interval region!");
        PRG.Check();
        if (r0->I.front()->IsCompatibleWith(r->I.front(),!sorted_by_strand)&&(r->I.front()->START-new_stop<=max_difference)) { new_stop = max(r->I.front()->STOP,new_stop); continue; }
      }
      r0->PrintModified((char*)"_",r0->I.front()->START,new_stop);
      break;
    }
    Release(r0);
  }
  PRG.Done();
}



//---------RunGlobalReverseOrder--------
//
void GenomicRegionSet::RunGlobalReverseOrder()
{
  if (load_in_memory==false) { fprintf(stderr, "Sorry, need to load entire region set into memory :-(\n"); exit(1); }
  Progress PRG("Reversing the order of reverse-strand regions...",n_regions);
  long int n = 0;
  while (n<n_regions) {
    if (R[n]->I.size()!=1) { fprintf(stderr, "Line %ld: not a single-interval region!\n", n+1); exit(1); }
    char *chromosome = R[n]->I.front()->CHROMOSOME;
    while ((n<n_regions)&&(R[n]->I.front()->STRAND=='+')&&(strcmp(R[n]->I.front()->CHROMOSOME,chromosome))==0) { R[n]->Print(); n++; PRG.Check(); }
    long int nn = n;
    while ((n<n_regions)&&(R[n]->I.front()->STRAND=='-')&&(strcmp(R[n]->I.front()->CHROMOSOME,chromosome))) n++;
    for (long int k=n-1; k>=nn; k--) { R[k]->Print(); PRG.Check(); }
  }
  PRG.Done();
}




//---------RunGlobalScan--------
//
void GenomicRegionSet::RunGlobalScan(StringLIntMap *bounds, long int win_step, long int win_size)
{
  if (format=="SEQ") { cerr << "Error: operation 'scan' does not accept SEQ format!\n"; exit(1); }
  if (n_regions==0) return;
  if (bounds==NULL) { cerr << "Error: operation 'scan' requires genomic bounds!\n"; exit(1); }
  cerr << "Warning: regions are assumed to be point-singleton forward-strand and sorted!\n"; 

  if (win_size%win_step!=0) { cerr << "Error: window size must be a multiple of window step in operation 'scan'!\n"; exit(1); }
  long int n_win_combine = win_size/win_step;
  Progress PRG("Scanning input intervals...",n_regions);
  string *v = new string[n_win_combine];
  GenomicRegion *r = Get();
  for (StringLIntMap::iterator chr=bounds->begin(); chr!=bounds->end(); chr++) {
    while ((r!=NULL)&&(chr->first>r->I.front()->CHROMOSOME)) { r = Next(); PRG.Check(); }
    if ((r==NULL)||(chr->first!=r->I.front()->CHROMOSOME)) cerr << "Warning: no overlap found in chromosome " << chr->first << "!\n";
    char strand = '+';
    for (long int k=0,start=1,stop=start+win_step-1; stop<=chr->second; start+=win_step,stop+=win_step,k=(k+1)%n_win_combine) {
      for (v[k]=""; (r!=NULL)&&(strcmp(r->I.front()->CHROMOSOME,chr->first.c_str())==0)&&(r->I.front()->STRAND==strand)&&(r->I.front()->START<=stop); r=Next()) {
        if (r->I.size()!=1) r->PrintError("single-interval regions expected for this operation!");
        v[k] = v[k] + " " + r->LABEL;
        PRG.Check();
      }
      if (stop>=win_size) {
        for (long int i=0; i<n_win_combine; i++) cout << v[i]; 
        cout << '\t' << chr->first << ' ' << strand << ' ' << stop-win_size+1 << ' ' << stop << '\n';
      }
    }
  }
  PRG.Done();
  //free(v);
}



//-----RunGlobalScanCount----------
//
void GenomicRegionSet::RunGlobalScanCount(StringLIntMap *bounds, char *ref_reg_file, long int win_step, long int win_size, bool ignore_strand, char preprocess, bool use_labels_as_values, long int min_reads)
{
  SortedGenomicRegionSetScanner input_scanner(this,bounds,win_step,win_size,use_labels_as_values,ignore_strand,preprocess);
  GenomicRegionSet *RefRegSet = strlen(ref_reg_file)>0?new GenomicRegionSet(ref_reg_file,10000,false,false,true):NULL;
  Progress PRG("Scanning...",1);
  for (long int v=input_scanner.Next(RefRegSet); v!=-1; v=input_scanner.Next(RefRegSet)) {
    if (v>=min_reads) {
      printf("%ld\t", v);
      input_scanner.PrintInterval();
      printf("\n");
    }
    PRG.Check();
  }
  PRG.Done();
  if (RefRegSet!=NULL) delete RefRegSet;
}




//---------RunGlobalPartition--------
//
void GenomicRegionSet::RunGlobalPartition()
{
  bool ignore_strand = false; 
  long int r = 0;
  while (r<n_regions) {
    if (R[r]->I.size()!=1) R[r]->PrintError("not a single-interval region!");
    //printf(">>> "); R[r]->Print();
    long int start = R[r]->I.front()->START;
    long int stop = R[r]->I.front()->STOP;
    printf("%s", R[r]->LABEL);
    long int rr = r+1;
    while ((rr<n_regions)&&(R[rr]->I.front()->IsCompatibleWith(R[r]->I.front(),ignore_strand))&&(R[rr]->I.front()->START==R[r]->I.front()->START)) printf(",%s", R[rr++]->LABEL);
    if ((rr<n_regions)&&(R[rr]->I.front()->IsCompatibleWith(R[r]->I.front(),ignore_strand))) {
      if (R[rr]->I.front()->START<=stop) stop = R[rr]->I.front()->START-1;
    }
    for (long int q=r; q<rr; q++) R[q]->I.front()->START = stop+1;
    printf("\t%s %c %ld %ld\n", R[r]->I.front()->CHROMOSOME, R[r]->I.front()->STRAND, start, stop);
    while ((r<n_regions)&&(R[r]->I.front()->STOP==stop)) r++;
    //for (long int q=r; q<rr; q++) { printf("-- "); R[q]->Print(); }
  }
}



//---------RunGlobalTest--------
//
void GenomicRegionSet::RunGlobalTest(bool sorted_by_strand)
{
  if (n_regions==0) return;
  Progress PRG("Testing regions...",n_regions);
  long int n_inclusions = 0;
  long int n_overlaps = 0;
  for (GenomicRegion *r=Get(),*r0=NULL; r!=NULL; r=Next()) {
    if (r->IsCompatibleSortedAndNonoverlapping()==false) r->PrintError("input regions must be compatible, sorted and non-overlapping!");
    if (r0!=NULL) {
      if (r->I.front()->IsBefore(r0->I.front(),sorted_by_strand)) r->PrintError("input regions are not sorted (sorted-by-strand = " + (string)(sorted_by_strand?"true":"false") + ")!");
      if (r->I.front()->IsCompatibleWith(r0->I.front(),!sorted_by_strand)) {
        if (r->I.front()->START<=r0->I.front()->STOP) {
          if (r->I.front()->STOP<=r0->I.front()->STOP) ++n_inclusions; 
          else ++n_overlaps; 
        }
      }
    }
    if (load_in_memory==false) { R[0] = NULL; if (r0!=NULL) delete r0; } 
    r0 = r;
    PRG.Check();
  }
  PRG.Done();
  fprintf(stderr, "* The file is sorted! Found %ld inclusions and %ld overlaps.\n", n_inclusions, n_overlaps);
}



//---------PrintBEDGraphFormat--------
//
void GenomicRegionSet::PrintBEDGraphFormat(char *title, char *color, char *position, bool convert_chromosome)
{
  if (n_regions==0) return;
  //NOTE: we should check on the fly if the regions are sorted
  Progress PRG("Printing in bedgraph format...",n_regions);
  if (strlen(position)>0) printf("browser position %s\n", position);
  if (strlen(title)>0) printf("track type=bedGraph name='%s' color=%s visibility=1\n", title, color);
  for (GenomicRegion *r=Get(); r!=NULL; r=Next()) {
    if (r->I.front()->STRAND=='-') r->PrintError("only forward strand accepted for this operation!\n"); 
    PrintChromosome(r->I.front()->CHROMOSOME,convert_chromosome);
    printf(" %ld %ld %s\n", r->I.front()->START, r->I.front()->START, r->LABEL);
    PRG.Check();
  }
  PRG.Done();
}



//---------------------------------------------------------------------------------------------//
// END CLASS: GenomicRegionSet                                                                 //
//---------------------------------------------------------------------------------------------//
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//
















  
  

//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//
//---------------------------------------------------------------------------------------------//
// CLASS: GenomicRegionSetScanner                                                              //
//---------------------------------------------------------------------------------------------//

//---------Constructor--------
//
GenomicRegionSetScanner::GenomicRegionSetScanner(GenomicRegionSet *R, StringLIntMap *bounds, long int win_step, long int win_size, bool use_labels_as_values, bool ignore_strand, char preprocess)
{
  if (R->format=="SEQ") { cerr << "Error: this operation does not accept SEQ format!\n"; exit(1); }
  if (R->n_regions==0) return;
  if (bounds==NULL) { cerr << "Error: this operation requires genomic bounds!\n"; exit(1); }

  this->R = R;
  this->bounds = bounds;
  this->win_step = win_step;
  this->win_size = win_size;
  this->use_labels_as_values = use_labels_as_values;
  this->ignore_strand = ignore_strand;
  this->preprocess = preprocess;
  if (win_size%win_step!=0) { cerr << "Error: window size must be a multiple of window step in 'GenomicRegionSetScanner'!\n"; exit(1); }
  n_win_combine = win_size/win_step;

}


//---------Destructor--------
//
GenomicRegionSetScanner::~GenomicRegionSetScanner()
{

}


//---------------------------------------------------------------------------------------------//
// END CLASS: GenomicRegionSetScanner                                                          //
//---------------------------------------------------------------------------------------------//
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//












//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//
//---------------------------------------------------------------------------------------------//
// CLASS: SortedGenomicRegionSetScanner                                                        //
//---------------------------------------------------------------------------------------------//

//---------Constructor--------
//
SortedGenomicRegionSetScanner::SortedGenomicRegionSetScanner(GenomicRegionSet *R, StringLIntMap *bounds, long int win_step, long int win_size, bool use_labels_as_values, bool ignore_strand, char preprocess)
  : GenomicRegionSetScanner(R,bounds,win_step,win_size,use_labels_as_values,ignore_strand,preprocess)
{
  chr = bounds->begin();
  r = R->Get();
  strand = '+';
  v = new long int[n_win_combine];
  v_sum = -1;
}


//---------Destructor--------
//
SortedGenomicRegionSetScanner::~SortedGenomicRegionSetScanner()
{
  delete v;
}


//---------Test--------
//
bool SortedGenomicRegionSetScanner::Test()
{
  int t = strcmp(chr->first.c_str(),r->I.front()->CHROMOSOME);
  return (t>0)||((t==0)&&(strand>r->I.front()->STRAND));
}


//---------PrintInterval--------
//
void SortedGenomicRegionSetScanner::PrintInterval(FILE *out_file)
{
  fprintf(out_file, "%s %c %ld %ld", chr->first.c_str(), strand, stop-win_size+1, stop);
}


//---------Next--------
//
long int SortedGenomicRegionSetScanner::Next()
{
  if (v_sum>=0) goto next;
  for (; chr!=bounds->end(); chr++,strand='+') {
    while (strand!=' ') {
      for (long int j=0; j<n_win_combine; j++) v[j] = 0;
      while ((r!=NULL)&&Test()) r = R->Next(!ignore_strand,false);
      for (k=0,v_sum=0,start=1,stop=start+win_step-1; stop<=chr->second; start+=win_step,stop+=win_step,k=(k+1)%n_win_combine) {
        v_sum -= v[k];
        for (v[k]=0; (r!=NULL)&&(strcmp(r->I.front()->CHROMOSOME,chr->first.c_str())==0)&&((ignore_strand==true)||(r->I.front()->STRAND==strand))&&(r->I.front()->START<=stop); ) {
          if (r->I.size()!=1) r->PrintError("single-interval regions expected for this operation!\n"); 
          if (preprocess=='p') {
            if (r->I.front()->STOP<=stop) { v[k] += r->I.front()->STOP-r->I.front()->START+1; r = R->Next(!ignore_strand,false); }
            else { v[k] += stop-r->I.front()->START+1; r->I.front()->START = stop+1; }
          }
          else if (preprocess=='1') {
            if (use_labels_as_values) v[k] += max(1L,atol(r->LABEL));
            else v[k]++;
          }
          else { fprintf(stderr, "Error: [SortedGenomicRegionSetScanner] preprocess operator '%c' not supported!\n", preprocess); exit(1); }
          r = R->Next(!ignore_strand,false);
        }
        v_sum += v[k];
        if (stop>=win_size) return v_sum;
next:
        while (0) {};
      }
      if ((strand=='+')&&(ignore_strand==false)) strand = '-';
      else strand = ' ';
    }
  }
  return -1;
}


//---------Next--------
//
long int SortedGenomicRegionSetScanner::Next(GenomicRegionSet *Ref)
{
  if (Ref==NULL) return Next();
  GenomicRegion *q = Ref->Get();
  bool sorted_by_strand = !ignore_strand;
  while (q!=NULL) {
    long int c = Next();
    if (c==-1) return -1;
    GenomicInterval w(chr->first.c_str(),strand,stop-win_size+1,stop);
	while (q!=NULL) {
	  int d = q->I.front()->CalcDirection(&w,sorted_by_strand);
	  if (d<0) q = Ref->Next(sorted_by_strand,false);
	  else if (d==0) return c;
	  else if (d>0) break;
	}
  }
  return -1;
}


//---------Next--------
//
long int SortedGenomicRegionSetScanner::Next(GenomicRegionSetIndex *index)
{
  if (index==NULL) return Next();
  while (true) {
    long int c = Next();
    if (c==-1) return -1;
    GenomicInterval w(chr->first.c_str(),strand,stop-win_size+1,stop);
    if (index->GetOverlap(&w,false,ignore_strand)!=NULL) return c;
  }
}


//---------------------------------------------------------------------------------------------//
// END CLASS: SortedGenomicRegionSetScanner                                                    //
//---------------------------------------------------------------------------------------------//
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//












//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//
//---------------------------------------------------------------------------------------------//
// CLASS: UnsortedGenomicRegionSetScanner                                                      //
//---------------------------------------------------------------------------------------------//

//---------Constructor--------
//
UnsortedGenomicRegionSetScanner::UnsortedGenomicRegionSetScanner(GenomicRegionSet *R, StringLIntMap *bounds, long int win_step, long int win_size, bool use_labels_as_values, bool ignore_strand, char preprocess)
  : GenomicRegionSetScanner(R,bounds,win_step,win_size,use_labels_as_values,ignore_strand,preprocess)
{
	// initialize window counts per chromosome
	Progress PRG1("Initializing structures...",bounds->size());
	for (StringLIntMap::iterator p=bounds->begin(); p!=bounds->end(); p++) {
		unsigned long int n = p->second/win_step;				// number of micro-windows
		count[p->first] = new unsigned long int*[2];
		unsigned long int *v_fwd = count[p->first][0] = new unsigned long int[n+1];
		unsigned long int *v_rev = count[p->first][1] = ignore_strand==true?NULL:new unsigned long int[n+1];
		v_fwd[0] = n;
		if (ignore_strand==false) v_rev[0] = n;
		for (unsigned long int k=1; k<=n; k++) { v_fwd[k] = 0; if (ignore_strand==false) v_rev[k] = 0; }
		PRG1.Check();
	}
	PRG1.Done();

	// scan genomic region set and determine counts
	Progress PRG2("Determining counts per window...",R->n_regions);
	for (GenomicRegion *r=R->Get(); r!=NULL; r=R->Next()) {
	    for (GenomicIntervalSet::iterator i=r->I.begin(); i!=r->I.end(); i++) {
	    	if (((*i)->START>(*i)->STOP)||((*i)->STOP<=0)) continue;
			CountMap::iterator p=count.find((*i)->CHROMOSOME);
			if (p!=count.end()) {
				long int start;
				if (preprocess=='1') start = (*i)->START;
				else if (preprocess=='c') start = (*i)->START+((*i)->STOP-(*i)->START)/2;
				else { fprintf(stderr, "Error: [UnsortedGenomicRegionSetScanner] preprocess operator '%c' not supported!\n", preprocess); exit(1); }
				long int w = (start-1)/win_step+1;
				unsigned long int *v = ((ignore_strand==true)||((*i)->STRAND=='+'))?p->second[0]:p->second[1];
				if ((start>=1)&&(w<=(long int)v[0])) {
					if (use_labels_as_values) v[w] += max(1L,atol(r->LABEL));
					else v[w]++;
				}
			}
	    }
	    PRG2.Check();
	}
	PRG2.Done();

	// sum up counts on sliding windows comprising n_win_combine micro-windows
	Progress PRG3("Summing counts in sliding windows...",count.size());
	for (CountMap::iterator p=count.begin(); p!=count.end(); p++) {
		for (int z=0; z<=(ignore_strand==true?0:1); z++) {
			unsigned long int *v = p->second[z];
			if (v[0]<(unsigned long int)n_win_combine) v[0] = 0;
			else {
				unsigned long int c, sum = 0;
				for (int k=1; k<=n_win_combine-1; k++) sum += v[k];
				v[0] = v[0] - n_win_combine + 1;
				for (unsigned long int k=1; k<=v[0]; k++) {
					sum = sum + v[k+n_win_combine-1];
					c = v[k];
					v[k] = sum;
					sum = sum - c;
				}
			}
		}
		PRG3.Check();
	}
	PRG3.Done();

	// initialize current window information
	Init();
}


//---------Destructor--------
//
UnsortedGenomicRegionSetScanner::~UnsortedGenomicRegionSetScanner()
{
	for (CountMap::iterator p=count.begin(); p!=count.end(); p++) {
		if (p->second[0]!=NULL) delete p->second[0];
		if (p->second[1]!=NULL) delete p->second[1];
		delete p->second;
	}
}


//---------SetCurrent--------
//
void UnsortedGenomicRegionSetScanner::Init()
{
	current_chr = count.begin();
	current_strand = '+';
	current_v = current_chr->second[0];
	current_n = current_chr!=count.end()?current_v[0]:0;
	current_win = 0;
}


//---------PrintInterval--------
//
void UnsortedGenomicRegionSetScanner::PrintInterval(FILE *out_file)
{
	fprintf(out_file, "%s %c %ld %ld", current_chr->first.c_str(), current_strand, win_step*(current_win-1)+1, win_step*(current_win-1)+win_size);
}


//---------Next--------
//
long int UnsortedGenomicRegionSetScanner::Next()
{
	if (current_chr==count.end()) return -1;
	current_win++;
	if (current_win>current_n) {
		if ((ignore_strand==true)||(current_strand=='-')) {
			current_chr++;
			if (current_chr==count.end()) return -1;
			current_strand = '+';
		}
		else current_strand = '-';
		current_v = current_strand=='+'?current_chr->second[0]:current_chr->second[1];
		current_n = current_chr!=count.end()?current_v[0]:0;
		current_win = 1;
	}
	return current_v[current_win];
}


//---------Next--------
//
long int UnsortedGenomicRegionSetScanner::Next(GenomicRegionSet *Ref)
{
	if (Ref==NULL) return Next();

	GenomicRegion *q = Ref->Get();
	bool sorted_by_strand = !ignore_strand;
	while (q!=NULL) {
		long int c = Next();
		if (c==-1) return -1;
		GenomicInterval w(current_chr->first.c_str(),current_strand,win_step*(current_win-1)+1,win_step*(current_win-1)+win_size);
		while (q!=NULL) {
		  int d = q->I.front()->CalcDirection(&w,sorted_by_strand);
		  if (d<0) q = Ref->Next(sorted_by_strand,false);
		  else if (d==0) return c;
		  else if (d>0) break;
		}
	}
	return -1;
}


//---------Next--------
//
long int UnsortedGenomicRegionSetScanner::Next(GenomicRegionSetIndex *index)
{
  if (index==NULL) return Next();
  while (true) {
    long int c = Next();
    if (c==-1) return -1;
	GenomicInterval w(current_chr->first.c_str(),current_strand,win_step*(current_win-1)+1,win_step*(current_win-1)+win_size);
    if (index->GetOverlap(&w,false,ignore_strand)!=NULL) return c;
  }
}


//---------------------------------------------------------------------------------------------//
// END CLASS: UnsortedGenomicRegionSetScanner                                                  //
//---------------------------------------------------------------------------------------------//
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//












//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//
//---------------------------------------------------------------------------------------------//
// CLASS: GenomicRegionSetOverlaps                                                             //
//---------------------------------------------------------------------------------------------//


//---------Constructor--------
//
GenomicRegionSetOverlaps::GenomicRegionSetOverlaps(GenomicRegionSet *QuerySet, GenomicRegionSet *IndexSet)
{
  this->QuerySet = QuerySet;
  this->IndexSet = IndexSet;
}



//---------Destructor--------
//
GenomicRegionSetOverlaps::~GenomicRegionSetOverlaps()
{

}



//---------GetOverlap--------
//
GenomicRegion *GenomicRegionSetOverlaps::GetOverlap(bool match_gaps, bool ignore_strand)
{
  for (GenomicRegion *r=GetMatch(); r!=NULL; r=NextMatch()) {
    if (match_gaps||current_qreg->OverlapsWith(r,ignore_strand)) {
      if (ignore_strand) return r;
      else if (current_qreg->I.front()->STRAND==r->I.front()->STRAND) return r;
    }
  }
  return NULL;
}



//---------NextOverlap--------
//
GenomicRegion *GenomicRegionSetOverlaps::NextOverlap(bool match_gaps, bool ignore_strand)
{
  for (GenomicRegion *r=NextMatch(); r!=NULL; r=NextMatch()) {
    if (match_gaps||current_qreg->OverlapsWith(r,ignore_strand)) {
      if (ignore_strand) return r;
      else if (current_qreg->I.front()->STRAND==r->I.front()->STRAND) return r;
    }
  }
  return NULL;
}



//---------CalcQueryCoverage--------
//
unsigned long int GenomicRegionSetOverlaps::CalcQueryCoverage(bool match_gaps, bool ignore_strand, bool use_labels_as_values)
{
  unsigned long int c = 0;
  for (GenomicRegion *r=GetOverlap(match_gaps,ignore_strand); r!=NULL; r=NextOverlap(match_gaps,ignore_strand)) {
    long int cc = match_gaps ? min(r->I.back()->STOP,current_qreg->I.back()->STOP)-max(r->I.front()->START,current_qreg->I.front()->START)+1 : current_qreg->CalcOverlap(r,ignore_strand);
    if (use_labels_as_values) cc *= atol(r->LABEL);
    c += cc;
  }
  return c;
}



//---------CalcIndexCoverage--------
//
unsigned long int *GenomicRegionSetOverlaps::CalcIndexCoverage(bool match_gaps, bool ignore_strand, bool use_labels_as_values)
{
  if (IndexSet->load_in_memory==false) { fprintf(stderr, "[GenomicRegionSetOverlaps::CalcIndexCoverage]: index set must be loaded in memory for this operation!\n"); exit(1); }
  Progress PRG("Processing queries...",QuerySet->n_regions);
  unsigned long int *coverage = new unsigned long int[IndexSet->n_regions];
  for (long int k=0; k<IndexSet->n_regions; k++) { IndexSet->R[k]->n_line = k; coverage[k] = 0; }
  for (GenomicRegion *qreg=GetQuery(); qreg!=NULL; qreg=NextQuery()) {
    for (GenomicRegion *ireg=GetOverlap(match_gaps,ignore_strand); ireg!=NULL; ireg=NextOverlap(match_gaps,ignore_strand)) {
      long int cc = match_gaps ? min(qreg->I.back()->STOP,ireg->I.back()->STOP)-max(qreg->I.front()->START,ireg->I.front()->START)+1 : ireg->CalcOverlap(qreg,ignore_strand);
      if (use_labels_as_values) cc *= atol(qreg->LABEL);
      coverage[ireg->n_line] += cc; 
    }
    PRG.Check();
  }
  PRG.Done();
  return coverage;
}




//---------CountQueryOverlaps--------
//
unsigned long int GenomicRegionSetOverlaps::CountQueryOverlaps(bool match_gaps, bool ignore_strand, bool use_labels_as_values)
{
  unsigned long int c = 0;
  for (GenomicRegion *r=GetOverlap(match_gaps,ignore_strand); r!=NULL; r=NextOverlap(match_gaps,ignore_strand)) c += use_labels_as_values?atol(r->LABEL):1;
  return c;
}




//---------CountIndexOverlaps--------
//
unsigned long int *GenomicRegionSetOverlaps::CountIndexOverlaps(bool match_gaps, bool ignore_strand, bool use_labels_as_values)
{
  if (IndexSet->load_in_memory==false) { fprintf(stderr, "[GenomicRegionSetOverlaps::CountIndexOverlaps]: index set must be loaded in memory for this operation!\n"); exit(1); }
  Progress PRG("Processing queries...",QuerySet->n_regions);
  unsigned long int *hits = new unsigned long int[IndexSet->n_regions];
  for (long int k=0; k<IndexSet->n_regions; k++) { IndexSet->R[k]->n_line = k; hits[k] = 0; }
  for (GenomicRegion *qreg=GetQuery(); qreg!=NULL; qreg=NextQuery()) {
    for (GenomicRegion *ireg=GetOverlap(match_gaps,ignore_strand); ireg!=NULL; ireg=NextOverlap(match_gaps,ignore_strand)) 
      hits[ireg->n_line] += use_labels_as_values?atol(qreg->LABEL):1; 
    PRG.Check();
  }
  PRG.Done();
  return hits;
}




//---------------------------------------------------------------------------------------------//
// END CLASS: GenomicRegionSetOverlaps                                                         //    
//---------------------------------------------------------------------------------------------//
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//
















//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//
//---------------------------------------------------------------------------------------------//
// CLASS: GenomicRegionSetIndex                                                                //
//---------------------------------------------------------------------------------------------//


//---------Constructor--------
//
GenomicRegionSetIndex::GenomicRegionSetIndex(GenomicRegionSet *regSet, const char *bin_bits)
{
  this->regSet = regSet;
  
  if (regSet->load_in_memory==false) { fprintf(stderr, "[GenomicRegionSetIndex] Error: index regions must be loaded in memory!\n"); exit(1); }

  // book-keeping for regions that fall in the same bin
  r_next = new long int[regSet->n_regions];
  for (long int k=0; k<regSet->n_regions; k++) r_next[k] = -1;

  // Calculate basic statistics
  map<string,long int> chrom_size;
  Progress PRG1("Checking index regions...",regSet->n_regions);
  for (long int k=0; k<regSet->n_regions; k++) {
    GenomicRegion *r = regSet->R[k];
    if (r->IsCompatibleSortedAndNonoverlapping()==false) r->PrintError("index regions should be compatible, sorted and non-overlapping!");
    long int start = r->I.front()->START;
    long int stop = r->I.back()->STOP;
    if ((start>stop)||(stop<=0)) continue; 	// ignore invalid intervals
    map<string,long int>::iterator it = chrom_size.find(r->I.front()->CHROMOSOME);
    if (it==chrom_size.end()) chrom_size[r->I.front()->CHROMOSOME] = stop;
    else it->second = max(it->second,stop);
    PRG1.Check();
  }
  PRG1.Done();

  // set bin parameters 
  if (bin_bits==NULL) {
    n_levels = 5;				// use Kent et al. (UCSC Genome Browser) settings as default
    n_bits = new int[n_levels];
    n_bits[0] = 17;
    n_bits[1] = 20;
    n_bits[2] = 23;
    n_bits[3] = 26;
    n_bits[n_levels-1] = 60;			// the last bin level should only have one bin, i think 60 bits will guarrantee that :-) 
  }
  else {
    char *p = StrCopy(bin_bits);
    n_levels = CountTokens(p,',')+1;
    n_bits = new int[n_levels];
    int l = 0;
    char *q = p; for (char *s=GetNextToken(&q,','); s[0]!=0; s=GetNextToken(&q,',')) n_bits[l++] = atoi(s);
    delete p;
    n_bits[n_levels-1] = 60;
  }

  // initialize bin structures for each chromosome
  Progress PRG1a("Initializing bin structures for each chromosome...",1); 
  for (map<string,long int>::iterator it=chrom_size.begin(); it!=chrom_size.end(); it++) {
    long int *chrom_n_bins = new long int[n_levels];
    long int **chrom_bins = new long int*[n_levels];
    index[it->first] = new BinSet(chrom_n_bins,chrom_bins);
    //cerr << "* chrom = " << it->first << "; size = " << it->second;
    for (int l=0; l<n_levels; l++) {
      chrom_n_bins[l] = (it->second>>n_bits[l])+1;
      //cerr << "; n_bins[" << l << "] = " << chrom_n_bins[l];
      chrom_bins[l] = new long int[chrom_n_bins[l]];
      for (long int b=0; b<chrom_n_bins[l]; b++) chrom_bins[l][b] = -1;
    }
    //cerr << '\n';
	PRG1a.Check();
  }
  PRG1a.Done();
  
  // process regions into the bins
  //double mean_bin_occupancy = 0;
  Progress PRG2("Creating index...",regSet->n_regions);
  for (long int k=0; k<regSet->n_regions; k++) {
    GenomicRegion *r = regSet->R[k];
    long int start = r->I.front()->START;
    long int stop = r->I.back()->STOP;
    if ((start>stop)||(stop<=0)) continue; 	// AddWarning("invalid interval (start>stop or stop<=0)!");
    if (start<=0) start = 1; 
    long int **chrom_bins = index[r->I.front()->CHROMOSOME]->second;
    for (int l=0; l<n_levels; l++) {
      long int b_start = start>>n_bits[l];
      long int b_stop = stop>>n_bits[l];
      if (b_start==b_stop) {
        long int z = chrom_bins[l][b_start]; 
        if (z!=-1) r_next[k] = z;
        chrom_bins[l][b_start] = k;
        //mean_bin_occupancy += min(1.0,(double)(r->I.back()->STOP-r->I.front()->START+1+100)/(1<<n_bits[l]));
        break;
      }
    }
    PRG2.Check();
  }
  //mean_bin_occupancy /= regSet->n_regions;
  PRG2.Done();
  //fprintf(stderr, "* mean bin occupancy = %.6f\n", mean_bin_occupancy);

  // double-check if all regions are stored in the bins
  /*
  Progress PRG3("Counting used...",index.size());
  long int n_used = 0;
  for (map<string,BinSet*>::iterator it=index.begin(); it!=index.end(); it++,PRG3.Check()) {
    long int *chrom_n_bins = it->second->first;
    long int **chrom_bins = it->second->second; 
    for (int l=0; l<n_levels; l++)
      for (long int b=0; b<chrom_n_bins[l]; b++) 
        for (long int z=chrom_bins[l][b]; z!=-1; z=r_next[z]) n_used++;
  }
  PRG3.Done();
  // The following test will fail in general because 'rogue' regions are not used!
  if (n_used!=regSet->n_regions) { fprintf(stderr, "Bug: [UnsortedGenomicRegionSetOverlaps] n_used should equal n_regions!\n"); exit(1); }
  */
}



//---------Destructor--------
//
GenomicRegionSetIndex::~GenomicRegionSetIndex()
{
  delete r_next;
  delete n_bits;
  for (map<string,BinSet*>::iterator it=index.begin(); it!=index.end(); it++) {
    delete it->second->first; 
    delete [] it->second->second;
  }
}



//---------GetMatch--------
//
GenomicRegion *GenomicRegionSetIndex::GetMatch(GenomicInterval *i)
{
  current_qint = i;
  map<string,BinSet*>::iterator it = index.find(i->CHROMOSOME);
  current_binset = it!=index.end()?it->second:NULL;
  new_query = true;
  return NextMatch();
}



//---------NextMatch--------
//
GenomicRegion *GenomicRegionSetIndex::NextMatch()
{
  if (current_binset==NULL) return (current_ireg=NULL); 
  static long int start, stop, l, b, b_stop, current_k;
  static long int *n_bins;
  if (new_query) {
    new_query = false;
    n_bins = current_binset->first;
    l = 0;
    start = current_qint->START;
    stop = current_qint->STOP;
    if (stop<=0) return (current_ireg=NULL);
    if (start>stop) return (current_ireg=NULL);
    if (start<=0) start = 1; //current_qint->PrintError("start position must be positive!");
    b = start>>n_bits[l];
    b_stop = min(stop>>n_bits[l],n_bins[l]-1);
    if (b>=n_bins[l]) return (current_ireg=NULL);  
    current_k = current_binset->second[l][b];
  }
  while (true) {
    while (current_k!=-1) {
      current_ireg = regSet->R[current_k];
      current_k = r_next[current_k];
      if ((start<=current_ireg->I.back()->STOP)&&(stop>=current_ireg->I.front()->START)) return current_ireg;  
    }
    b++;
    if (b>b_stop) { 
      l++;
      if (l>=n_levels) break; 
      b = start>>n_bits[l];
      b_stop = min(stop>>n_bits[l],n_bins[l]-1);
    }
    current_k = current_binset->second[l][b];
  }
  return (current_ireg=NULL);
}



//---------GetOverlap--------
//
GenomicRegion *GenomicRegionSetIndex::GetOverlap(GenomicInterval *i, bool match_gaps, bool ignore_strand)
{
  current_qint = i;
  for (GenomicRegion *r=GetMatch(i); r!=NULL; r = NextMatch()) {
    if (match_gaps||current_qint->OverlapsWith(r,ignore_strand)) {
      if (ignore_strand) return r;
      else if (current_qint->STRAND==r->I.front()->STRAND) return r;
    }
  }
  return NULL;
}



//---------NextOverlap--------
//
GenomicRegion *GenomicRegionSetIndex::NextOverlap(bool match_gaps, bool ignore_strand)
{
  for (GenomicRegion *r=NextMatch(); r!=NULL; r=NextMatch()) {
    if (match_gaps||current_qint->OverlapsWith(r,ignore_strand)) {
      if (ignore_strand) return r;
      else if (current_qint->STRAND==r->I.front()->STRAND) return r;
    }
  }
  return NULL;
}




//---------------------------------------------------------------------------------------------//
// END CLASS: GenomicRegionSetIndex                                                            //    
//---------------------------------------------------------------------------------------------//
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//





















//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//
//---------------------------------------------------------------------------------------------//
// CLASS: UnsortedGenomicRegionSetOverlaps                                                     //
//---------------------------------------------------------------------------------------------//


//---------Constructor--------
//
UnsortedGenomicRegionSetOverlaps::UnsortedGenomicRegionSetOverlaps(GenomicRegionSet *QuerySet, GenomicRegionSet *IndexSet, const char *bin_bits)
 : GenomicRegionSetOverlaps(QuerySet,IndexSet)
{
  if (IndexSet->load_in_memory==false) { fprintf(stderr, "Error: [UnsortedGenomicRegionSetOverlaps] index regions must be loaded in memory!\n"); exit(1); }

  // book-keeping for regions that fall in the same bin
  r_next = new long int[IndexSet->n_regions];
  for (long int k=0; k<IndexSet->n_regions; k++) r_next[k] = -1;

  // Calculate basic statistics
  map<string,long int> chrom_size;
  Progress PRG1("Checking index regions...",IndexSet->n_regions);
  for (long int k=0; k<IndexSet->n_regions; k++) {
    GenomicRegion *r = IndexSet->R[k];
    if (r->IsCompatibleSortedAndNonoverlapping()==false) r->PrintError("index regions should be compatible, sorted and non-overlapping!");
    long int start = r->I.front()->START;
    long int stop = r->I.back()->STOP;
    if ((start>stop)||(stop<=0)) continue; 
    map<string,long int>::iterator it = chrom_size.find(r->I.front()->CHROMOSOME);
    if (it==chrom_size.end()) chrom_size[r->I.front()->CHROMOSOME] = stop;
    else it->second = max(it->second,stop);
    PRG1.Check();
  }
  PRG1.Done();

  // set bin parameters 
  if (bin_bits==NULL) {
    n_levels = 5;				// use Kent et al. (UCSC Genome Browser) settings as default
    n_bits = new int[n_levels];
    n_bits[0] = 17;
    n_bits[1] = 20;
    n_bits[2] = 23;
    n_bits[3] = 26;
    n_bits[n_levels-1] = 60;			// the last bin level should only have one bin, i think 60 bits will guarrantee that :-) 
  }
  else {
    char *p = StrCopy(bin_bits);
    n_levels = CountTokens(p,',')+1;
    n_bits = new int[n_levels];
    int l = 0;
    char *q = p; for (char *s=GetNextToken(&q,','); s[0]!=0; s=GetNextToken(&q,',')) n_bits[l++] = atoi(s);
    delete p;
    n_bits[n_levels-1] = 60;
  }

  // initialize bin structures for each chromosome
  Progress PRG1a("Initializing bin structures for each chromosome...",1); 
  for (map<string,long int>::iterator it=chrom_size.begin(); it!=chrom_size.end(); it++) {
    long int *chrom_n_bins = new long int[n_levels];
    long int **chrom_bins = new long int*[n_levels];
    index[it->first] = new BinSet(chrom_n_bins,chrom_bins);
    for (int l=0; l<n_levels; l++) {
      chrom_n_bins[l] = (it->second>>n_bits[l])+1;
      chrom_bins[l] = new long int[chrom_n_bins[l]];
      for (long int b=0; b<chrom_n_bins[l]; b++) chrom_bins[l][b] = -1;
    }
	PRG1a.Check();
  }
  PRG1a.Done();
  
  // process regions into the bins
  Progress PRG2("Creating index...",IndexSet->n_regions);
  for (long int k=0; k<IndexSet->n_regions; k++) {
    GenomicRegion *r = IndexSet->R[k];
    long int start = r->I.front()->START;
    long int stop = r->I.back()->STOP;
    if ((start>stop)||(stop<=0)) continue; 	// AddWarning("invalid interval (start>stop or stop<=0)!");
    if (start<=0) start = 1; 
    long int **chrom_bins = index[r->I.front()->CHROMOSOME]->second;
    for (int l=0; l<n_levels; l++) {
      long int b_start = start>>n_bits[l];
      long int b_stop = stop>>n_bits[l];
      if (b_start==b_stop) {
        long int z = chrom_bins[l][b_start]; 
        if (z!=-1) r_next[k] = z;
        chrom_bins[l][b_start] = k;
        break;
      }
    }
    PRG2.Check();
  }
  PRG2.Done();
}



//---------Destructor--------
//
UnsortedGenomicRegionSetOverlaps::~UnsortedGenomicRegionSetOverlaps()
{
  delete r_next;
  delete n_bits;
  for (map<string,BinSet*>::iterator it=index.begin(); it!=index.end(); it++) {
    delete it->second->first; 
    delete [] it->second->second;
  }
}



//---------GetQuery--------
//
GenomicRegion *UnsortedGenomicRegionSetOverlaps::GetQuery()
{
  current_qreg = QuerySet->Get();
  if (current_qreg&&(current_qreg->IsCompatibleSortedAndNonoverlapping()==false)) current_qreg->PrintError("query regions should be compatible, sorted and non-overlapping!");
  return current_qreg;
}



//---------NextQuery--------
//
GenomicRegion *UnsortedGenomicRegionSetOverlaps::NextQuery()
{
  current_qreg = QuerySet->Next();
  if (current_qreg&&(current_qreg->IsCompatibleSortedAndNonoverlapping()==false)) current_qreg->PrintError("query regions should be compatible, sorted and non-overlapping!");
  return current_qreg;
}



//---------GetMatch--------
//
GenomicRegion *UnsortedGenomicRegionSetOverlaps::GetMatch()
{
  map<string,BinSet*>::iterator it = index.find(current_qreg->I.front()->CHROMOSOME);
  current_binset = it!=index.end()?it->second:NULL;
  new_query = true;
  return NextMatch();
}



//---------NextMatch--------
//
GenomicRegion *UnsortedGenomicRegionSetOverlaps::NextMatch()
{
  if (current_binset==NULL) return (current_ireg=NULL); 
  static long int start, stop, l, b, b_stop, current_k;
  static long int *n_bins;
  if (new_query) {
    new_query = false;
    n_bins = current_binset->first;
    l = 0;
    start = current_qreg->I.front()->START;
    stop = current_qreg->I.back()->STOP;
    if (stop<=0) current_qreg->PrintError("stop position must be positive!");
    if (start>stop) current_qreg->PrintError("start position cannot be greater than stop position!");
    if (start<=0) start = 1; //current_qreg->PrintError("start position must be positive!");
    b = start>>n_bits[l];
    b_stop = min(stop>>n_bits[l],n_bins[l]-1);
    if (b>=n_bins[l]) return (current_ireg=NULL);  
    current_k = current_binset->second[l][b];
  }
  while (true) {
    while (current_k!=-1) {
      current_ireg = IndexSet->R[current_k];
      current_k = r_next[current_k];
      if ((start<=current_ireg->I.back()->STOP)&&(stop>=current_ireg->I.front()->START)) return current_ireg;  
    }
    b++;
    if (b>b_stop) { 
      l++;
      if (l>=n_levels) break; 
      b = start>>n_bits[l];
      b_stop = min(stop>>n_bits[l],n_bins[l]-1);
    }
    current_k = current_binset->second[l][b];
  }
  return (current_ireg=NULL);
}



//---------Done--------
//
bool UnsortedGenomicRegionSetOverlaps::Done()
{
  return current_qreg==NULL;
}



//---------------------------------------------------------------------------------------------//
// END CLASS: UnsortedGenomicRegionSetOverlaps                                                 //    
//---------------------------------------------------------------------------------------------//
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//



















//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//
//---------------------------------------------------------------------------------------------//
// CLASS: SortedGenomicRegionSetOverlaps                                                       //
//---------------------------------------------------------------------------------------------//

//---------Constructor--------
//
SortedGenomicRegionSetOverlaps::SortedGenomicRegionSetOverlaps(GenomicRegionSet *QuerySet, GenomicRegionSet *IndexSet, bool sorted_by_strand)
 : GenomicRegionSetOverlaps(QuerySet,IndexSet)
{
  this->sorted_by_strand = sorted_by_strand;
  this->ireg_buffer_interval = NULL;
  this->max_ireg_buffer_size = 0;
  
  current_qreg = QuerySet->Get();
  current_ireg = IndexSet->Get();
}



//---------Destructor--------
//
SortedGenomicRegionSetOverlaps::~SortedGenomicRegionSetOverlaps()
{
  ClearIndexBuffer();
  if (_MESSAGES_) fprintf(stderr, "Max used memory = %lu regions.\n", max_ireg_buffer_size);
}



//---------ClearIndexBuffer--------
//
void SortedGenomicRegionSetOverlaps::ClearIndexBuffer()
{
  if (IndexSet->load_in_memory==false) for (GenomicRegionList::iterator it=IRegBuffer.begin(); it!=IRegBuffer.end(); it++) delete *it;
  IRegBuffer.clear();
  if (ireg_buffer_interval!=NULL) delete ireg_buffer_interval;
  ireg_buffer_interval = NULL;
}



//---------LoadIndexBuffer--------
//
void SortedGenomicRegionSetOverlaps::LoadIndexBuffer()
{
  if (current_qreg==NULL) return;
  if (IRegBuffer.size()>max_ireg_buffer_size) max_ireg_buffer_size = IRegBuffer.size();
  if ((IRegBuffer.empty()==false)&&(current_qreg->CalcDirection(ireg_buffer_interval,sorted_by_strand)>0)) ClearIndexBuffer();
  current_ireg = IndexSet->Get();
  if ((current_ireg==NULL)&&IRegBuffer.empty()) return;
  while (current_ireg!=NULL) {
    if (current_ireg->IsCompatibleSortedAndNonoverlapping()==false) current_ireg->PrintError("index regions should be compatible, sorted and non-overlapping!");
    int d = current_qreg->CalcDirection(current_ireg,sorted_by_strand);
    if (d<0) break;
    else if (d==0) { 
      if (IRegBuffer.empty()==true) {
        if (ireg_buffer_interval!=NULL) delete ireg_buffer_interval;
        ireg_buffer_interval = new GenomicInterval(current_ireg->I.front()->CHROMOSOME,current_ireg->I.front()->STRAND,current_ireg->I.front()->START,current_ireg->I.back()->STOP);
      }
      else {
        ireg_buffer_interval->START = min(ireg_buffer_interval->START,current_ireg->I.front()->START);
        ireg_buffer_interval->STOP = max(ireg_buffer_interval->STOP,current_ireg->I.back()->STOP);
      }
      IRegBuffer.push_back(current_ireg); 
    }
    GenomicRegion *ireg_prev = IndexSet->GetRetain();
    current_ireg = IndexSet->Next();
    if ((current_ireg!=NULL)&&(current_ireg->IsBefore(ireg_prev,sorted_by_strand)==true)) current_ireg->PrintError("index regions are not sorted (sorted-by-strand = " + (string)(sorted_by_strand?"true":"false") + ")!");
    if ((d>0)&&(IndexSet->load_in_memory==false)) delete ireg_prev;
  }
  IRegBufferIterator = IRegBuffer.begin();
  current_ireg = IRegBufferIterator==IRegBuffer.end()?NULL:*IRegBufferIterator;
}


//---------GetQuery--------
//
GenomicRegion *SortedGenomicRegionSetOverlaps::GetQuery()
{
  current_qreg = QuerySet->Get();
  if ((current_qreg!=NULL)&&(current_qreg->IsCompatibleSortedAndNonoverlapping()==false)) current_qreg->PrintError("query regions should be compatible, sorted and non-overlapping!");
  LoadIndexBuffer();
  return current_qreg;
}


//---------NextQuery--------
//
GenomicRegion *SortedGenomicRegionSetOverlaps::NextQuery()
{
  GenomicRegion *qreg_prev = QuerySet->GetRetain();
  current_qreg = QuerySet->Next();
  if ((current_qreg!=NULL)&&(current_qreg->IsCompatibleSortedAndNonoverlapping()==false)) current_qreg->PrintError("query regions should be compatible, sorted and non-overlapping!");
  if ((current_qreg!=NULL)&&(current_qreg->IsBefore(qreg_prev,sorted_by_strand)==true)) current_qreg->PrintError("query regions are not sorted (sorted-by-strand = " + (string)(sorted_by_strand?"true":"false") + ")!");
  if (QuerySet->load_in_memory==false) if (qreg_prev!=NULL) delete qreg_prev;
  LoadIndexBuffer();
  return current_qreg;
}


//---------GetMatch--------
//
GenomicRegion *SortedGenomicRegionSetOverlaps::GetMatch()
{
  if ((current_qreg==NULL)||(current_ireg==NULL)) return NULL; 
  while (true) {
    int d = current_qreg->CalcDirection(current_ireg,sorted_by_strand);
    if (d>0) { 
      if (IndexSet->load_in_memory==false) delete *IRegBufferIterator;
      IRegBufferIterator = IRegBuffer.erase(IRegBufferIterator); 
      if (IRegBufferIterator==IRegBuffer.end()) return NULL; 
      else current_ireg = *IRegBufferIterator; 
    }
    else if (d<0) return NULL;
    else return current_ireg;
  }
  return NULL;
}


//---------NextMatch--------
//
GenomicRegion *SortedGenomicRegionSetOverlaps::NextMatch()
{
  IRegBufferIterator++; 
  current_ireg = IRegBufferIterator==IRegBuffer.end()?NULL:*IRegBufferIterator;
  if (current_ireg==NULL) return NULL;
  return GetMatch();
}


//---------Done--------
//
bool SortedGenomicRegionSetOverlaps::Done()
{
  return (current_qreg==NULL)||((IndexSet->Get()==NULL)&&IRegBuffer.empty());
}


//---------------------------------------------------------------------------------------------//
// END CLASS: SortedGenomicRegionSetOverlaps                                                   //
//---------------------------------------------------------------------------------------------//
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//










//---------ProcessStrand-----------
//
char ProcessStrand(char *strand)
{
  if ((strcmp(strand,"1")==0)||(strcmp(strand,"+")==0)) return '+';
  else if ((strcmp(strand,"-1")==0)||(strcmp(strand,"-")==0)) return '-';
  else if (strcmp(strand,".")==0) return '.';
  else { cerr << "Error: invalid strand '" << strand << "'!\n"; exit(1); }
}



//---------PrintChromosome-----------
//
void PrintChromosome(char *chromosome, bool convert)
{
  if (convert) printf("chr%s", strcmp(chromosome,"MT")==0?"M":chromosome);
  else printf("%s", chromosome);
}



//---------RegCenter--------
//
GenomicRegion *RegCenter(GenomicRegion *r)
{
  if (r==NULL) return NULL;

  if (r->I.size()!=1) r->PrintError("single-interval regions expected for this operation!\n"); 

  long int start = r->I.front()->START+(r->I.front()->STOP-r->I.front()->START)/2;
  long int stop = r->I.front()->STOP-(r->I.front()->STOP-r->I.front()->START)/2;
  r->I.front()->START = start;
  r->I.front()->STOP = stop;

  return r;
}




//---------ReadBounds-----------
//
StringLIntMap *ReadBounds(char *genome_reg_file, bool verbose) 
{
  if ((genome_reg_file!=NULL)&&(strlen(genome_reg_file)>0)) {
    StringLIntMap *bounds = new StringLIntMap();
    GenomicRegionSet RegSet(genome_reg_file,10000,verbose,false,true);
    long int line = 1;
    for (GenomicRegion *r=RegSet.Get(); r!=NULL; r=RegSet.Next(),line++) {
      if (r->I.size()!=1) { cerr << "label = " << r->LABEL << '\n'; r->PrintError("genome regions should be single-interval regions!\n"); }
      string chr = r->I.front()->CHROMOSOME;
      if (bounds->find(chr)==bounds->end()) (*bounds)[chr] = r->I.front()->STOP;
      else if ((*bounds)[chr]!=r->I.front()->STOP) { cerr << "Error: chromosome " << chr << " has multiple lengths in genome file '" << genome_reg_file << "' line " << line << "!\n"; exit(1); }
    }
    return bounds;
  }
  else {
    cerr << "Error: genome region file is necessary for this operation!\n"; exit(1);
    return NULL;
  } 
}



//-------CalcBoundSize-----------
//
unsigned long int CalcBoundSize(StringLIntMap *bounds)
{
  unsigned long int y = 0;
  for (StringLIntMap::iterator x=bounds->begin(); x!=bounds->end(); x++) y += (unsigned long int)x->second; 
  return y;
}



//-------CalcRegSize-----------
//
unsigned long int CalcRegSize(char *reg_file)
{
  GenomicRegionSet RegSet(reg_file,100000,false,false,true);
  unsigned long int N = 0;
  Progress PRG("Calculating region set size...",1);
  for (GenomicRegion *r=RegSet.Get(); r!=NULL; r=RegSet.Next(),PRG.Check()) N += r->GetSize(true);
  PRG.Done();
  return N;
}


//---------CompareBinnedGenomicRegions-----------
//
bool CompareBinnedGenomicRegions(GenomicRegion *r, GenomicRegion *q)
{
  if (r->I.front()->START!=q->I.front()->START) return r->I.front()->START<q->I.front()->START;
  return r->I.back()->STOP>q->I.back()->STOP;
}



//---------CompareGenomicIntervals-----------
//
bool CompareGenomicIntervals(GenomicInterval *I, GenomicInterval *J)
{
  if (strcmp(I->CHROMOSOME,J->CHROMOSOME)!=0) return strcmp(I->CHROMOSOME,J->CHROMOSOME)<0;
  if (I->STRAND!=J->STRAND) return I->STRAND<J->STRAND;
  return I->START<J->START;
}



//---------SortGenomicIntervals-----------
//
void SortGenomicIntervals(GenomicIntervalSetAsList *L) 
{ 
  L->sort(CompareGenomicIntervals); 
}



//---------SortGenomicIntervals-----------
//
void SortGenomicIntervals(GenomicIntervalSetAsVector *V)
{ 
  sort(V->begin(),V->end(),CompareGenomicIntervals);
}



//---------SortGenomicIntervals-----------
//
void SortGenomicIntervals(GenomicIntervalSetAsArray *A)
{ 
  cerr << "Error: [SortGenomicIntervals] not implemented for arrays yet!\n";
  exit(1);
  //sort(V->begin(),V->end(),CompareGenomicIntervals);
}



//-------BinGenomicRegions-------------
//
GenomicRegionSetBins *BinGenomicRegions(GenomicRegionSet *rset, bool sorted_by_strand, int bin_bits)
{
	// initialize
	GenomicRegionSetBins *bins = new GenomicRegionSetBins();

	// determine chromosome limits
	if (rset->load_in_memory==false) { fprintf(stderr, "Error: [BinGenomicRegions] region set must be loaded in memory!\n"); exit(1); }
	Progress PRG0("[BinGenomicRegions] Determining chromosome limits...",rset->n_regions);
	for (GenomicRegion *r=rset->Begin(); r!=NULL; r=rset->Next()) {
		r->Sort();
		GenomicInterval *i = r->I.front();
		GenomicRegionSetBins::iterator p = bins->find(i->CHROMOSOME);
		if (p==bins->end()) {
			ChromosomeBins *b = new ChromosomeBins(i->START);
			bins->insert(pair<string,ChromosomeBins*>(i->CHROMOSOME,b));
		}
		else {
			p->second->min_start = min(p->second->min_start,i->START);
			p->second->max_start = max(p->second->max_start,i->START);
		}
		PRG0.Check();
	}
	PRG0.Done();

	// initialize window counts per chromosome
	Progress PRG1("[BinGenomicRegions] Initializing structures...",bins->size());
	for (GenomicRegionSetBins::iterator p=bins->begin(); p!=bins->end(); p++) {
		ChromosomeBins *chrom_bins = p->second;
		chrom_bins->n_bins = ((chrom_bins->max_start-chrom_bins->min_start)>>bin_bits) + 1;
		//cerr << p->first << '\t' << chrom_bins->min_start << ' ' << chrom_bins->max_start << ' ' << chrom_bins->n_bins << '\n';
		chrom_bins->b_fwd = new GenomicRegionList[chrom_bins->n_bins];
		chrom_bins->b_rev = sorted_by_strand==true?new GenomicRegionList[chrom_bins->n_bins]:NULL;
		PRG1.Check();
	}
	PRG1.Done();

	// assign genomic regions to bins
	Progress PRG2("[BinGenomicRegions] Assigning genomic regions to bins...",rset->n_regions);
	for (GenomicRegion *r=rset->Begin(); r!=NULL; r=rset->Next()) {
		GenomicInterval *i = r->I.front();
    	GenomicRegionSetBins::iterator p = bins->find(i->CHROMOSOME);
		if (p!=bins->end()) {
			ChromosomeBins *chrom_bins = p->second;
			unsigned long int w = (i->START-chrom_bins->min_start)>>bin_bits;
			GenomicRegionList *b = ((sorted_by_strand==false)||(i->STRAND=='+'))?chrom_bins->b_fwd:chrom_bins->b_rev;
			if ((w<0)||(w>=chrom_bins->n_bins)) { fprintf(stderr, "Error: [BinGenomicRegions] this must be a bug, please contact atsirigo@us.ibm.com.\n"); exit(1); }
			b[w].push_back(r);
	    }
	    PRG2.Check();
	}
	PRG2.Done();


	return bins;
}



//---------CalcOffsetsWithoutGaps-----------
//
long int *GetGapSizes(GenomicRegion *r, char *offset_op)
{
	long int *gap_size = new long int[r->I.size()];
	char strand = r->I.front()->STRAND;
	if ((strcmp(offset_op,"1")==0)||((strand=='+')&&(strcmp(offset_op,"5p")==0))||((strand=='-')&&(strcmp(offset_op,"3p")==0))) {
		gap_size[0] = 0;
		int k = 1;
		GenomicIntervalSet::iterator ref_int_prev = r->I.begin();
		for (GenomicIntervalSet::iterator ref_int=ref_int_prev+1; ref_int!=r->I.end(); ref_int_prev=ref_int,ref_int++,k++)
			gap_size[k] = gap_size[k-1] + (*ref_int)->START - (*ref_int_prev)->STOP - 1;
	}
	else {
		int k = r->I.size()-1;
		gap_size[k--] = 0;
		GenomicIntervalSet::reverse_iterator ref_int_prev = r->I.rbegin();
		for (GenomicIntervalSet::reverse_iterator ref_int=ref_int_prev+1; ref_int!=r->I.rend(); ref_int_prev++,ref_int++,k--)
			gap_size[k] = gap_size[k+1] + (*ref_int_prev)->START - (*ref_int)->STOP - 1;
	}
	return gap_size;
}



//---------CalcOffsetsWithoutGaps-----------
//
OffsetList *CalcOffsetsWithoutGaps(GenomicRegion *query_reg, GenomicRegion *ref_reg, char *offset_op, bool ignore_strand)
{
    if (query_reg->IsCompatibleSortedAndNonoverlapping()==false) { fprintf(stderr, "Warning: [CalcOffsetsWithoutGaps] query regions should be compatible, sorted and non-overlapping!"); return NULL; }
    if (ref_reg->IsCompatibleSortedAndNonoverlapping()==false) { fprintf(stderr, "Warning: [CalcOffsetsWithoutGaps] reference regions should be compatible, sorted and non-overlapping!"); return NULL; }
	OffsetList *offsets = new OffsetList();
	long int *gap_size = GetGapSizes(ref_reg,offset_op);
	int k = 0;
	for (GenomicIntervalSet::iterator ref_int=ref_reg->I.begin(); ref_int!=ref_reg->I.end(); ref_int++,k++) {
		for (GenomicIntervalSet::iterator query_int=query_reg->I.begin(); query_int!=query_reg->I.end(); query_int++) {
			if ((*query_int)->IsContained(*ref_int,ignore_strand)) {
				long int start_offset, stop_offset;
				(*query_int)->GetOffsetFrom(ref_reg,offset_op,ignore_strand,&start_offset,&stop_offset);
				offsets->push_back(pair<long int,long int>(start_offset-gap_size[k],stop_offset-gap_size[k]));
			}
		}
	}
	delete gap_size;
	return offsets;
}







