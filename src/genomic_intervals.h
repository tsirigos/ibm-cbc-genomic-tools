//
// Copyright (c) 2011 IBM Corporation. 
// All rights reserved. This program and the accompanying materials are made available under the terms of the GNU General Public License v2.0 
// which accompanies this distribution, and is available at http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
//

const string VERSION = "genomic_tools 2.1.0 [under development; revision 3]";


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


using namespace std;

// forward declarations
class Chromosomes;
class GenomicInterval;
typedef vector<GenomicInterval *> GenomicIntervalSetAsVector;
typedef list<GenomicInterval *> GenomicIntervalSetAsList;
class GenomicIntervalSetAsArray;
class GenomicRegion;
class GenomicRegionSet;

// types
typedef map<string,string> StringMap;
typedef map<string,long int> StringLIntMap;
typedef map<string,vector<long int> > StringVecLIntMap;

// choose container type for interval set (vector, list or array)
typedef GenomicIntervalSetAsVector GenomicIntervalSet;
//typedef GenomicIntervalSetAsList GenomicIntervalSet;
//typedef GenomicIntervalSetAsArray GenomicIntervalSet;

// for buffer in SortedGenomicRegionSetOverlaps
typedef list<GenomicRegion *> GenomicRegionList;


// constants
#define ALIGNMENT_GAP '*'					//! character for alignment gaps 
const char refCIGARop[] = "MD=X";			//! set of CIGAR operations that map to the reference sequence (SAM format)
const char fragCIGARop[] = "MIS=X";			//! set of CIGAR operations that map to the fragment sequence (SAM format)



  
  
//---------------------------------------------------------------------------------------------//
// CLASS DECLARATION: GenomicInterval                                                          //
//---------------------------------------------------------------------------------------------//
//!  This class implements the genomic interval.
/*!
  Genomic intervals can be constructed from strings that follow this format:

  CHROMOSOME \a \<SPACE\> STRAND \a \<SPACE\> START \a \<SPACE\> STOP

  Example:
  \code
    GenomicInterval x("chr1",'+',132034,135932);
    GenomicInterval y("chr1 + 135832 140102");
    cout << "size of interval x = " << x.GetSize() << '\n';
    cout << "size of interval y = " << y.GetSize() << '\n';
    cout << "overlap of x and y = " << x.GetOverlap(&y) << '\n';    
  \endcode
*/
//---------------------------------------------------------------------------------------------//
class GenomicInterval
{
 public:
  //! Class constructor, see example below.
  GenomicInterval(char *chromosome, char strand, long int start, long int stop, unsigned long int n_line=0);

  

  //! Class constructor, see example below.
  GenomicInterval(string chromosome, char strand, long int start, long int stop, unsigned long int n_line=0);

  

  //! Class constructor, see example below.
  GenomicInterval(char *inp, char delimiter=' ', unsigned long int n_line=0);

  

  //! Class constructor, see example below.
  GenomicInterval(GenomicInterval *I);

  
  
  //! Class Destructor.
  ~GenomicInterval();

  

  //------------------------------------------------------//
  //   Read & Print                                       //
  //------------------------------------------------------//


  //! Print interval (do not print newline character)
  void PrintInterval();

  

  //! Print interval in file (do not print newline character)
  void PrintInterval(FILE *file_ptr);

  

  //! Print interval with modified start and stop positions (do not print newline character)
  void PrintInterval(long int start, long int stop);


  //! Print interval with a newline character at the end. If <b>point</b> is 'true', then print every point in the interval separately. 
  void Print(bool point=false);

  

  //! Print interval with a label in front and newline character at the end. 
  void PrintWithLabel(char *label, bool point=false);

  
  
  //------------------------------------------------------//
  //   Check                                              //
  //------------------------------------------------------//


  //! Check integrity of interval start/stop positions (start>=1 and stop>=start)
  bool IsValid();

  

  //! Same as \ref GenomicInterval::IsValid, but it also prints error messages
  bool CheckValid(bool quiet=false);

  

  //! Compares with input interval <b>I</b> and returns true if its order is before <b>I</b>. The order is determined first by chromosome, then by strand, and finally by start position.
  bool IsBefore(GenomicInterval *I, bool sorted_by_strand);



  //! Compares with input interval <b>I</b> and returns true if its position order is before <b>I</b>. The position order is determined first by start position and then by stop position.
  bool IsPosBefore(GenomicInterval *I);



  //! Compares with input interval <b>I</b> and returns true if intervals are in the same chromosome and strand. 
  bool IsCompatibleWith(GenomicInterval *I, bool ignore_strand);



  //! Returns 'true' if it overlaps with input interval <b>I</b>.
  bool OverlapsWith(GenomicInterval *I, bool ignore_strand=false);



  //------------------------------------------------------//
  //   Set & Get                                          //
  //------------------------------------------------------//


  //! Set start/stop positions to new values
  void SetCoordinates(long int start, long int stop);

  

  //! If: op='1' returns start position, op='2' returns stop position, op='5p' returns 5-prime, op='3p' returns 3-prime. 
  long int GetCoordinate(char *op);
  


  //! Returns interval size
  size_t GetSize();
  


  //! Extracts genomic interval sequence from chromosomal sequences supplied in <b>C</b>. If parameter <b>replace</b> is true, it replaces "N" with lowercase "a". 
  char *GetSeq(Chromosomes *C, bool replace=false);
  


  //! Returns sequence length excluding 'N' characters.
  size_t GetSeqLength(Chromosomes *C);
  


  //------------------------------------------------------//
  //   Operations                                         //
  //------------------------------------------------------//


  //! Calculates the overlap with input interval <b>I</b>.
  long int CalcOverlap(GenomicInterval *I, bool ignore_strand);



  //! Calculates the distance from input interval <b>I</b> based on position operators <b>op</b> and <b>I_op</b>.
  long int CalcDistanceFrom(GenomicInterval *I, char *op, char *I_op);
  


  //! Returns -1 if before interval <b>i</b>, +1 if after, and 0 if the intervals overlap. The order is determined first by chromosome, then (optionally) by strand and finally by start position.
  int CalcDirection(GenomicInterval *i, bool sorted_by_strand);
  


  //! Prints genomic interval sequence from chromosomal sequences supplied in <b>C</b>. If parameter <b>replace</b> is true, it replaces "N" with lowercase "a". 
  void PrintSeq(Chromosomes *C, bool replace=false);
  


  //! If interval strand is '+', it is changed to '-' and vice versa. 
  void ReverseStrand();	
  


  //! Shifts interval start/stop positions. First it chooses reference position, and then shifts start and stop positions separately
  /*!
    \param start_shift		determines the shift of the start position
    \param stop_shift		determines the shift of the stop position
    \param strand_aware		if 'true', then start=5-prime and stop=3-prime
  */
  void ShiftPos(long int start_shift, long int stop_shift, bool strand_aware);
  


  //! Modifies interval start/stop positions. First it chooses which position to keep fixed, and then it shifts the non-fixed position. 
  /*!
    \param position_op 		selects fixed position as follows: '1'=start position, 'c'=center position, '5p'=5-prime, '3p'=3-prime
    \param position_shift	determines the shift of the non-fixed position (in the 'c' case, both positions are shifted symmetrically around the center)
  */
  void ModifyPos(char *position_op, long int position_shift);           



  //! Reverses interval start/stop positions with respect to chromosomal bounds, as if the reference point were in the end of the chromosome as opposed to the beginning. 
  void ReversePos(StringLIntMap *bounds);



  //! Corrects interval start/stop positions so as to comply with chromosomal bounds
  bool ApplyBounds(StringLIntMap *bounds);



  //! Prints in BED format. If <b>convert_chromosome</b> is 'true', then convert from ENSEMBL to UCSC names (by adding 'chr' prefix to each chromosome name, and by converting 'MT' to 'chrM'). 
  void PrintBEDFormat(char *color, bool convert_chromosome);



  //! Randomizes interval position within chromosomal bounds
  void Randomize(gsl_rng *random_generator, StringLIntMap *bounds);



  //! Computes start/stop offset distances from <b>ReferenceI</b>. 
  /*!
    \param ReferenceI 		pointer to reference interval
    \param op			selects reference point: '1'=start position, '2'=stop position, '5p'=5-prime position, '3p'=3-prime position
    \param ignore_strand	ignore mismatch between test and reference strand
    \param start_offset 	return value: the start offset
    \param stop_offset 		return value: the stop offset
  */
  void GetOffsetFrom(GenomicInterval *ReferenceI, char *op, bool ignore_strand, long int *start_offset, long int *stop_offset);


  

  //------------------------------------------------------//
  //   Data                                               //
  //------------------------------------------------------//
  
  unsigned long int n_line;		//!< file line number 
  char *CHROMOSOME;				//!< chromosome name (space and tab characters are not allowed)
  char STRAND;					//!< orientation ('+' or '-')
  long int START;				//!< interval start position in the chromosome (>=1)
  long int STOP;				//!< interval stop position in the chromosome (<=chromosome_size)
};

//---------------------------------------------------------------------------------------------//
// END CLASS DECLARATION: GenomicIntervals                                                     //
//---------------------------------------------------------------------------------------------//








//---------------------------------------------------------------------------------------------//
// CLASS: Chromosomes                                                                          //
//---------------------------------------------------------------------------------------------//
//!  This class is used to store chromosome DNA sequences. 
//---------------------------------------------------------------------------------------------//
class Chromosomes
{
 public:

  typedef map<string,pair<char*,size_t> > ChromosomeSeqData;		//!< for every chromosome, store a pointer to its sequence, and its size. 

  //! Class constructor.
  /*!
    \param chrom_map_dir 		the directory where the chromosome map file is located
    \param chrom_map_name 		the name of the chromosome map file

  Chromosome sequences are loaded on the fly from the *.dna files. For efficient sequence extraction, genomic intervals must be grouped by chromosome. 
  If not, use the alternative constructor described below. 

  Example of map file (stored in \ref map_name):

	1 \a \<TAB\> Mus_musculus.NCBIM37.52.dna.chromosome.1.fa.dna		\n
	10 \a \<TAB\> Mus_musculus.NCBIM37.52.dna.chromosome.10.fa.dna	\n
	11 \a \<TAB\> Mus_musculus.NCBIM37.52.dna.chromosome.11.fa.dna	\n
	12 \a \<TAB\> Mus_musculus.NCBIM37.52.dna.chromosome.12.fa.dna	\n
	13 \a \<TAB\> Mus_musculus.NCBIM37.52.dna.chromosome.13.fa.dna	\n
	14 \a \<TAB\> Mus_musculus.NCBIM37.52.dna.chromosome.14.fa.dna	\n
	15 \a \<TAB\> Mus_musculus.NCBIM37.52.dna.chromosome.15.fa.dna	\n
	16 \a \<TAB\> Mus_musculus.NCBIM37.52.dna.chromosome.16.fa.dna	\n
	17 \a \<TAB\> Mus_musculus.NCBIM37.52.dna.chromosome.17.fa.dna	\n
	18 \a \<TAB\> Mus_musculus.NCBIM37.52.dna.chromosome.18.fa.dna	\n
	19 \a \<TAB\> Mus_musculus.NCBIM37.52.dna.chromosome.19.fa.dna	\n
	2 \a \<TAB\> Mus_musculus.NCBIM37.52.dna.chromosome.2.fa.dna	\n
	3 \a \<TAB\> Mus_musculus.NCBIM37.52.dna.chromosome.3.fa.dna	\n
	4 \a \<TAB\> Mus_musculus.NCBIM37.52.dna.chromosome.4.fa.dna	\n
	5 \a \<TAB\> Mus_musculus.NCBIM37.52.dna.chromosome.5.fa.dna	\n
	6 \a \<TAB\> Mus_musculus.NCBIM37.52.dna.chromosome.6.fa.dna	\n
	7 \a \<TAB\> Mus_musculus.NCBIM37.52.dna.chromosome.7.fa.dna	\n
	8 \a \<TAB\> Mus_musculus.NCBIM37.52.dna.chromosome.8.fa.dna	\n
	9 \a \<TAB\> Mus_musculus.NCBIM37.52.dna.chromosome.9.fa.dna	\n
	MT \a \<TAB\> Mus_musculus.NCBIM37.52.dna.chromosome.MT.fa.dna	\n
	X \a \<TAB\> Mus_musculus.NCBIM37.52.dna.chromosome.X.fa.dna	\n
	Y \a \<TAB\> Mus_musculus.NCBIM37.52.dna.chromosome.Y.fa.dna	\n

  Each *.dna file contains the chromosome sequence in one line. 
  */
  Chromosomes(char *chrom_map_dir, char *chrom_map_name);


  //! Class constructor.
  /*!
    \param fasta_file 	FASTA file containing chromosome sequences
    \param verbose 		if 'true', print messages while loading chromosome sequences

   Chromosome sequences are loaded in memory from the FASTA file. Example FASTA file:

   >chromosome_1										\n
   TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTA	\n
   >chromosome_2										\n
   CGTATCCCACACACCACACCCACACACCACACCCACACACACCCACACCC	\n
   >chromosome_3										\n
   TCCTCCCAGAATCTGGAGAGGTCAACCTGTTCTTCAAAGCAGTGGTGGAT	\n
   >chromosome_4										\n
   ACCCTAACCCTAACCCTAACCCTAACCCTAACCCTACCCTAACCCTAACC	\n
   
  */
  Chromosomes(char *fasta_file, bool verbose);

  
  //! Class destructor.
  ~Chromosomes();
  

  //! Returns 'true' if chromosome name is found in the data. 
  bool FindChromosome(char *chromosome_name);
  
  
  //! Returns interval sequence.
  char *GetSeq(GenomicInterval *I, bool replace=false);
  
  
  //! Prints interval sequence.
  void PrintSeq(GenomicInterval *I, bool replace=false);
  
  
  //! Calculates interval sequence length excluding 'N' characters.
  size_t GetSeqLength(GenomicInterval *I);
  
  
  //! Loads chromosome sequence in memory.
  void LoadChromosome(string chromosome_name);

  
  // data

  bool verbose;								//!< if 'true', print messages while loading chromosome sequences
  bool load_in_memory;						//!< 'true' if reading from FASTA file
  StringMap chrom_map;						//!< maps chromosome names to chromosome sequence file names
  string chrom_map_dir;						//!< the directory where the chromosome map file is located
  string chrom_map_name;					//!< the name of the chromosome map file
  string current_chromosome_name;			//!< name of the currently used chromosome
  size_t current_chromosome_size;			//!< size of the currently used chromosome 
  char *current_chromosome_seq;				//!< pointer to the currently used chromosome sequence
  ChromosomeSeqData chromosome_seq_data;	//!< all chromosomes sequences and their sizes are stored here
};


//---------------------------------------------------------------------------------------------//
// END CLASS: Chromosomes                                                                      //
//---------------------------------------------------------------------------------------------//















//---------------------------------------------------------------------------------------------//
// CLASS: GenomicIntervalSetAsArray                                                            //
//---------------------------------------------------------------------------------------------//
//!  [UNDER DEVELOPMENT] This class implements an array of genomic intervals. 
//---------------------------------------------------------------------------------------------//

class GenomicIntervalSetAsArray
{
 public:
  typedef GenomicInterval** iterator;				//!< iterator defined as a pointer in the array of GenomicInterval
  typedef GenomicInterval** reverse_iterator;		//!< reverse_iterator defined as a pointer in the array of GenomicInterval


  //! Empty class constructor.
  GenomicIntervalSetAsArray();


  //! Class destructor.
  ~GenomicIntervalSetAsArray();
  
  
  //! Implementation of begin() for this container.
  iterator begin() { return &I[0]; }
  

  //! Implementation of end() for this container.
  iterator end() { return &I[n_intervals]; }
  

  //! Implementation of rbegin() for this container.
  iterator rbegin() { return &I[n_intervals-1]; }
  

  //! Implementation of rend() for this container.
  iterator rend() { return &I[-1]; }
  

  //! Implementation of front() for this container.
  GenomicInterval *front() { return I[0]; }
  

  //! Implementation of back() for this container.
  GenomicInterval *back() { return I[n_intervals-1]; }
  

  //! Implementation of size() for this container.
  size_t size() { return n_intervals; }
  

  //! Implementation of push_back() for this container.
  void push_back(GenomicInterval *j) { I[n_intervals++] = j; }
  

  //! Implementation of clear() for this container.
  void clear() {  }
  

  //! Implementation of erase() for this container (exits with error message).
  iterator erase(iterator i, iterator j) { fprintf(stderr, "Error: [GenomicIntervalArray]: array entries cannot be erased!"); exit(1); return i; }


  //! Implementation of erase() for this container (exits with error message).
  iterator erase(iterator i) { fprintf(stderr, "Error: [GenomicIntervalArray]: array entries cannot be erased!"); exit(1); return i; }
  

  //! Implementation of insert() for this container (exits with error message).
  iterator insert(iterator i, GenomicInterval *u) { fprintf(stderr, "Error: [GenomicIntervalArray]: array entries cannot be inserted!"); exit(1); return i; }
  

  
  //------------------------------------------------------//
  //   Data                                               //
  //------------------------------------------------------//
  
  long int n_intervals;					//!< number of intervals
  long int max_intervals;				//!< number of maximum allowed intervals (this will be removed)
  GenomicInterval **I;					//!< array of GenomicInterval
  
};





//---------------------------------------------------------------------------------------------//
// END CLASS: GenomicIntervalArray                                                             //
//---------------------------------------------------------------------------------------------//








//---------------------------------------------------------------------------------------------//
// CLASS: GenomicRegion                                                                        //
//---------------------------------------------------------------------------------------------//
//                                                                                             //
// Input formats:                                                                              //
//  1. LABEL <TAB> CHROMOSOME STRAND START[,START]* STOP[,STOP]*                               //
//  2. LABEL <TAB> [CHROMOSOME STRAND START STOP]+                                             //
//                                                                                             //
//---------------------------------------------------------------------------------------------//
//!  This class implements the genomic regions (i.e. ordered sets of genomic intervals).
/*!
  Example:
  \code
    GenomicRegionREG x("Region#1\tchr1 + 132034 135932 chr1 + 135832 140102");		// standard REG format
    GenomicRegionREG y("Region#2\tchr1 + 132034,135832 135932,140102");				// compact REG format
    cout << "x = \n"; x.Print();
    cout << "y = \n"; y.Print(); 
  \endcode
*/
//---------------------------------------------------------------------------------------------//

class GenomicRegion
{
 public:
  //! Empty class constructor.
  GenomicRegion();



  //! Class constructor, see example below.
  GenomicRegion(FileBuffer *B);



  //! Class constructor, see example below.
  GenomicRegion(char *inp, long int n_line=0); 



  //! Class constructor.
  GenomicRegion(char *label, GenomicInterval *i);



  //! Class destructor.
  virtual ~GenomicRegion();
  


  //------------------------------------------------------//
  //   Read & Print                                       //
  //------------------------------------------------------//


  //! Reads label and genomic region intervals from input string <b>inp</b>
  void Read(char *inp, long int n_line=0, char del1=' ', char del2=',');



  //! Prints label and genomic region intervals with a newline character at the end
  virtual void Print(FILE *file_ptr=stdout);


  
  //! Same as Print(), but the label and interval's coordinates are modified.
  virtual void PrintModified(char *label, long int start, long int stop);



  //! Prints all pairwise intersection with the intervals in region <b>r</b> (compatible, sorted and non-overlapping). 
  virtual void PrintIntersection(GenomicRegion *r, bool ignore_strand, bool merge_labels);



  //! Same as Print(), but the interval is constrained between <b>r->I.front()->START</b> and <b>r->I.back()->STOP</b>.
  virtual void PrintConstrained(GenomicRegion *r, bool merge_labels=false);



  //! Prints label and genomic region intervals with a newline character at the end
  void PrintREG(bool compact=false);



  //! Prints sequences (SEQ format)
  void PrintQ(bool reverse=false);



  //! Prints error message and exits
  void PrintError(string error_msg);




  //------------------------------------------------------//
  //   Get & Set                                          //
  //------------------------------------------------------//


  //! Deletes intervals from iterator <b>i</b> to <b>j</b> 
  void DeleteIntervals(GenomicIntervalSet::iterator i, GenomicIntervalSet::iterator j);	



  //! Changes the label
  void SetLabel(char *label);



  //! Interval strands are set to <b>strand</b>. 
  void SetStrand(char strand);



  //! Returns pointer to chromosome name
  char *GetChromosome();



  //! Returns sum of interval's sizes
  size_t GetSize();



  //! Returns genomic region sequence length excluding 'N' characters
  size_t GetSeqLength(Chromosomes *C);



  //! Returns genomic region sequence
  char *GetSeq(Chromosomes *C, bool replace=false);



  //------------------------------------------------------//
  //   Check                                              //
  //------------------------------------------------------//


  //! Returns 'true' if all intervals are in the same chromosome and strand.
  bool IsCompatible(bool ignore_strand);



  //! Returns 'true' if all intervals are compatible with the intervals in <b>r</b>.
  bool IsCompatibleWith(GenomicRegion *r, bool ignore_strand);



  //! Returns 'true' if intervals are compatible and sorted by chromosome, strand (if <b>ignore_strand</b> is false) and start position.
  bool IsCompatibleSorted(bool ignore_strand);



  //! Returns 'true' if intervals are compatible, sorted and non-overlapping. 
  bool IsCompatibleSortedAndNonoverlapping();



  //! Returns 'true' if it overlaps with input region <b>r</b>.
  bool OverlapsWith(GenomicRegion *r, bool ignore_strand);



  //! Compares with input region <b>r</b> and returns true if its order is before <b>r</b>. The order is determined first by chromosome, second by strand, then by start position and finally by stop position.
  bool IsBefore(GenomicRegion *r, bool sorted_by_strand);



  //! Compares with input region <b>r</b> and returns true if its position order is before <b>r</b>. The position order is determined first by start position and then by stop position.
  bool IsPosBefore(GenomicRegion *r);



  //! Returns the total number of overlapping nucleotides (gaps are not matched) with input region <b>r</b>.
  long int CalcOverlap(GenomicRegion *r, bool ignore_strand);



  //! Returns -1 if before interval <b>i</b>, +1 if after, otherwise returns 0. The order is determined first by chromosome, then (optionally) by strand and finally by start position.
  int CalcDirection(GenomicInterval *i, bool sorted_by_strand);

  

  //! Returns -1 if before region <b>r</b>, +1 if after, otherwise returns 0. The order is determined first by chromosome, then (optionally) by strand and finally by start position.
  int CalcDirection(GenomicRegion *r, bool sorted_by_strand);

  

  //------------------------------------------------------//
  //   Line-based (horizontal) operations                 //
  //------------------------------------------------------//


  //! Prints alignments of input sequences to reference genome 
  virtual void RunAlign(Chromosomes *C);


  
  //! Prints in BED format. If <b>convert_chromosome</b> is 'true', then convert from ENSEMBL to UCSC names (by adding 'chr' prefix to each chromosome name, and by converting 'MT' to 'chrM'). 
  virtual void PrintBEDFormat(char *color, bool convert_chromosome);



  //! Corrects interval start/stop positions so as to comply with chromosomal bounds.
  virtual bool ApplyBounds(StringLIntMap *bounds);



  //! Replace each interval in the region with the corresponding center interval
  virtual void Center();



  //! Replace region with a single interval from minimum start to maximum stop position
  virtual void Connect();



  //! Replace regions with the corresponding gap intervals between successive intervals (e.g. used to compute intron boundaries)
  virtual void Diff();



  //! Print distance between successive intervals 
  virtual void RunCalcDistances(char *op1, char *op2);



  //! Divide intervals divided in half (i.e. two intervals for each original interval in the genomic region)
  virtual void Divide();



  //! Print intervals divided in half (i.e. two intervals for each original interval in the genomic region)
  virtual void RunDivide();



  //! Removes rogue intervals (e.g. start>stop)
  virtual void Fix();



  //! Print intersection among all intervals, i.e. maximum start to minimum stop position
  virtual void Intersect();



  //! Prints region label and sum of intervals' sizes 
  void RunSize();



  //! Modifies interval start/stop positions. First it chooses which position to keep fixed, and then it shifts the non-fixed position. 
  /*!
    \param position_op 		selects fixed position as follows: '1'=start position, 'c'=center position, '5p'=5-prime, '3p'=3-prime
    \param position_shift	determines the shift of the non-fixed position (in the 'c' case, both positions are shifted symmetrically around the center)
  */
  virtual void ModifyPos(char *position_op, long int position_shift);           



  //! Randomizes interval position within chromosomal bounds
  virtual void Randomize(gsl_rng *random_generator, StringLIntMap *bounds);



  //! Reverses interval start/stop positions with respect to chromosomal bounds, as if the reference point were in the end of the chromosome as opposed to the beginning. 
  void ReversePos(StringLIntMap *bounds);


  
  //! Selects a subset of intervals according to their relative positions
  virtual void Select(bool first, bool last, bool from5p, bool from3p);



  //! Shifts interval start/stop positions. First it chooses reference position, and then shifts start and stop positions separately
  /*!
    \param start_shift		determines the shift of the start position
    \param stop_shift		determines the shift of the stop position
    \param strand_aware		if 'true', then start=5-prime and stop=3-prime
  */
  virtual void ShiftPos(long int start_shift, long int stop_shift, bool strand_aware);



  //! Prints shuffled region within specified reference regions in <b>refReg</b>.
  virtual void RunShuffle(gsl_rng *random_generator, GenomicRegionSet *refReg, StringLIntMap *index, StringVecLIntMap *loc);



  //! Sort intervals according to start position (only for compatible intervals)
  virtual void Sort(); 


  
  //! Print region intervals on separate lines. 
  virtual void RunSplit(); 


  
  //! Modified strand information: '+'=only forward, '-'=only reverse, 'r'=reverse strands from '+' to '-' and vice versa, 'b'=both strands
  virtual void ModifyStrand(char *strand_op);


  
  //! Compute intervals' union
  virtual void Union();



  //! Print intervals' union
  virtual void RunUnion();



  //! Print sliding windows (implemented only for single-interval regions)
  virtual void PrintWindows(long int win_step, long int win_size);



  //! Remove sub-intervals that correspond to sequences of 'N' characters
  void PrintRemoveN(Chromosomes *C);



  //! Search sequence for short pattern (requires SEQ input)
  virtual void PrintSearch(char *pattern, bool header, bool summary);

  
  //! Prints region label, intervals and extracted sequence in SEQ format
  void PrintSeq(Chromosomes *C, bool replace=false);



  //! Calculates the offset distances of interval start and stop position with respect to the first interval in the set. 
  /*!
    \param op		determines reference point as follows: '1'=start position, '2'=stop position, '5p'=5-prime position, '3p'=3-prime position
    \param fraction 	if 'true', offsets are reported as a fraction of the total region length
  */
  void PrintOffsetFormat(char *op, bool fraction);



  //! Prints region label and sequence length (excluding 'N' characters)
  void PrintSeqLength(Chromosomes *C);



  //! Verifies extracted sequence against region label
  void PrintVerifySeq(Chromosomes *C, bool ignore);



  //! If interval strand is '+', it is changed to '-' and vice versa. 
  void ReverseStrand();



  //------------------------------------------------------//
  //   Operations returning new regions                   //
  //------------------------------------------------------//

  
  //! Returns a version of this genomic region constrained by the start/stop coordinates.
  GenomicRegion *Constrain(GenomicRegion *r, char *label=NULL);

  

  //! Returns the difference of this genomic region and region <b>r</b>.
  virtual GenomicRegion *Diff(GenomicRegion *r);

  

  //------------------------------------------------------//
  //   Data                                               //
  //------------------------------------------------------//

  unsigned long int n_line;				//!< keeps track of line number
  char *LABEL;							//!< genomic region label
  GenomicIntervalSet I;					//!< set of intervals
  
  
};



//---------------------------------------------------------------------------------------------//
// END CLASS: GenomicRegion                                                                    //
//---------------------------------------------------------------------------------------------//













//---------------------------------------------------------------------------------------------//
// CLASS: GenomicRegionSEQ                                                                     //
//---------------------------------------------------------------------------------------------//
//!  This class implements the genomic regions REG format with attached sequence.
//---------------------------------------------------------------------------------------------//

class GenomicRegionSEQ : public GenomicRegion
{
 public:
  //! Class constructor, see example below.
  GenomicRegionSEQ(FileBuffer *B);

  //! Class destructor.
  ~GenomicRegionSEQ();
  


  //------------------------------------------------------//
  //   Operations                                         //
  //------------------------------------------------------//
  
  //! Search sequence for short pattern (requires SEQ input)
  void PrintSearch(char *pattern, bool header, bool summary);

  
  //------------------------------------------------------//
  //   Data                                               //
  //------------------------------------------------------//

  char *SEQ;							//!< sequence

  
/*
  // NOTE: utilize all data & code related to the SEQ format in a new class: GenomicRegionSEQ.
  //! Print sequence with no alignment gaps (requires SEQ input)
  void PrintNoGaps();											// alignments

  //! Print sequence labels in the header (only for SEQ format)
  void PrintModifyLabel();
  
  vector<string> Q;						//!< sequences (e.g. from a multiple sequence alignment)
  vector<string> L;						//!< sequence labels (e.g. names of species in multiple sequence alignment)
*/

};



//---------------------------------------------------------------------------------------------//
// END CLASS: GenomicRegionSEQ                                                                 //
//---------------------------------------------------------------------------------------------//











//---------------------------------------------------------------------------------------------//
// CLASS: GenomicRegionBEDToREG                                                                //
//---------------------------------------------------------------------------------------------//
//!  This class reads BED files but converts them internally to the REG format, so all operations are performed on the REG format.
/*!
  Example:
  \code
    GenomicRegionBEDToREG x("chr1 132033 140102 Region#1 1000 + 132033 140102 0 2 3899,4104 0,3965");
    cout << "x = \n"; x.Print();
  \endcode
*/
//---------------------------------------------------------------------------------------------//
class GenomicRegionBEDToREG : public GenomicRegion
{
 public:
  //! Class constructor, see example below.
  GenomicRegionBEDToREG(FileBuffer *B);



  //! Class constructor, see example below.
  GenomicRegionBEDToREG(char *inp, long int n_line=0);



  //! Class destructor.
  ~GenomicRegionBEDToREG();
  


  //! Reads BED format from input string <b>inp</b>
  void Read(char *inp, long int n_line=0);
  


  // data


};

//---------------------------------------------------------------------------------------------//
// END CLASS: GenomicRegionBEDToREG                                                            //
//---------------------------------------------------------------------------------------------//








//---------------------------------------------------------------------------------------------//
// CLASS: GenomicRegionBED                                                                     //
//---------------------------------------------------------------------------------------------//
//!  This class processes genomic regions in BED format.
/*!
  Example:
  \code
    GenomicRegionBED x("chr1 132033 140102 Region#1 1000 + 132033 140102 0 2 3899,4104 0,3965");
    cout << "x = \n"; x.Print();
  \endcode
*/
//---------------------------------------------------------------------------------------------//
class GenomicRegionBED : public GenomicRegion
{
 public:
  //! Class constructor, see example below.
  GenomicRegionBED(FileBuffer *B);



  //! Class constructor, see example below.
  GenomicRegionBED(char *inp, long int n_line=0);



  //! Class constructor.
  GenomicRegionBED(char *label, GenomicInterval *i);



  //! Class destructor.
  ~GenomicRegionBED();
  


  //------------------------------------------------------//
  //   Read & Print                                       //
  //------------------------------------------------------//


  //! Reads BED format from input string <b>inp</b>
  void Read(char *inp, long int n_line=0);
  


  //! Prints label and genomic region intervals with a newline character at the end
  void Print(FILE *file_ptr=stdout);



  //! Same as generic Print() except that it modifies the color and the chromosome names
  void Print(FILE *file_ptr, char *color, bool convert_chromosome);



  //! Prints all pairwise intersection with the intervals in region <b>r</b> (compatible, sorted and non-overlapping). 
  void PrintIntersection(GenomicRegion *r, bool ignore_strand, bool merge_labels);



  //! Same as Print(), but the interval is constrained between <b>start</b> and <b>stop</b>
  void PrintConstrained(GenomicRegion *r, bool merge_labels=false);



  //! Same as Print(), but the label and interval's coordinates are modified.
  void PrintModified(char *label, long int start, long int stop);



  //------------------------------------------------------//
  //   Line-based (horizontal) operations                 //
  //------------------------------------------------------//


  //! Auxiliary function: updates thickStart/thickEnd variables to honor changes in coordinates
  void UpdateThick();
  
  
  
  //! Prints in BED format. If <b>convert_chromosome</b> is 'true', then convert from ENSEMBL to UCSC names (by adding 'chr' prefix to each chromosome name, and by converting 'MT' to 'chrM'). 
  void PrintBEDFormat(char *color, bool convert_chromosome);



  //! Corrects interval start/stop positions so as to comply with chromosomal bounds.
  bool ApplyBounds(StringLIntMap *bounds);



  //! Replace each interval in the region with the corresponding center interval
  void Center();



  //! Replace region with a single interval from minimum start to maximum stop position
  void Connect();



  //! Replace regions with the corresponding gap intervals between successive intervals (e.g. used to compute intron boundaries)
  void Diff();



  //! Print distance between successive intervals 
  void RunCalcDistances(char *op1, char *op2);



  //! Divide intervals divided in half (i.e. two intervals for each original interval in the genomic region)
  void Divide();



  //! Print intervals divided in half (i.e. two intervals for each original interval in the genomic region)
  void RunDivide();



  //! Removes rogue intervals (e.g. start>stop)
  void Fix();



  //! Print intersection among all intervals, i.e. maximum start to minimum stop position
  void Intersect();



  //! Modifies interval start/stop positions. First it chooses which position to keep fixed, and then it shifts the non-fixed position. 
  /*!
    \param position_op 		selects fixed position as follows: '1'=start position, 'c'=center position, '5p'=5-prime, '3p'=3-prime
    \param position_shift	determines the shift of the non-fixed position (in the 'c' case, both positions are shifted symmetrically around the center)
  */
  void ModifyPos(char *position_op, long int position_shift);           



  //! Randomizes interval position within chromosomal bounds
  void Randomize(gsl_rng *random_generator, StringLIntMap *bounds);



  //! Selects a subset of intervals according to their relative positions
  void Select(bool first, bool last, bool from5p, bool from3p);



  //! Shifts interval start/stop positions. First it chooses reference position, and then shifts start and stop positions separately
  /*!
    \param start_shift		determines the shift of the start position
    \param stop_shift		determines the shift of the stop position
    \param strand_aware		if 'true', then start=5-prime and stop=3-prime
  */
  void ShiftPos(long int start_shift, long int stop_shift, bool strand_aware);



  //! Prints shuffled region within specified reference regions in <b>refReg</b>.
  void RunShuffle(gsl_rng *random_generator, GenomicRegionSet *refReg, StringLIntMap *index, StringVecLIntMap *loc);



  //! Sort intervals according to start position (only for compatible intervals)
  void Sort(); 


  
  //! Print region intervals on separate lines. 
  void RunSplit(); 


  
  //! Compute intervals' union
  void Union();



  //! Print sliding windows (implemented only for single-interval regions)
  void PrintWindows(long int win_step, long int win_size);



  //------------------------------------------------------//
  //   Operations returning new regions                   //
  //------------------------------------------------------//

  
  //! Returns a version of this genomic region constrained by the start/stop coordinates.
  GenomicRegionBED *Constrain(GenomicRegion *r, char *label=NULL);

  

  //! Returns the difference of this genomic region and region <b>r</b>.
  GenomicRegion *Diff(GenomicRegion *r);

  

  //------------------------------------------------------//
  //   Data                                               //
  //------------------------------------------------------//

  long int n_tokens; 			//!< number of tokens
  long int score;				//!< score field
  char *itemRgb;				//!< color field
  long int thickStart; 			//!< thickStart field  
  long int thickEnd; 			//!< thickEnd field  

};

//---------------------------------------------------------------------------------------------//
// END CLASS: GenomicRegionBED                                                                 //
//---------------------------------------------------------------------------------------------//










//---------------------------------------------------------------------------------------------//
// CLASS: GenomicRegionSAMToREG                                                                //
//---------------------------------------------------------------------------------------------//
//!  This class reads SAM files but converts them internally to the REG format, so all operations are performed on the REG format.
//---------------------------------------------------------------------------------------------//
class GenomicRegionSAMToREG : public GenomicRegion
{
 public:
  
  //! Class constructor, see example below.
  GenomicRegionSAMToREG(FileBuffer *B);



  //! Class constructor, see example below.
  GenomicRegionSAMToREG(char *inp, long int n_line=0);



  //! Class destructor.
  ~GenomicRegionSAMToREG();
  


  //! Reads SAM format from input string <b>inp</b>
  void Read(char *inp, long int n_line=0);
  


  //! Extracts the strand from the FLAG field. 
  char CalcStrandFromFlag(unsigned long int flag);


  
  //! Returns next token in CIGAR string.
  long int GetNextTokenOfCIGAR(char **cigar, char *type);

  
  
  // data


};

//---------------------------------------------------------------------------------------------//
// END CLASS: GenomicRegionSAMToREG                                                            //
//---------------------------------------------------------------------------------------------//









//---------------------------------------------------------------------------------------------//
// CLASS: GenomicRegionSAM                                                                     //
//---------------------------------------------------------------------------------------------//
//!  This class processes genomic regions in SAM format.
//---------------------------------------------------------------------------------------------//
class GenomicRegionSAM : public GenomicRegion
{
 public:
  typedef vector<pair<long int,char> > CIGARTokens;			//!< This type is used for storing CIGAR tokens, i.e. (token_type, token_length)
  
  //! Class constructor.
  GenomicRegionSAM(FileBuffer *B);



  //! Class constructor.
  GenomicRegionSAM(char *inp, long int n_line=0);



  //! Class destructor.
  ~GenomicRegionSAM();
  


  //------------------------------------------------------//
  //   Read & Print                                       //
  //------------------------------------------------------//


  //! Reads SAM format from input string <b>inp</b>
  void Read(char *inp, long int n_line=0);



  //! Prints label and genomic region intervals with a newline character at the end
  void Print(FILE *file_ptr=stdout);



  //! Same as Print(), but the interval is constrained between <b>start</b> and <b>stop</b>
  void PrintConstrained(GenomicRegion *r, bool merge_labels=false);



  //! Same as Print(), but the label and interval's coordinates are modified.
  void PrintModified(char *label, long int start, long int stop);


  
  //------------------------------------------------------//
  //   Auxiliary functions                                //
  //------------------------------------------------------//

  //! Extracts the strand from the FLAG field. 
  char CalcStrandFromFlag(unsigned long int flag);


  
  //! Updates the FLAG field given a new strand.
  void UpdateFlagFromStrand(char strand);
  
  
  
  //! Prints gapped sequence.
  void PrintGappedSequence(char *seq, char *gap_token_types, bool lowercase_for_sort_clipping=true);

  
  
  //! Returns gapped sequence.
  char *GetGappedSequence(char *seq, char *gap_token_types, bool lowercase_for_sort_clipping=true);

  
  
  //! Sets CIGAR, SEQ and QUAL fields to "*" and removes OPTIONAL field; also \ref GenomicRegionSAM::n_tokens is set to 11.
  void InvalidateSeqData();


  
  //! Computes and stores CIGAR tokens in \ref GenomicRegionSAM::T.
  void TokenizeCIGAR();

  
  
  //! Returns next token in CIGAR string.
  long int GetNextTokenOfCIGAR(char **cigar, char *type);

  
  
  //! Update CIGAR string based only on interval information (chooses simplest possible CIGAR).
  void UpdateCigarFromIntervals();
  
  
 
  //! Calculate interval length from CIGAR string.
  long int CalcFragmentLengthFromCIGAR();
  
  
  
  //! Calculate interval length from CIGAR string.
  long int CalcReferenceLengthFromCIGAR();
  
  
  
  //! Trim CIGAR string from left and right by a given distance (measured on the reference sequence). 
  char *TrimCIGAR(long int d_ref_start, long int d_ref_stop); 
 
 
 
  //! Trim input <b>seq</b> from left and right (used for trimming both SEQ and QUAL fields). 
  char *TrimSequence(char *seq, long int d_ref_start, long int d_ref_stop); 
 
 
 
  //------------------------------------------------------//
  //   Line-based (horizontal) operations                 //
  //------------------------------------------------------//


  //! Prints alignments of input sequences to reference genome 
  void RunAlign(Chromosomes *C);


  
  //! Corrects interval start/stop positions so as to comply with chromosomal bounds.
  bool ApplyBounds(StringLIntMap *bounds);



  //! Replace each interval in the region with the corresponding center interval
  void Center();



  //! Replace region with a single interval from minimum start to maximum stop position
  void Connect();



  //! Replace regions with the corresponding gap intervals between successive intervals (e.g. used to compute intron boundaries)
  void Diff();



  //! Print distance between successive intervals 
  void RunCalcDistances(char *op1, char *op2);



  //! Divide intervals divided in half (i.e. two intervals for each original interval in the genomic region)
  void Divide();



  //! Print intervals divided in half (i.e. two intervals for each original interval in the genomic region)
  void RunDivide();



  //! Removes rogue intervals (e.g. start>stop)
  void Fix();



  //! Print intersection among all intervals, i.e. maximum start to minimum stop position
  void Intersect();



  //! Randomizes interval position within chromosomal bounds
  void Randomize(gsl_rng *random_generator, StringLIntMap *bounds);



  //! Modifies interval start/stop positions. First it chooses which position to keep fixed, and then it shifts the non-fixed position. 
  /*!
    \param position_op 		selects fixed position as follows: '1'=start position, 'c'=center position, '5p'=5-prime, '3p'=3-prime
    \param position_shift	determines the shift of the non-fixed position (in the 'c' case, both positions are shifted symmetrically around the center)
  */
  void ModifyPos(char *position_op, long int position_shift);           



  //! Selects a subset of intervals according to their relative positions
  void Select(bool first, bool last, bool from5p, bool from3p);



  //! Shifts interval start/stop positions. First it chooses reference position, and then shifts start and stop positions separately
  /*!
    \param start_shift		determines the shift of the start position
    \param stop_shift		determines the shift of the stop position
    \param strand_aware		if 'true', then start=5-prime and stop=3-prime
  */
  void ShiftPos(long int start_shift, long int stop_shift, bool strand_aware);



  //! Prints shuffled region within specified reference regions in <b>refReg</b>.
  void RunShuffle(gsl_rng *random_generator, GenomicRegionSet *refReg, StringLIntMap *index, StringVecLIntMap *loc);



  //! Sort intervals according to start position (only for compatible intervals)
  void Sort(); 


  
  //! Print region intervals on separate lines. 
  void RunSplit(); 


  
  //! Modified strand information: '+'=only forward, '-'=only reverse, 'r'=reverse strands from '+' to '-' and vice versa, 'b'=both strands
  void ModifyStrand(char *strand_op);


  
  //! Compute intervals' union
  void Union();



  //! Print sliding windows (implemented only for single-interval regions)
  void PrintWindows(long int win_step, long int win_size);



  //------------------------------------------------------//
  //   Data                                               //
  //------------------------------------------------------//

  long int n_tokens;					//!< number of TAB-separated fields in SAM entry (>= 12)
  //char *QNAME;						//!< Query template NAME: mapped to GenomicRegion::LABEL
  unsigned long int FLAG;				//!< bitwise FLAG
  //char *RNAME;						//!< Reference sequence NAME: mapped to GenomicRegion::I.front()->CHROMOSOME
  //long int POS;						//!< 1-based leftmost mapping POSition: mapped to GenomicRegion::I.front()->START
  long int MAPQ;						//!< MAPping Quality
  char *CIGAR;							//!< CIGAR string
  char *RNEXT;							//!< Reference name of the mate/next fragment
  long int PNEXT;						//!< Position of the mate/next fragment
  long int TLEN;						//!< observed Template LENgth
  char *SEQ;							//!< fragment SEQuence
  char *QUAL;							//!< ASCII of Phred-scaled base QUALity+33
  char *OPTIONAL;						//!< all SAM fields after the 12th mandatory field
  CIGARTokens T;						//!< CIGAR tokens
  
};

//---------------------------------------------------------------------------------------------//
// END CLASS: GenomicRegionSAM                                                                 //
//---------------------------------------------------------------------------------------------//















//---------------------------------------------------------------------------------------------//
// CLASS: GenomicRegionGFFToREG                                                                //
//---------------------------------------------------------------------------------------------//
//!  This class reads GFF files but converts them internally to the REG format, so all operations are performed on the REG format.
//---------------------------------------------------------------------------------------------//
class GenomicRegionGFFToREG : public GenomicRegion
{
 public:
  //! Class constructor, see example below.
  GenomicRegionGFFToREG(FileBuffer *B);



  //! Class constructor, see example below.
  GenomicRegionGFFToREG(char *inp, long int n_line=0);



  //! Class destructor.
  ~GenomicRegionGFFToREG();
  


  //! Reads GFF format from input string <b>inp</b>
  void Read(char *inp, long int n_line=0);
  


  // data


};

//---------------------------------------------------------------------------------------------//
// END CLASS: GenomicRegionGFFToREG                                                            //
//---------------------------------------------------------------------------------------------//











//---------------------------------------------------------------------------------------------//
// CLASS: GenomicRegionGFF                                                                     //
//---------------------------------------------------------------------------------------------//
//!  This class processes genomic regions in SAM format.
//---------------------------------------------------------------------------------------------//
class GenomicRegionGFF : public GenomicRegion
{
 public:
  //! Class constructor.
  GenomicRegionGFF(FileBuffer *B);



  //! Class constructor.
  GenomicRegionGFF(char *inp, long int n_line=0);



  //! Class destructor.
  ~GenomicRegionGFF();
  


  //------------------------------------------------------//
  //   Read & Print                                       //
  //------------------------------------------------------//


  //! Reads GFF format from input string <b>inp</b>
  void Read(char *inp, long int n_line=0);



  //! Prints label and genomic region intervals with a newline character at the end
  void Print(FILE *file_ptr=stdout);



  //! Same as Print(), but the interval is constrained between <b>start</b> and <b>stop</b>
  void PrintConstrained(GenomicRegion *r, bool merge_labels=false);



  //! Same as Print(), but the label and interval's coordinates are modified.
  void PrintModified(char *label, long int start, long int stop);



  //! Replace regions with the corresponding gap intervals between successive intervals (e.g. used to compute intron boundaries)
  void Diff();



  //! Print distance between successive intervals 
  void RunCalcDistances(char *op1, char *op2);



  //! Divide intervals divided in half (i.e. two intervals for each original interval in the genomic region)
  void Divide();



  //! Print intervals divided in half (i.e. two intervals for each original interval in the genomic region)
  void RunDivide();



  //! Print intersection among all intervals, i.e. maximum start to minimum stop position
  void Intersect();



  //! Selects a subset of intervals according to their relative positions
  void Select(bool first, bool last, bool from5p, bool from3p);



  //! Print region intervals on separate lines. 
  void RunSplit(); 


  
  //! Print sliding windows (implemented only for single-interval regions)
  void PrintWindows(long int win_step, long int win_size);



  //---------------------------------------//
  //  Data                                 //
  //---------------------------------------//
  
  long int n_tokens;			//!< number of TAB-separated fields in GFF entry (= 11)
  //char *SEQNAME; 				//!< sequence name: mapped to GenomicRegion::I->front()->CHROMOSOME
  char *SOURCE;					//!< feature source
  char *FEATURE;				//!< feture type
  //long int START;				//!< 1-based start position: mapped to GenomicRegion::I->front()->START
  //long int END;				//!< 1-based end position: mapped to GenomicRegion::I->front()->STOP
  char *SCORE;					//!< score
  //char STRAND;				//!< strand: mapped to GenomicRegion::I->front()->STRAND
  char FRAME;					//!< coding frame
  //char *ATTRIBUTE;			//!< attribute: mapped to GenomicRegion::LABEL if non-empty
  char *COMMENT;				//!< comment

};

//---------------------------------------------------------------------------------------------//
// END CLASS: GenomicRegionGFF                                                                 //
//---------------------------------------------------------------------------------------------//







//---------------------------------------------------------------------------------------------//
// CLASS: GenomicRegionSet                                                                     //
//---------------------------------------------------------------------------------------------//
//!  This class is used to read and manipulate a set of genomic regions from a file or from standard input.
/*!
  Example:
  \code
    GenomicRegionSet RegSet("test.reg",10000,true,false,false);
    Progress PRG("Printing center of intervals...",1);
    for (GenomicRegion *r=RegSet.Get(); r!=NULL; r=RegSet.Next()) {
      r->Center();
      r->Print();
      PRG.Check();
    }
    PRG.Done();
  \endcode
*/
//---------------------------------------------------------------------------------------------//
class GenomicRegionSet
{
 public:
  //! Class constructor, see example below.
  GenomicRegionSet(char *file, unsigned long int buffer_size, bool verbose, bool load_in_memory, bool hide_header); 


  //! Class constructor, see example below.
  GenomicRegionSet(FILE *file_ptr, unsigned long int buffer_size, bool verbose, bool load_in_memory, bool hide_header);


  //! Class destructor.
  virtual ~GenomicRegionSet();


  //-------------------------------------//
  //  Read & Print                       //
  //-------------------------------------//

  //! Processes file header. 
  void ProcessFileHeader(bool hide);

  
  //! Detects input format.
  void DetectFileFormat();

  
  //! Creates new genomic region based on input format.
  GenomicRegion *CreateGenomicRegion(FileBuffer *buffer);


  //! Initializes data structures. 
  void Init();


  //! Resets file; note, that this is not possible if the regions are being read from the standard input. 
  void Reset();


  //! Returns pointer to the first region in the set; note, that this is not possible if the regions are being read from the standard input.
  GenomicRegion *Begin(); 


  //! Returns pointer to current region in the set.
  virtual GenomicRegion *Get();


  //! Moves pointer to next region in the set; if \ref GenomicRegionSet::load_in_memory is 'false', it deletes the current region.
  virtual GenomicRegion *Next();


  //! Same as \ref Next() but allows user to check sort order and/or retain the current region into memory. 
  /*!
    \param sorted_by_strand 	if 'true', it takes into account strand information while checking the sort order
    \param retain_current	if 'true', the current region is not deleted from memory; this must be done 'manually' using the \ref Release() method (do not use 'delete'). 
  */
  GenomicRegion *Next(bool sorted_by_strand, bool retain_current=false);


  //! Deletes genomic region from memory. 
  void Release(GenomicRegion *r);


  //! [DEPRECATED] Returns current region, and protects its from deletion by \ref GenomicRegionSet::Next().
  GenomicRegion *GetRetain();


  //! Prints error message
  void PrintError(string error_msg);
  


  //-------------------------------------//
  //  File-based (vertical) operations   //
  //-------------------------------------//


  //! Counts the number of regions in the set. If <b>use_label_counts</b>, the region labels are assumed to be numbers, and the sum of these numbers is reported. 
  long int CountRegions(bool use_label_counts);


  //! Annotate input region neighborhoods as upstream, downstream, etc. 
  void RunGlobalAnnotate(StringLIntMap *bounds);


  //! Cluster input regions based on pair-wise overlaps (note that this is computationally intensive).
  void RunGlobalCluster(bool merge=false);


  //! Print distance between successive regions 
  void RunGlobalCalcDistances(char *op1, char *op2);


  //! Prints sorted region set.
  void RunGlobalSort();


  //! Prints the difference between chromosome intervals and input region set (input region set must be sorted).
  void RunGlobalInvert(StringLIntMap *bounds);


  //! Links successive regions if they overlap (input region set must be sorted).
  void RunGlobalLink(bool sorted_by_strand, long int max_difference);


  //! Partitions overlapping regions into a non-overlapping region set (input region set must be sorted).
  void RunGlobalPartition();


  //! Tests if input region set is sorted.
  void RunGlobalTest(bool sorted_by_strand);


  //! Merges successive regions if they are separated by a distance of <b>merge_win</b> or less (input region set must be sorted).
  void RunGlobalMerge(long int merge_win);


  //! Reverses the order of reverse-strand regions (for this operation, the entire region set must be loaded in memory).
  void RunGlobalReverseOrder();


  //! Computes densities of input regions in index regions in <b>Index</b>. This is OBSOLETE, use class SortedGenomicRegionSetOverlaps instead. 
  void PrintDensities(GenomicRegionSet &Index, bool report_exact_boundaries);


  //! Scans input regions in sliding windows and reports count distribution (input regions must be sorted).
  void RunGlobalScan(StringLIntMap *bounds, long int win_step, long int win_size);


  //! Scans input regions in sliding windows (input regions must be sorted). See also class GenomicRegionSetScanner.  
  /*!
    \param bounds 					chromosome sizes
    \param ref_reg_file 			only windows that overlap with regions in this file will be reported
    \param win_step					sliding window step
    \param win_size					sliding window size (must be a multiple of window step)
    \param ignore_reverse_strand	if true, no sliding windows on the negative strand are reported
    \param preprocess				if '1', only start position is counted; if 'p', all positions are counted; if 'c', center of interval is counted
    \param use_labels_as_values		if true, genomic region labels are assumed to be integers and are included in the counting
    \param min_reads				report windows only if value is greater that this parameter
  */
  void RunGlobalScanCount(StringLIntMap *bounds, char *ref_reg_file, long int win_step, long int win_size, bool ignore_reverse_strand, char preprocess, bool use_labels_as_values, long int min_reads);


  //! Prints in Wiggle format. If <b>convert_chromosome</b> is 'true', then convert from ENSEMBL to UCSC names (by adding 'chr' prefix to each chromosome name, and by converting 'MT' to 'chrM'). 
  void RunConvertToWIG(char *title, char *color, char *position, char *options, long int span, bool convert_chromosome);




  //---------------------------------------//
  //  Line-based (horizontal) operations   //
  //---------------------------------------//

  //! Prints alignments of input sequences to reference genome 
  void RunAlign(Chromosomes *C);


  //! Converts to BED format. If <b>convert_chromosome</b> is 'true', then convert from ENSEMBL to UCSC names (by adding 'chr' prefix to each chromosome name, and by converting 'MT' to 'chrM'). 
  void RunConvertToBED(char *title, char *color, char *position, bool convert_chromosome);


  //! Corrects interval start/stop positions so as to comply with chromosomal bounds 
  void RunBounds(StringLIntMap *bounds);


  //! Print intervals' center
  void RunCenter();


  //! Print connected intervals as a single interval from minimum start to maximum stop position
  void RunConnect();


  //! Print gap intervals between successive intervals (e.g. for computing intron boundaries)
  void RunDiff();


  //! Print distances between successive intervals of each region (i.e. local operation)
  void RunCalcDistances(char *op1, char *op2);


  //! Print intervals divided in half (i.e. two intervals for each original interval in the genomic region)
  void RunDivide();


  //! Removes rogue intervals (e.g. start>stop)
  void RunFix();



  //! Print intersection among all intervals, i.e. maximum start to minimum stop position
  void RunIntersection();



  //! Prints region label and sum of intervals' sizes 
  void RunSize();



  //! Modifies interval start/stop positions. First it chooses which position to keep fixed, and then it shifts the non-fixed position. 
  /*!
    \param position_op 		selects fixed position as follows: '1'=start position, 'c'=center position, '5p'=5-prime, '3p'=3-prime
    \param position_shift	determines the shift of the non-fixed position (in the 'c' case, both positions are shifted symmetrically around the center)
  */
  void RunModifyPos(char *position_op, long int position_shift);           



  //! Prints label and genomic region intervals with a newline character at the end for all regions in the set
  void RunConvertToREG(bool compact=false);



  //! Randomizes interval position within chromosomal bounds
  void RunRandomize(gsl_rng *random_generator, StringLIntMap *bounds);



  //! Selects a subset of intervals according to their relative positions
  void RunSelect(bool first, bool last, bool from5p, bool from3p);



  //! Shifts interval start/stop positions. First it chooses reference position, and then shifts start and stop positions separately
  /*!
    \param start_shift		determines the shift of the start position
    \param stop_shift		determines the shift of the stop position
    \param strand_aware		if 'true', then start=5-prime and stop=3-prime
  */
  void RunShiftPos(long int start_shift, long int stop_shift, bool strand_aware);


  //! Print sorted intervals according to start position (only for compatible intervals)
  void RunSort(); 

  
  //! Split regions into intervals which are printed on separate lines. 
  void RunSplit(); 

  
  //! Print region with modified strands: '+'=only forward, '-'=only reverse, 'r'=reverse strands from '+' to '-' and vice versa; if <b>sorted</b> is 'true', it calls \ref GenomicRegionSet::RunModifyStrandSorted instead.
  void RunModifyStrand(char *strand_op, bool sorted);


  //! Same as \ref GenomicRegionSet::RunModifyStrand but returns a sorted region set (note: it only works if the original region set is sorted).
  void RunModifyStrandSorted(char *strand_op);


  //! Shuffles input regions within specified reference regions loaded from file <b>ref_reg_file</b>
  void RunShuffle(gsl_rng *random_generator, char *ref_reg_file);


  //! Print intervals' union
  void RunUnion();


  //! Print sliding windows (implemented only for single-interval regions)
  void RunSlidingWindows(long int win_step, long int win_size);


  //! Prints region label, intervals and extracted sequence in SEQ format
  void RunExtractSeq(Chromosomes *C, bool replace=false);



  //---------------------------------------//
  //  UNDER DEVELOPMENT/TESTING            //
  //---------------------------------------//


  //! Prints in BEDGraph format. If <b>convert_chromosome</b> is 'true', then convert from ENSEMBL to UCSC names (by adding 'chr' prefix to each chromosome name, and by converting 'MT' to 'chrM'). 
  void PrintBEDGraphFormat(char *title, char *color, char *position, bool convert_chromosome);


  //! Prints region label and sequence length excluding 'N' characters
  void PrintSeqLength(Chromosomes *C);


  //! Print sequence with no alignment gaps (requires SEQ input)
  //void PrintNoGaps();


  //! Remove sub-intervals that correspond to sequences of 'N' characters
  void PrintRemoveN(Chromosomes *C);


  //! Calculates the offset distances of interval start and stop position with respect to the first interval in the set. 
  /*!
    \param op		determines reference point as follows: '1'=start position, '2'=stop position, '5p'=5-prime position, '3p'=3-prime position
    \param fraction 	if 'true', offsets are reported as a fraction of the total region length
  */
  void PrintOffsetFormat(char *op, bool fraction);


  //! Reverses interval start/stop positions with respect to chromosomal bounds, as if the reference point were in the end of the chromosome as opposed to the beginning. 
  void PrintReversePos(StringLIntMap *bounds);


  //! Search sequence for short pattern (requires SEQ input)
  void PrintSearch(char *pattern, bool header, bool summary);


  //! Verifies extracted sequence against region label
  void PrintVerifySeq(Chromosomes *C, bool ignore);



 private:
  //! Creates a temporary file from this instance of GenomicRegionSet which contains the current chromosome's positive strand regions (used in \ref GenomicRegionSet::RunModifyStrandSorted)
  GenomicRegionSet *StoreInTempFile();


  //! Uses mergesort between <b>tempSet</b> and this instance of GenomicRegionSet (used in \ref GenomicRegionSet::RunModifyStrandSorted)
  GenomicRegion *PrintMergeSort(GenomicRegionSet *rtempSet, char strand);


  //! Reverses the order between <b>tempSet</b> and this instance of GenomicRegionSet (used in \ref GenomicRegionSet::RunModifyStrandSorted)
  GenomicRegion *PrintReverse(GenomicRegionSet *rtempSet);



  //---------------------------------------//
  //  Data                                 //
  //---------------------------------------//

 public:
  char *file;						//!< file name ('NULL' if reading from standard input)
  FILE *file_ptr;					//!< file pointer ('NULL' if reading from standard input)
  unsigned long int buffer_size;	//!< line buffer size (automatically adjusted during execution to ensure the entire line is read)
  bool verbose;						//!< verbose mode
  bool from_stdin;					//!< 'true' if reading from standard input
  bool load_in_memory;				//!< 'true' if the entire region set is loaded in memory
  bool hide_header;					//!< if 'true', then do not print file header
  bool use_interval_array; 			//!< if 'true', use an array for storing genomic regions' intervals (for time efficiency)
  string format;					//!< input format
  Progress progress;	 			//!< keeps track of progress of computation
  FileBuffer *buffer;				//!< pointer to the buffer class used to read the region sets
  long int n_regions;				//!< number of regions ('1' if the regions are not loaded in memory)
  GenomicRegion **R;				//!< array where the regions are stored; if the region set is not entirely loaded in memory, only <b>R[0]</b> is used
  long int r_index;					//!< index pointing to the current region in <b>R</b>
};



//---------------------------------------------------------------------------------------------//
// END CLASS: GenomicRegionSet                                                                 //
//---------------------------------------------------------------------------------------------//
















//---------------------------------------------------------------------------------------------//
// CLASS: GenomicRegionSetScanner                                                              //
//---------------------------------------------------------------------------------------------//
//!  This class is used to scan a set of regions by sliding windows. See detailed description below for an example.  
/*!
     Example:
  \code
    // initialize (note: you need to set inputs and parameters, such as genome_reg_file, input_reg_file, WIN_DIST and WIN_SIZE)
    StringLIntMap *bounds = ReadBounds(genome_reg_file);
    GenomicRegionSet InputRegSet(input_reg_file,10000,true,false);
    GenomicRegionSetScanner input_scanner(&InputRegSet,bounds,WIN_DIST,WIN_SIZE,false,false,'c');

    // run
    Progress PRG("Scanning...",1);
    for (long int v=input_scanner.Next(); v!=-1; v=input_scanner.Next()) {
      if (v>=MIN_READS) {
        cout << v << '\t';
        input_scanner.PrintInterval();
        cout << '\n';
      }
      PRG.Check();
    }
    PRG.Done();
  
    // cleanup
    if (bounds!=NULL) delete bounds;
  \endcode
*/
//---------------------------------------------------------------------------------------------//
class GenomicRegionSetScanner
{
 public:
  //! Class constructor.
  /*!
    \param R 				the GenomicRegionSet which will be scanned by sliding windows
    \param bounds 			chromosome sizes
    \param win_step			sliding window step
    \param win_size			sliding window size (must be a multiple of window step)
    \param use_labels_as_values		if true, genomic region labels are assumed to be integers and are included in the counting
    \param ignore_reverse_strand	if true, no sliding windows on the negative strand are reported
    \param preprocess			if '1', only start position is counted; if 'p', all positions are counted; if 'c', center of interval is counted
  */
  GenomicRegionSetScanner(GenomicRegionSet *R, StringLIntMap *bounds, long int win_step, long int win_size, bool use_labels_as_values, bool ignore_reverse_strand, char preprocess);
  ~GenomicRegionSetScanner();
  
  // operations
  long int Next();				//!< computes value in the next window
  long int Next(GenomicRegionSet *Ref);		//!< computes value in the next window that overlaps with <b>Ref</b>
  void PrintInterval();				//!< prints current window's interval
  
  // data
 private: 
  bool Test();					//!< tests whether current input region should be skipped because it does not match any chromosome name in the provided bounds

  GenomicRegionSet *R;				//!< pointer to the GenomicRegionSet which will be scanned by sliding windows
  StringLIntMap *bounds;			//!< chromosome sizes
  long int win_step;				//!< sliding window step
  long int win_size;				//!< sliding window size (must be a multiple of window step)
  bool use_labels_as_values;			//!< if true, genomic region labels are assumed to be integers and are included in the counting
  bool ignore_reverse_strand;			//!< if true, no sliding windows on the negative strand are reported
  char preprocess;				//!< if '1', only start position is counted; if 'p', all positions are counted; if 'c', center of interval is counted
  long int n_win_combine;			//!< win_size divided by win_step (note: each window is comprised of non-overlapping sub-windows of size win_step)
  GenomicRegion *r;				//!< pointer to the most recently scanned region 
  StringLIntMap::iterator chr;			//!< keeps track of current window's chromosome information
  char strand;					//!< keeps track of current window's strand information
  long int start;				//!< keeps track of current window's start information
  long int stop;				//!< keeps track of current window's stop information
  long int *v;					//!< a vector that keeps track of values in the sub-windows (of size win_step) that comprise the window (of size win_size)
  long int k;					//!< the k-th element of v to be updated next (in a round-robin manner)
  long int v_sum;				//!< keeps track of current window's value information (this is the sum of elements in vector v); it is set to '-1' if no more windows are available.
};

//---------------------------------------------------------------------------------------------//
// END CLASS: GenomicRegionSetScanner                                                          //
//---------------------------------------------------------------------------------------------//











//---------------------------------------------------------------------------------------------//
// CLASS: GenomicRegionSetOverlaps                                                             //
//---------------------------------------------------------------------------------------------//
//!  Abstract class for computing overlaps between genomic regions.
//---------------------------------------------------------------------------------------------//
class GenomicRegionSetOverlaps
{
 public:
  //! Class constructor.
  /*!
    \param QuerySet 				pointer to query region set
    \param IndexSet 				pointer to index region set
  */
  GenomicRegionSetOverlaps(GenomicRegionSet *QuerySet, GenomicRegionSet *IndexSet);


  //! Class destructor.
  virtual ~GenomicRegionSetOverlaps();

  
  //! Returns a pointer to the current query region ('NULL' if no more query regions are available).
  virtual GenomicRegion *GetQuery() = 0;


  //! Returns a pointer to the next query region ('NULL' if no more query regions are available).
  virtual GenomicRegion *NextQuery() = 0;


  //! Returns a pointer to the current index region that overlaps the current query region (even if the overlap is only in the gaps between intervals).
  virtual GenomicRegion *GetMatch() = 0;

 
  //! Returns a pointer to the next index region that overlaps the current query region (even if the overlap is only in the gaps between intervals).
  virtual GenomicRegion *NextMatch() = 0;

 
  //! Returns true, if no more query or index regions are available and the \ref SortedGenomicRegionSetOverlaps::IRegBuffer buffer is empty.
  virtual bool Done() = 0;
 

  //! Returns a pointer to the current index region that overlaps the current query region.
  /*!
    \param match_gaps 			if 'true', overlaps are defined as in \ref GetMatch.
    \param ignore_strand		if 'true', overlaps are strand-ignorant
  */
  GenomicRegion *GetOverlap(bool match_gaps, bool ignore_strand);


  //! Returns a pointer to the next index region that overlaps the current query region (even if the overlap is only in the gaps between intervals).
  /*!
    \param match_gaps 			if 'true', overlaps are defined as in \ref GetMatch.
    \param ignore_strand		if 'true', overlaps are strand-ignorant
  */
  GenomicRegion *NextOverlap(bool match_gaps, bool ignore_strand);

  
  //! Calculates the total overlap between the current query region and the index region set.
  /*!
    \param match_gaps 			if 'true', overlaps are defined as in \ref GetMatch.
    \param ignore_strand		if 'true', overlaps are strand-ignorant
    \param use_labels_as_values 	region labels contain number to be used in calculation
  */
  unsigned long int CalcQueryCoverage(bool match_gaps, bool ignore_strand, bool use_labels_as_values);


  //! Calculates the total overlap for each index region. The index region set must be loaded in memory.  
  /*!
    \param match_gaps 			if 'true', overlaps are defined as in \ref GetMatch.
    \param ignore_strand		if 'true', overlaps are strand-ignorant
    \param use_labels_as_values 	region labels contain number to be used in calculation
  */
  unsigned long int *CalcIndexCoverage(bool match_gaps, bool ignore_strand, bool use_labels_as_values);


  //! Calculates the total number of matches between the current query region and the index region set.
  /*!
    \param match_gaps 			if 'true', overlaps are defined as in \ref GetMatch.
    \param ignore_strand		if 'true', overlaps are strand-ignorant
    \param use_labels_as_values 	region labels contain number to be used in calculation
   */
  unsigned long int CountQueryOverlaps(bool match_gaps, bool ignore_strand, bool use_labels_as_values);


  //! Calculates the total number of matches for each index region. The index region set must be loaded in memory. 
  /*!
    \param match_gaps 			if 'true', overlaps are defined as in \ref GetMatch.
    \param ignore_strand		if 'true', overlaps are strand-ignorant
    \param use_labels_as_values 	region labels contain number to be used in calculation
   */
  unsigned long int *CountIndexOverlaps(bool match_gaps, bool ignore_strand, bool use_labels_as_values);

 
  
  // data
 public: 
  GenomicRegionSet *QuerySet;				//!< pointer to query region set
  GenomicRegionSet *IndexSet;				//!< pointer to index region set
  GenomicRegion *current_qreg;				//!< pointer to the current query region 
  GenomicRegion *current_ireg;				//!< pointer to the current index region 
   
};

//---------------------------------------------------------------------------------------------//
// END CLASS: GenomicRegionSetOverlaps                                                         //
//---------------------------------------------------------------------------------------------//











//---------------------------------------------------------------------------------------------//
// CLASS: UnsortedGenomicRegionSetOverlaps                                                     //
//---------------------------------------------------------------------------------------------//
//! Class for computing overlaps between unsorted regions.
/*!
    A simple example for computing RNAseq read densities in known exons (this is actually implemented in the <b>genomic_overlaps</b> command-line tool as 'density' operation):
  \code
    #include "core.h"
    #include "genomic_intervals.h"
    
    // open region sets
    char *REF_REG_FILE = argv[next_arg];
    char *TEST_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet RefRegSet(REF_REG_FILE,BUFFER_SIZE,VERBOSE,true);
    GenomicRegionSet TestRegSet(TEST_REG_FILE,BUFFER_SIZE,VERBOSE,false);

    // process overlaps
    GenomicRegionSetOverlaps *overlaps = new UnsortedGenomicRegionSetOverlaps(&TestRegSet,&RefRegSet);
    unsigned long int *coverage = overlaps->CalcIndexCoverage(MATCH_GAPS,IGNORE_STRAND,USE_VALUES); 
    Progress PRG("Printing densities...",RefRegSet.n_regions);
    for (long int k=0; k<RefRegSet.n_regions; k++) {
      GenomicRegion *qreg = RefRegSet.R[k];
      long int qreg_size = MATCH_GAPS ? (qreg->I.back()->STOP-qreg->I.front()->START+1) : qreg->GetSize();
      double density = (double)coverage[k]/qreg_size;
      if (density>=MIN_DENSITY) printf("%s\t%.4e\n", qreg->LABEL, density);
      PRG.Check();
    }
    PRG.Done();
    delete coverage;
    delete overlaps;
  \endcode
*/
//---------------------------------------------------------------------------------------------//
class UnsortedGenomicRegionSetOverlaps : public GenomicRegionSetOverlaps
{
 public:
  //! Class constructor.
  /*!
    \param QuerySet 			pointer to query region set
    \param IndexSet 			pointer to index region set
    \param bin_bits			number of shift-bits per bin level, e.g. "10,15,18"
  */
  UnsortedGenomicRegionSetOverlaps(GenomicRegionSet *QuerySet, GenomicRegionSet *IndexSet, char *bin_bits=NULL);


  //! Class destructor.
  ~UnsortedGenomicRegionSetOverlaps();

  
  //! Returns a pointer to the current query region ('NULL' if no more query regions are available).
  GenomicRegion *GetQuery();


  //! Returns a pointer to the next query region ('NULL' if no more query regions are available).
  GenomicRegion *NextQuery();


  //! Returns a pointer to the current index region that overlaps the current query region (even if the overlap is only in the gaps between intervals).
  GenomicRegion *GetMatch();

 
  //! Returns a pointer to the next index region that overlaps the current query region (even if the overlap is only in the gaps between intervals).
  GenomicRegion *NextMatch();

 
  //! Returns true, if no more query or index regions are available and the \ref SortedGenomicRegionSetOverlaps::IRegBuffer buffer is empty.
  bool Done();
 

  
  // data
 public: 
  typedef pair<long int*,long int**> BinSet;		//!< a two-dimensional grid of bins
  long int *r_next;					//!< for each region store, store a pointer to the "next" one in the same bin
  int n_levels;						//!< number of levels in the two-dimensional grid of bins
  int *n_bits;						//!< number of bits (per level) by which the start position is right-shifted to compute its bin
  map<string,BinSet*> index;				//!< the two-dimensional grid of bins per chromosome
  BinSet *current_binset;				//!< pointer to the current chromosome set of bins
  bool new_query; 					//!< 'true', if a new query is about to be processed
};

//---------------------------------------------------------------------------------------------//
// END CLASS: UnsortedGenomicRegionSetOverlaps                                                 //
//---------------------------------------------------------------------------------------------//











//---------------------------------------------------------------------------------------------//
// CLASS: SortedGenomicRegionSetOverlaps                                                       //
//---------------------------------------------------------------------------------------------//
//!  This class is used to find overlaps between two genomic regions sets. See detailed description below for an example.  
/*!
     A simple example for computing RNAseq read densities in known exons (this is actually implemented in the <b>genomic_overlaps</b> command-line tool as 'density' operation):
  \code
    #include "core.h"
    #include "genomic_intervals.h"
 
    // open region sets
    char *REF_REG_FILE = argv[next_arg];
    char *TEST_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet RefRegSet(REF_REG_FILE,BUFFER_SIZE,VERBOSE,true);
    GenomicRegionSet TestRegSet(TEST_REG_FILE,BUFFER_SIZE,VERBOSE,false);

    // process overlaps
    GenomicRegionSetOverlaps *overlaps = new SortedGenomicRegionSetOverlaps(&TestRegSet,&RefRegSet,SORTED_BY_STRAND);
    unsigned long int *coverage = overlaps->CalcIndexCoverage(MATCH_GAPS,IGNORE_STRAND,USE_VALUES); 
    Progress PRG("Printing densities...",RefRegSet.n_regions);
    for (long int k=0; k<RefRegSet.n_regions; k++) {
      GenomicRegion *qreg = RefRegSet.R[k];
      long int qreg_size = MATCH_GAPS ? (qreg->I.back()->STOP-qreg->I.front()->START+1) : qreg->GetSize();
      double density = (double)coverage[k]/qreg_size;
      if (density>=MIN_DENSITY) printf("%s\t%.4e\n", qreg->LABEL, density);
      PRG.Check();
    }
    PRG.Done();
    delete coverage;
    delete overlaps;
  \endcode
  A more complex example for creating ChIPseq read profiles in TSS regions (this is actually implemented in the <b>genomic_overlaps</b> command-line tool as 'offset' operation):
  \code
    #include "core.h"
    #include "genomic_intervals.h"
 
    // open region sets
    char *QUERY_REG_FILE = "chipseq.reads.reg";
    char *INDEX_REG_FILE = "TSS.flank10kb.reg";
    char *OFFSET_OP = "5p";
    bool SORTED_BY_STRAND = false; 
    bool MATCH_GAPS = true;
    bool IGNORE_STRAND = true;
    GenomicRegionSet QueryRegSet(QUERY_REG_FILE,10000,true,false);
    GenomicRegionSet IndexRegSet(INDEX_REG_FILE,10000,true,false);

    // process overlaps
    SortedGenomicRegionSetOverlaps Overlaps(&QueryRegSet,&IndexRegSet,SORTED_BY_STRAND);
    Progress PRG("Processing queries...",1);
    for (GenomicRegion *qreg=Overlaps.GetQuery(); Overlaps.Done()==false; qreg=Overlaps.NextQuery()) {
      size_t Qsize = qreg->GetSize();
      for (GenomicRegion *ireg=Overlaps.GetOverlap(MATCH_GAPS,IGNORE_STRAND); ireg!=NULL; ireg=Overlaps.NextOverlap(MATCH_GAPS,IGNORE_STRAND)) {
        cout << qreg->LABEL << '\t';
        long int start_offset, stop_offset;
        ireg->I.front()->GetOffsetFrom(qreg->I.front(),OFFSET_OP,IGNORE_STRAND,&start_offset,&stop_offset);
        printf("%f", ((float)start_offset/Qsize+(float)stop_offset/Qsize)/2); 
        cout << '\n';
      }
      PRG.Check();
    }
    PRG.Done();
  \endcode
*/
//---------------------------------------------------------------------------------------------//
class SortedGenomicRegionSetOverlaps : public GenomicRegionSetOverlaps
{
 public:
  //! Class constructor.
  /*!
    \param QuerySet 				pointer to query region set
    \param IndexSet 				pointer to index region set
    \param sorted_by_strand			if true, query and index regions are sorted by strand
  */
  SortedGenomicRegionSetOverlaps(GenomicRegionSet *QuerySet, GenomicRegionSet *IndexSet, bool sorted_by_strand);


  //! Class destructor.
  ~SortedGenomicRegionSetOverlaps();

  
  //! Returns a pointer to the current query region ('NULL' if no more query regions are available); also, \ref SortedGenomicRegionSetOverlaps::IRegBuffer is automatically updated via \ref SortedGenomicRegionSetOverlaps::LoadIndexBuffer.
  GenomicRegion *GetQuery();


  //! Returns a pointer to the next query region ('NULL' if no more query regions are available); also, \ref SortedGenomicRegionSetOverlaps::IRegBuffer is automatically updated via \ref SortedGenomicRegionSetOverlaps::LoadIndexBuffer.
  GenomicRegion *NextQuery();


  //! Returns a pointer to the current index region that overlaps the current query region (even if the overlap is only in the gaps between intervals).
  GenomicRegion *GetMatch();

 
  //! Returns a pointer to the next index region that overlaps the current query region (even if the overlap is only in the gaps between intervals).
  GenomicRegion *NextMatch();

 
  //! Returns true, if no more query or index regions are available and the \ref SortedGenomicRegionSetOverlaps::IRegBuffer buffer is empty.
  bool Done();
 

 private:
  void LoadIndexBuffer();				//!< manages \ref SortedGenomicRegionSetOverlaps::IRegBuffer
  void ClearIndexBuffer();				//!< clears \ref SortedGenomicRegionSetOverlaps::IRegBuffer

  // data
  bool sorted_by_strand;				//!< if true, query and index regions sorted by strand
  GenomicRegionList IRegBuffer;				//!< temporary buffer containing index regions that overlap with current query; the buffer is cleared only if the next query has no overlap 
  GenomicRegionList::iterator IRegBufferIterator;	//!< iterator on \ref SortedGenomicRegionSetOverlaps::IRegBuffer
  GenomicInterval *ireg_buffer_interval;		//!< keeps tracks of the (linked) interval stored in \ref SortedGenomicRegionSetOverlaps::IRegBuffer
  unsigned long int max_ireg_buffer_size;		//!< keeps track of the maximum size of \ref SortedGenomicRegionSetOverlaps::IRegBuffer 
};

//---------------------------------------------------------------------------------------------//
// END CLASS: SortedGenomicRegionSetOverlaps                                                   //
//---------------------------------------------------------------------------------------------//










char ProcessStrand(char *strand);
void PrintChromosome(char *chromosome, bool convert);
GenomicRegion *RegCenter(GenomicRegion *r);

StringLIntMap *ReadBounds(char *genome_reg_file, bool verbose=false);		//!< reads chromosome sizes from REG file <b>genome_reg_file</b>
unsigned long int CalcBoundSize(StringLIntMap *bounds);
unsigned long int CalcRegSize(char *reg_file);
bool CompareGenomicIntervals(GenomicInterval *I, GenomicInterval *J);
void SortGenomicIntervals(GenomicIntervalSetAsList *L);
void SortGenomicIntervals(GenomicIntervalSetAsVector *V);
void SortGenomicIntervals(GenomicIntervalSetAsArray *A);





