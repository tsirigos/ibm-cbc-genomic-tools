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
#include <map>
#include <string>
#include <list>
#include <vector>
#include <iostream>
#include <algorithm>
#include "core.h"
#include "genomic_intervals.h"



using namespace std;


//---------------------------------------------------------------------------------//
// Constants                                                                       //
//---------------------------------------------------------------------------------//

const string PROGRAM = "genomic_overlaps";
const long int BUFFER_SIZE = 10000;



  
//---------------------------------------------------------------------------------//
// Command-line options                                                            //
//---------------------------------------------------------------------------------//

bool HELP;
bool VERBOSE;
char *BIN_BITS;		// NOTE: for future use in case we want to control bin shift-bits
bool IS_SORTED;
bool SORTED_BY_STRAND;
bool IGNORE_STRAND;
bool MERGE_LABELS;
bool PRINT_LABELS;
bool MATCH_GAPS;
bool USE_VALUES;
unsigned long int MIN_COUNT;
double MIN_DENSITY;
char *OFFSET_OP;
bool OFFSET_FRACTION;
bool CENTER;
bool SUBSET_NONOVERLAPS;




//-------InitCmdLine-----------
//
CmdLineWithOperations *InitCmdLine(int argc, char *argv[], int *next_arg)
{
  // initialize
  CmdLineWithOperations *cmd_line = new CmdLineWithOperations(); 
  cmd_line->SetProgramName(PROGRAM,VERSION);

  // set operations
  cmd_line->AddOperation("count", "[OPTIONS] REFERENCE-REGION-FILE <TEST-REGION-FILE>", \
  "Counts the number of overlapping test regions per reference region.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operands: region, region-set\n\
  * Region requirements: chromosome/strand-compatible, sorted, non-overlapping\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("coverage", "[OPTIONS] REFERENCE-REGION-FILE <TEST-REGION-FILE>", \
  "Calculates the depth coverage (i.e. the total number of overlapping nucleotides) per reference region.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operands: region, region-set\n\
  * Region requirements: chromosome/strand-compatible, sorted, non-overlapping\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("density", "[OPTIONS] REFERENCE-REGION-FILE <TEST-REGION-FILE>", \
  "Computes the density (i.e. the coverage divided by the size of the reference region) of overlaps per reference region.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operands: region, region-set\n\
  * Region requirements: chromosome/strand-compatible, sorted, non-overlapping\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("intersect", "[OPTIONS] REFERENCE-REGION-FILE <TEST-REGION-FILE>", \
  "Computes the intersection between all pairs of test and reference regions. Results are grouped by test region.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operands: region, region-set\n\
  * Test region requirements: chromosome/strand-compatible, sorted, non-overlapping\n\
  * Reference region requirements: single-interval regions\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("offset", "[OPTIONS] REFERENCE-REGION-FILE <TEST-REGION-FILE>", \
  "Computes the distances of test regions from their overlapping reference regions.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operands: region, region-set\n\
  * Region requirements: chromosome/strand-compatible, sorted, non-overlapping\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("overlap", "[OPTIONS] REFERENCE-REGION-FILE <TEST-REGION-FILE>", \
  "Finds the overlaps between all pairs of test and reference regions. Results are grouped by test region.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operands: region, region-set\n\
  * Region requirements: chromosome/strand-compatible, sorted, non-overlapping\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("subset", "[OPTIONS] REFERENCE-REGION-FILE <TEST-REGION-FILE>", \
  "Picks a subset of test regions depending on their overlap with reference regions. Results are grouped by test region.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operands: region, region-set\n\
  * Region requirements: chromosome/strand-compatible, sorted, non-overlapping\n\
  * Region-set requirements: none"\
  );

  if (argc<2) { cmd_line->OperationSummary("OPERATION [OPTIONS] REFERENCE-REGION-FILE <TEST-REGION-FILE>","Performs overlap operations between a test and a reference set of genomic regions."); exit(1); }

  // current operation
  string op = argv[1];
  if (op[0]=='-') op = op.substr(1);		// ensure compatibility with previous version
  cmd_line->SetCurrentOperation(op);
  
  // common options
  cmd_line->AddOption("--help", &HELP, false, "help");
  cmd_line->AddOption("-h", &HELP, false, "help");
  cmd_line->AddOption("-v", &VERBOSE, false, "verbose mode");
  cmd_line->AddOption("-B", &BIN_BITS, "17,20,23,26", "number of shift-bits for each bin level");
  cmd_line->AddOption("-S", &IS_SORTED, false, "test and reference regions are sorted by chromosome and start position");
  cmd_line->AddOption("-s", &SORTED_BY_STRAND, false, "test and reference regions are also sorted by strand (-S must be set)");
  cmd_line->AddOption("-i", &IGNORE_STRAND, false, "ignore strand while finding overlaps");
   
  // Main options
  if (op=="count") {
    cmd_line->AddOption("-gaps", &MATCH_GAPS, false, "matching gaps between intervals are considered overlaps");
    cmd_line->AddOption("-val", &USE_VALUES, false, "use values contained in the labels of index intervals");
    cmd_line->AddOption("-min", &MIN_COUNT, 0, "minimum count");
  }
  else if (op=="coverage") {
    cmd_line->AddOption("-gaps", &MATCH_GAPS, false, "matching gaps between intervals are considered overlaps");
    cmd_line->AddOption("-val", &USE_VALUES, false, "use values contained in the labels of index intervals");
    cmd_line->AddOption("-min", &MIN_COUNT, 0, "minimum coverage");
  }
  else if (op=="density") {
    cmd_line->AddOption("-gaps", &MATCH_GAPS, false, "matching gaps between intervals are considered overlaps");
    cmd_line->AddOption("-val", &USE_VALUES, false, "use values contained in the labels of index intervals");
    cmd_line->AddOption("-min", &MIN_DENSITY, 0.0, "minimum density");
  }
  else if (op=="intersect") {
    cmd_line->AddOption("-label", &MERGE_LABELS, false, "print query label for each match");
  }
  else if (op=="offset") {
    cmd_line->AddOption("-gaps", &MATCH_GAPS, false, "matching gaps between intervals are considered overlaps");
    cmd_line->AddOption("-label", &PRINT_LABELS, false, "print test region labels");
    cmd_line->AddOption("-op", &OFFSET_OP, "1", "reference point (1=start, 2=stop, 5p=5'-end, 3p=3'-end)");
    cmd_line->AddOption("-a", &OFFSET_FRACTION, false, "print distances as a fraction of total size");
    cmd_line->AddOption("-c", &CENTER, false, "print center of interval only");
  }
  else if (op=="overlap") {
    cmd_line->AddOption("-gaps", &MATCH_GAPS, false, "matching gaps between intervals are considered overlaps");
    cmd_line->AddOption("-label", &MERGE_LABELS, false, "print query label for each match");
  }
  else if (op=="subset") {
    cmd_line->AddOption("-gaps", &MATCH_GAPS, false, "matching gaps between intervals are considered overlaps");
    cmd_line->AddOption("-inv", &SUBSET_NONOVERLAPS, false, "print test regions that do *not* overlap with reference regions");
  }
  else {
    cerr << "Unknown operation '" << op << "'!\n";
    delete cmd_line;
    exit(1);
  }

  // process command line
  *next_arg = cmd_line->Read(argv+1,argc-1) + 1;
  if (HELP||(argc-*next_arg<1)) { cmd_line->OperationUsage(); delete cmd_line; exit(1); }
  
  return cmd_line;
}






//---------------------------------------------------------------------------------//
// MAIN                                                                            //
//---------------------------------------------------------------------------------//
int main(int argc, char* argv[]) 
{
  // process command-line arguments
  int next_arg;
  CmdLineWithOperations *cmd_line = InitCmdLine(argc,argv,&next_arg); 
  _MESSAGES_ = VERBOSE;

  if (IS_SORTED&&SORTED_BY_STRAND&&IGNORE_STRAND) { fprintf(stderr, "[Error]: the input is sorted by chromosome/strand/start (i.e. -S and -s are set), therefore the overlap algorithm can only report strand-specific results (i.e. -i cannot be set)!\n"); exit(1); }

  //--------------------
  // count/sorted
  //--------------------
  /*if ((cmd_line->current_cmd_operation=="count")&&(IS_SORTED==true)) {
    // open region sets
    char *REF_REG_FILE = argv[next_arg];
    char *TEST_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet RefRegSet(REF_REG_FILE,BUFFER_SIZE,VERBOSE,false,true);
    GenomicRegionSet TestRegSet(TEST_REG_FILE,BUFFER_SIZE,VERBOSE,false,true);

    SortedGenomicRegionSetOverlaps Overlaps(&RefRegSet,&TestRegSet,SORTED_BY_STRAND);
    Progress PRG("Processing queries...",1);
    for (GenomicRegion *qreg=Overlaps.GetQuery(); qreg!=NULL; qreg=Overlaps.NextQuery()) {
      unsigned long int hits = Overlaps.CountQueryOverlaps(MATCH_GAPS,IGNORE_STRAND,USE_VALUES); 
      if (hits>=MIN_COUNT) {
        printf("%s", qreg->LABEL); printf("\t");
        printf("%ld", hits); printf("\n");
      }
      PRG.Check();
    }
    PRG.Done();
  }*/


  //--------------------
  // count
  //--------------------
  if (cmd_line->current_cmd_operation=="count") {
    // open region sets
    char *REF_REG_FILE = argv[next_arg];
    char *TEST_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet RefRegSet(REF_REG_FILE,BUFFER_SIZE,VERBOSE,true,true);
    GenomicRegionSet TestRegSet(TEST_REG_FILE,BUFFER_SIZE,VERBOSE,false,true);

    GenomicRegionSetOverlaps *overlaps;
    if (IS_SORTED) overlaps = new SortedGenomicRegionSetOverlaps(&TestRegSet,&RefRegSet,SORTED_BY_STRAND);
    else overlaps = new UnsortedGenomicRegionSetOverlaps(&TestRegSet,&RefRegSet,BIN_BITS);
    unsigned long int *hits = overlaps->CountIndexOverlaps(MATCH_GAPS,IGNORE_STRAND,USE_VALUES); 
    Progress PRG("Printing counts...",RefRegSet.n_regions);
    for (long int k=0; k<RefRegSet.n_regions; k++) {
      if (hits[k]>=MIN_COUNT) {
        printf("%s", RefRegSet.R[k]->LABEL); printf("\t");
        printf("%lu", hits[k]); printf("\n");
      }
      PRG.Check();
    }
    PRG.Done();
    delete hits;
    delete overlaps;
  }


  //--------------------
  // coverage/sorted
  //--------------------
/*
  else if ((cmd_line->current_cmd_operation=="coverage")&&(IS_SORTED==true)) {
    // open region sets
    char *REF_REG_FILE = argv[next_arg];
    char *TEST_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet RefRegSet(REF_REG_FILE,BUFFER_SIZE,VERBOSE,false,true);
    GenomicRegionSet TestRegSet(TEST_REG_FILE,BUFFER_SIZE,VERBOSE,false,true);

    SortedGenomicRegionSetOverlaps Overlaps(&RefRegSet,&TestRegSet,SORTED_BY_STRAND);
    Progress PRG("Processing queries...",1);
    for (GenomicRegion *qreg=Overlaps.GetQuery(); qreg!=NULL; qreg=Overlaps.NextQuery()) {
      unsigned long int coverage = Overlaps.CalcQueryCoverage(MATCH_GAPS,IGNORE_STRAND,USE_VALUES); 
      if (coverage>=MIN_COUNT) {
        printf("%s", qreg->LABEL); printf("\t");
        printf("%ld", coverage); printf("\n");
      }
      PRG.Check();
    }
    PRG.Done();
  }
*/


  //--------------------
  // coverage
  //--------------------
  else if (cmd_line->current_cmd_operation=="coverage") {
    // open region sets
    char *REF_REG_FILE = argv[next_arg];
    char *TEST_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet RefRegSet(REF_REG_FILE,BUFFER_SIZE,VERBOSE,true,true);
    GenomicRegionSet TestRegSet(TEST_REG_FILE,BUFFER_SIZE,VERBOSE,false,true);

    GenomicRegionSetOverlaps *overlaps;
    if (IS_SORTED) overlaps = new SortedGenomicRegionSetOverlaps(&TestRegSet,&RefRegSet,SORTED_BY_STRAND);
    else overlaps = new UnsortedGenomicRegionSetOverlaps(&TestRegSet,&RefRegSet,BIN_BITS);
    unsigned long int *coverage = overlaps->CalcIndexCoverage(MATCH_GAPS,IGNORE_STRAND,USE_VALUES); 
    Progress PRG("Printing coverages...",RefRegSet.n_regions);
    for (long int k=0; k<RefRegSet.n_regions; k++) {
      if (coverage[k]>=MIN_COUNT) {
        printf("%s", RefRegSet.R[k]->LABEL); printf("\t");
        printf("%lu", coverage[k]); printf("\n");
      }
      PRG.Check();
    }
    PRG.Done();
    delete coverage;
    delete overlaps;
  }


/*
  //--------------------
  // density/sorted
  //--------------------
  else if ((cmd_line->current_cmd_operation=="density")&&(IS_SORTED==true)) {    
    // open region sets
    char *REF_REG_FILE = argv[next_arg];
    char *TEST_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet RefRegSet(REF_REG_FILE,BUFFER_SIZE,VERBOSE,false,true);
    GenomicRegionSet TestRegSet(TEST_REG_FILE,BUFFER_SIZE,VERBOSE,false,true);

    SortedGenomicRegionSetOverlaps Overlaps(&RefRegSet,&TestRegSet,SORTED_BY_STRAND);
    Progress PRG("Processing queries...",1);
    for (GenomicRegion *qreg=Overlaps.GetQuery(); qreg!=NULL; qreg=Overlaps.NextQuery()) {
      long int qreg_size = MATCH_GAPS ? (qreg->I.back()->STOP-qreg->I.front()->START+1) : qreg->GetSize();
      double density = (double)Overlaps.CalcQueryCoverage(MATCH_GAPS,IGNORE_STRAND,USE_VALUES)/qreg_size; 
      if (density>=MIN_DENSITY) printf("%s\t%.4e\n", qreg->LABEL, density);
      PRG.Check();
    }
    PRG.Done();
  }
*/


  //--------------------
  // density
  //--------------------
  else if (cmd_line->current_cmd_operation=="density") {
    // open region sets
    char *REF_REG_FILE = argv[next_arg];
    char *TEST_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet RefRegSet(REF_REG_FILE,BUFFER_SIZE,VERBOSE,true,true);
    GenomicRegionSet TestRegSet(TEST_REG_FILE,BUFFER_SIZE,VERBOSE,false,true);

    // process overlaps
    GenomicRegionSetOverlaps *overlaps;
    if (IS_SORTED) overlaps = new SortedGenomicRegionSetOverlaps(&TestRegSet,&RefRegSet,SORTED_BY_STRAND);
    else overlaps = new UnsortedGenomicRegionSetOverlaps(&TestRegSet,&RefRegSet,BIN_BITS);
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
  }


  //--------------------
  // offset/sorted
  //--------------------
  else if ((cmd_line->current_cmd_operation=="offset")&&(IS_SORTED==true)) {
    // open region sets
    char *REF_REG_FILE = argv[next_arg];
    char *TEST_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet RefRegSet(REF_REG_FILE,BUFFER_SIZE,VERBOSE,false,true);
    GenomicRegionSet TestRegSet(TEST_REG_FILE,BUFFER_SIZE,VERBOSE,false,true);

    // compute offsets
    SortedGenomicRegionSetOverlaps Overlaps(&RefRegSet,&TestRegSet,SORTED_BY_STRAND);
    Progress PRG("Processing queries...",1);
    for (GenomicRegion *qreg=Overlaps.GetQuery(); Overlaps.Done()==false; qreg=Overlaps.NextQuery()) {
      size_t Qsize = qreg->GetSize();
      for (GenomicRegion *ireg=Overlaps.GetOverlap(MATCH_GAPS,IGNORE_STRAND); ireg!=NULL; ireg=Overlaps.NextOverlap(MATCH_GAPS,IGNORE_STRAND)) {
        cout << qreg->LABEL << '\t';
        if (PRINT_LABELS) cout << ireg->LABEL << ' ';
        long int start_offset, stop_offset;
        ireg->I.front()->GetOffsetFrom(qreg->I.front(),OFFSET_OP,IGNORE_STRAND,&start_offset,&stop_offset);
        if (start_offset>stop_offset) { fprintf(stderr, "Error: start offset is greater than stop offest (this must be a bug)!\n"); exit(1); }
        if (CENTER) {
          if (OFFSET_FRACTION) printf("%f", ((float)start_offset/Qsize+(float)stop_offset/Qsize)/2);
          else printf("%ld", (start_offset+stop_offset)/2);
        }
        else {
          if (OFFSET_FRACTION) printf("%f %f", (float)start_offset/Qsize, (float)stop_offset/Qsize);
          else printf("%ld %ld", start_offset, stop_offset);
        }
        cout << '\n';
      }
      PRG.Check();
    }
    PRG.Done();
  }


  //--------------------
  // offset/unsorted
  //--------------------
  else if ((cmd_line->current_cmd_operation=="offset")&&(IS_SORTED==false)) {
    // open region sets
    char *REF_REG_FILE = argv[next_arg];
    char *TEST_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet RefRegSet(REF_REG_FILE,BUFFER_SIZE,VERBOSE,true,true);
    GenomicRegionSet TestRegSet(TEST_REG_FILE,BUFFER_SIZE,VERBOSE,false,true);

    // compute offsets
    UnsortedGenomicRegionSetOverlaps Overlaps(&TestRegSet,&RefRegSet,BIN_BITS);
    Progress PRG("Processing queries...",1);
    for (GenomicRegion *qreg=Overlaps.GetQuery(); Overlaps.Done()==false; qreg=Overlaps.NextQuery()) {
      for (GenomicRegion *ireg=Overlaps.GetOverlap(MATCH_GAPS,IGNORE_STRAND); ireg!=NULL; ireg=Overlaps.NextOverlap(MATCH_GAPS,IGNORE_STRAND)) {
        size_t Isize = ireg->GetSize();
        printf("%s\t", ireg->LABEL);
        if (PRINT_LABELS) printf("%s ", qreg->LABEL);
        long int start_offset, stop_offset;
        qreg->I.front()->GetOffsetFrom(ireg->I.front(),OFFSET_OP,IGNORE_STRAND,&start_offset,&stop_offset);
        if (start_offset>stop_offset) { fprintf(stderr, "Error: start offset is greater than stop offest (this must be a bug)!\n"); exit(1); }
        if (CENTER) {
          if (OFFSET_FRACTION) printf("%f", ((float)start_offset/Isize+(float)stop_offset/Isize)/2);
          else printf("%ld", (start_offset+stop_offset)/2);
        }
        else {
          if (OFFSET_FRACTION) printf("%f %f", (float)start_offset/Isize, (float)stop_offset/Isize);
          else printf("%ld %ld", start_offset, stop_offset);
        }
        printf("\n");
      }
      PRG.Check();
    }
    PRG.Done();
  }


  //--------------------
  // intersect
  //--------------------
  else if (cmd_line->current_cmd_operation=="intersect") {
    // open region sets
    char *REF_REG_FILE = argv[next_arg];
    char *TEST_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet RefRegSet(REF_REG_FILE,BUFFER_SIZE,VERBOSE,IS_SORTED?false:true,true);
    GenomicRegionSet TestRegSet(TEST_REG_FILE,BUFFER_SIZE,VERBOSE,false,false);

    // process
    GenomicRegionSetOverlaps *overlaps;
    if (IS_SORTED) overlaps = new SortedGenomicRegionSetOverlaps(&TestRegSet,&RefRegSet,SORTED_BY_STRAND);
    else overlaps = new UnsortedGenomicRegionSetOverlaps(&TestRegSet,&RefRegSet,BIN_BITS);
    Progress PRG("Processing queries...",1);
    for (GenomicRegion *qreg=overlaps->GetQuery(); overlaps->Done()==false; qreg=overlaps->NextQuery()) {   
      GenomicRegion *ireg = overlaps->GetOverlap(false,IGNORE_STRAND);
      if (ireg!=NULL) {
        do { 
          qreg->PrintConstrained(ireg,MERGE_LABELS);
          ireg = overlaps->NextOverlap(false,IGNORE_STRAND);
        } while (ireg!=NULL);
      }
      PRG.Check();
    }
    PRG.Done();
    delete overlaps;
  }


  //--------------------
  // overlap
  //--------------------
  else if (cmd_line->current_cmd_operation=="overlap") {
    // open region sets
    char *REF_REG_FILE = argv[next_arg];
    char *TEST_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet RefRegSet(REF_REG_FILE,BUFFER_SIZE,VERBOSE,IS_SORTED?false:true,true);
    GenomicRegionSet TestRegSet(TEST_REG_FILE,BUFFER_SIZE,VERBOSE,false,false);

    // process
    GenomicRegionSetOverlaps *overlaps;
    if (IS_SORTED) overlaps = new SortedGenomicRegionSetOverlaps(&TestRegSet,&RefRegSet,SORTED_BY_STRAND);
    else overlaps = new UnsortedGenomicRegionSetOverlaps(&TestRegSet,&RefRegSet,BIN_BITS);
    Progress PRG("Computing query overlaps...",1);
    for (GenomicRegion *qreg=overlaps->GetQuery(); overlaps->Done()==false; qreg=overlaps->NextQuery()) {
      GenomicRegion *ireg = overlaps->GetOverlap(MATCH_GAPS,IGNORE_STRAND);
      if (ireg!=NULL) {
        do { 
          if (MERGE_LABELS) {
            char *old_label = qreg->LABEL;
            char *new_label = new char[strlen(qreg->LABEL)+1+strlen(ireg->LABEL)+1];
            sprintf(new_label, "%s:%s", qreg->LABEL, ireg->LABEL);
            qreg->LABEL = new_label;
            qreg->Print();
            qreg->LABEL = old_label;
            delete new_label;
          }
          else qreg->Print();
          ireg = overlaps->NextOverlap(MATCH_GAPS,IGNORE_STRAND);
        } while (ireg!=NULL);
      }
      PRG.Check();
    }
    PRG.Done();
    delete overlaps;
  }

  

  //--------------------
  // subset
  //--------------------
  else if (cmd_line->current_cmd_operation=="subset") {
    // open region sets
    char *REF_REG_FILE = argv[next_arg];
    char *TEST_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet RefRegSet(REF_REG_FILE,BUFFER_SIZE,VERBOSE,IS_SORTED?false:true,true);
    GenomicRegionSet TestRegSet(TEST_REG_FILE,BUFFER_SIZE,VERBOSE,false,false);

    // process
    GenomicRegionSetOverlaps *overlaps;
    if (IS_SORTED) overlaps = new SortedGenomicRegionSetOverlaps(&TestRegSet,&RefRegSet,SORTED_BY_STRAND);
    else overlaps = new UnsortedGenomicRegionSetOverlaps(&TestRegSet,&RefRegSet,BIN_BITS);
    Progress PRG("Creating query subset...",1);
    for (GenomicRegion *qreg=overlaps->GetQuery(); (SUBSET_NONOVERLAPS&&(qreg!=NULL))||(overlaps->Done()==false); qreg=overlaps->NextQuery()) {   
      if ((overlaps->GetOverlap(MATCH_GAPS,IGNORE_STRAND)==NULL)==SUBSET_NONOVERLAPS) qreg->Print(); 
      PRG.Check();
    }
    PRG.Done();
    delete overlaps;
  }



  // clean up
  delete cmd_line;
  
  return 0;
}



