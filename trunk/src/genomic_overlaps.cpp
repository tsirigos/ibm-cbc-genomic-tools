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
const string VERSION = "2.0.0beta";
const long int BUFFER_SIZE = 10000;



  
//---------------------------------------------------------------------------------//
// Command-line options                                                            //
//---------------------------------------------------------------------------------//

bool VERBOSE;
bool HELP;
char *OFFSET_OP;
bool OFFSET_FRACTION;
bool USE_VALUES;
bool CENTER;
unsigned long int MIN_COUNT;
double MIN_DENSITY;
bool SORTED_BY_STRAND;
bool IGNORE_STRAND;
bool PRINT_QUERY_LABEL;
bool IS_UNSORTED;
bool SUBSET_NONOVERLAPS;
bool MATCH_GAPS;




//-------InitCmdLine-----------
//
CmdLineWithOperations *InitCmdLine(int argc, char *argv[], int *next_arg)
{
  // initialize
  string program = PROGRAM + " (version " + VERSION + ")";
  CmdLineWithOperations *cmd_line = new CmdLineWithOperations(); 
  cmd_line->SetProgramName(program);

  // set operations
  cmd_line->AddOperation("count",     "count [OPTIONS] QUERY-REG-FILE <INDEX-REG-FILE>",     "Counts number of matching index regions per query region.");
  cmd_line->AddOperation("coverage",  "coverage [OPTIONS] QUERY-REG-FILE <INDEX-REG-FILE>",  "Calculates (breadth) coverage per query region.");
  cmd_line->AddOperation("density",   "density [OPTIONS] QUERY-REG-FILE <INDEX-REG-FILE>",   "Computes density of matches per query region.");
  cmd_line->AddOperation("intersect", "intersect [OPTIONS] INDEX-REG-FILE <QUERY-REG-FILE>", "Computes intersection between all pairs of query and index regions.");
  cmd_line->AddOperation("offset",    "offset [OPTIONS] QUERY-REG-FILE <INDEX-REG-FILE>",    "Computes distances of index regions from their overlapping query regions.");
  cmd_line->AddOperation("overlap",   "overlap [OPTIONS] INDEX-REG-FILE <QUERY-REG-FILE>",   "Finds overlaps between all pairs of query and index regions.");
  cmd_line->AddOperation("subset",    "subset [OPTIONS] INDEX-REG-FILE <QUERY-REG-FILE>",    "Picks a subset of query regions depending on their overlap with index regions.");
  if (argc<2) { cmd_line->OperationSummary("OPERATION [OPTIONS] INPUT-FILES","Performs overlap operations between a query and an index set of genomic regions. Results are always grouped by query."); exit(1); }

  // current operation
  string op = argv[1];
  if (op[0]=='-') op = op.substr(1);		// ensure compatibility with previous version
  cmd_line->SetCurrentOperation(op);
  
  // common options
  cmd_line->AddOption("-v", &VERBOSE, false, "verbose mode");
  cmd_line->AddOption("--help", &HELP, false, "help");
  cmd_line->AddOption("-h", &HELP, false, "help");
  cmd_line->AddOption("-s", &SORTED_BY_STRAND, false, "true, if input regions are sorted by strand");
  cmd_line->AddOption("-i", &IGNORE_STRAND, false, "ignore strand while finding overlaps");
   
  // Main options
  if (op=="count") {
    cmd_line->AddOption("-val", &USE_VALUES, false, "use values contained in the labels of index intervals");
    cmd_line->AddOption("-min", &MIN_COUNT, 0, "minimum count");
  }
  else if (op=="coverage") {
    cmd_line->AddOption("-val", &USE_VALUES, false, "use values contained in the labels of index intervals");
    cmd_line->AddOption("-min", &MIN_COUNT, 0, "minimum coverage");
  }
  else if (op=="density") {
    cmd_line->AddOption("-val", &USE_VALUES, false, "use values contained in the labels of index intervals");
    cmd_line->AddOption("-min", &MIN_DENSITY, 0.0, "minimum density");
  }
  else if (op=="intersect") {
    cmd_line->AddOption("-label", &PRINT_QUERY_LABEL, false, "print query label for each match");
  }
  else if (op=="offset") {
    cmd_line->AddOption("-gaps", &MATCH_GAPS, false, "matching gaps between intervals are considered overlaps");
    cmd_line->AddOption("-val", &USE_VALUES, false, "use values contained in the labels of index intervals");
    cmd_line->AddOption("-op", &OFFSET_OP, "1", "reference point (1=start, 2=stop, 5p=5'-end, 3p=3'-end)");
    cmd_line->AddOption("-a", &OFFSET_FRACTION, false, "print distances as a fraction of total size");
    cmd_line->AddOption("-c", &CENTER, false, "print center of interval only");
  }
  else if (op=="overlap") {
    //cmd_line->AddOption("-U", &IS_UNSORTED, false, "use algorithm for unsorted query region set");
    cmd_line->AddOption("-gaps", &MATCH_GAPS, false, "matching gaps between intervals are considered overlaps");
    cmd_line->AddOption("-label", &PRINT_QUERY_LABEL, false, "print query label for each match");
  }
  else if (op=="subset") {
    cmd_line->AddOption("-gaps", &MATCH_GAPS, false, "matching gaps between intervals are considered overlaps");
    cmd_line->AddOption("-no", &SUBSET_NONOVERLAPS, false, "print query regions that do *not* overlap with index regions");
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


  //--------------------
  // count
  //--------------------
  if (cmd_line->current_cmd_operation=="count") {
    // open region sets
    char *QUERY_REG_FILE = argv[next_arg];
    char *INDEX_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet QueryRegSet(QUERY_REG_FILE,BUFFER_SIZE,VERBOSE,false);
    GenomicRegionSet IndexRegSet(INDEX_REG_FILE,BUFFER_SIZE,VERBOSE,false);

    GenomicRegionSetOverlapScanner Overlaps(&QueryRegSet,&IndexRegSet,SORTED_BY_STRAND);
    Progress PRG("Processing queries...",1);
    for (GenomicRegion *qreg=Overlaps.GetQuery(); qreg!=NULL; qreg=Overlaps.NextQuery()) {
      unsigned long int hits = Overlaps.CountMatches(USE_VALUES,IGNORE_STRAND); 
      if (hits>=MIN_COUNT) {
        printf("%s", qreg->LABEL); printf("\t");
        printf("%ld", hits); printf("\n");
      }
      PRG.Check();
    }
    PRG.Done();
  }


  //--------------------
  // coverage
  //--------------------
  else if (cmd_line->current_cmd_operation=="coverage") {
    // open region sets
    char *QUERY_REG_FILE = argv[next_arg];
    char *INDEX_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet QueryRegSet(QUERY_REG_FILE,BUFFER_SIZE,VERBOSE,false);
    GenomicRegionSet IndexRegSet(INDEX_REG_FILE,BUFFER_SIZE,VERBOSE,false);

    GenomicRegionSetOverlapScanner Overlaps(&QueryRegSet,&IndexRegSet,SORTED_BY_STRAND);
    Progress PRG("Processing queries...",1);
    for (GenomicRegion *qreg=Overlaps.GetQuery(); qreg!=NULL; qreg=Overlaps.NextQuery()) {
      unsigned long int hits = Overlaps.CalcCoverage(USE_VALUES,IGNORE_STRAND); 
      if (hits>=MIN_COUNT) {
        printf("%s", qreg->LABEL); printf("\t");
        printf("%ld", hits); printf("\n");
      }
      PRG.Check();
    }
    PRG.Done();
  }


  //--------------------
  // density
  //--------------------
  else if (cmd_line->current_cmd_operation=="density") {    
    // open region sets
    char *QUERY_REG_FILE = argv[next_arg];
    char *INDEX_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet QueryRegSet(QUERY_REG_FILE,BUFFER_SIZE,VERBOSE,false);
    GenomicRegionSet IndexRegSet(INDEX_REG_FILE,BUFFER_SIZE,VERBOSE,false);

    GenomicRegionSetOverlapScanner Overlaps(&QueryRegSet,&IndexRegSet,SORTED_BY_STRAND);
    Progress PRG("Processing queries...",1);
    for (GenomicRegion *qreg=Overlaps.GetQuery(); qreg!=NULL; qreg=Overlaps.NextQuery()) {
      double density = (double)Overlaps.CalcCoverage(USE_VALUES,IGNORE_STRAND)/qreg->GetSize(); 
      if (density>=MIN_DENSITY) printf("%s\t%.4e\n", qreg->LABEL, density);
      PRG.Check();
    }
    PRG.Done();
  }


  //--------------------
  // offset
  //--------------------
  else if (cmd_line->current_cmd_operation=="offset") {
    // open region sets
    char *QUERY_REG_FILE = argv[next_arg];
    char *INDEX_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet QueryRegSet(QUERY_REG_FILE,BUFFER_SIZE,VERBOSE,false);
    GenomicRegionSet IndexRegSet(INDEX_REG_FILE,BUFFER_SIZE,VERBOSE,false);

    // compute offsets
    GenomicRegionSetOverlapScanner Overlaps(&QueryRegSet,&IndexRegSet,SORTED_BY_STRAND);
    Progress PRG("Processing queries...",1);
    for (GenomicRegion *qreg=Overlaps.GetQuery(); Overlaps.Done()==false; qreg=Overlaps.NextQuery()) {
      size_t Qsize = qreg->GetSize();
      for (GenomicRegion *ireg=Overlaps.GetOverlap(MATCH_GAPS,IGNORE_STRAND); ireg!=NULL; ireg=Overlaps.NextOverlap(MATCH_GAPS,IGNORE_STRAND)) {
        cout << qreg->LABEL << '\t';
        if (USE_VALUES) cout << ireg->LABEL << ' ';
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
  // intersect
  //--------------------
  else if (cmd_line->current_cmd_operation=="intersect") {
    // open region sets
    char *INDEX_REG_FILE = argv[next_arg];
    char *QUERY_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet IndexRegSet(INDEX_REG_FILE,BUFFER_SIZE,VERBOSE,false);
    GenomicRegionSet QueryRegSet(QUERY_REG_FILE,BUFFER_SIZE,VERBOSE,false);

    // process
    GenomicRegionSetOverlapScanner Overlaps(&QueryRegSet,&IndexRegSet,SORTED_BY_STRAND);
    Progress PRG("Processing queries...",1);
    for (GenomicRegion *qreg=Overlaps.GetQuery(); Overlaps.Done()==false; qreg=Overlaps.NextQuery()) {   
      GenomicRegion *ireg = Overlaps.GetOverlap(false,IGNORE_STRAND);
      if (ireg!=NULL) {
        do { 
          qreg->PrintConstrained(ireg,PRINT_QUERY_LABEL);
          ireg = Overlaps.NextOverlap(false,IGNORE_STRAND);
        } while (ireg!=NULL);
      }
      PRG.Check();
    }
    PRG.Done();
  }


  //--------------------
  // overlap/unsorted
  //--------------------
  /*else if ((cmd_line->current_cmd_operation=="overlap")&&(IS_UNSORTED)) {
    // open region sets
    char *INDEX_REG_FILE = argv[next_arg];
    char *QUERY_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    IndexedGenomicRegionSet IndexRegSet(INDEX_REG_FILE,BUFFER_SIZE,VERBOSE);
    GenomicRegionSet QueryRegSet(QUERY_REG_FILE,BUFFER_SIZE,VERBOSE,false);
    Progress PRG("Searching for overlaps...",1);
    //for (GenomicRegion *r=IndexRegSet.Get(); r!=NULL; r=IndexRegSet.Next(),PRG.Check()) QueryRegSet.Find(r);
    PRG.Done();
  }*/


  //--------------------
  // overlap
  //--------------------
  else if (cmd_line->current_cmd_operation=="overlap") {
    // open region sets
    char *INDEX_REG_FILE = argv[next_arg];
    char *QUERY_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet IndexRegSet(INDEX_REG_FILE,BUFFER_SIZE,VERBOSE,false);
    GenomicRegionSet QueryRegSet(QUERY_REG_FILE,BUFFER_SIZE,VERBOSE,false);

    // process
    GenomicRegionSetOverlapScanner Overlaps(&QueryRegSet,&IndexRegSet,SORTED_BY_STRAND);
    Progress PRG("Computing query overlaps...",1);
    for (GenomicRegion *qreg=Overlaps.GetQuery(); Overlaps.Done()==false; qreg=Overlaps.NextQuery()) {
      GenomicRegion *ireg = Overlaps.GetOverlap(MATCH_GAPS,IGNORE_STRAND);
      if (ireg!=NULL) {
        do { 
          if (PRINT_QUERY_LABEL) {
            char *old_label = qreg->LABEL;
            char *new_label = new char[strlen(qreg->LABEL)+1+strlen(ireg->LABEL)+1];
            sprintf(new_label, "%s:%s", qreg->LABEL, ireg->LABEL);
            qreg->LABEL = new_label;
            qreg->Print();
            qreg->LABEL = old_label;
            delete new_label;
          }
          else qreg->Print();
          ireg = Overlaps.NextOverlap(MATCH_GAPS,IGNORE_STRAND);
        } while (ireg!=NULL);
      }
      PRG.Check();
    }
    PRG.Done();
  }

  

  //--------------------
  // subset
  //--------------------
  else if (cmd_line->current_cmd_operation=="subset") {
    // open region sets
    char *INDEX_REG_FILE = argv[next_arg];
    char *QUERY_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet IndexRegSet(INDEX_REG_FILE,BUFFER_SIZE,VERBOSE,false);
    GenomicRegionSet QueryRegSet(QUERY_REG_FILE,BUFFER_SIZE,VERBOSE,false);

    // process
    GenomicRegionSetOverlapScanner Overlaps(&QueryRegSet,&IndexRegSet,SORTED_BY_STRAND);	
    Progress PRG("Creating query subset...",1);
    for (GenomicRegion *qreg=Overlaps.GetQuery(); Overlaps.Done()==false; qreg=Overlaps.NextQuery()) {   
      if ((Overlaps.GetOverlap(MATCH_GAPS,IGNORE_STRAND)==NULL)==SUBSET_NONOVERLAPS) qreg->Print(); 
      PRG.Check();
    }
    PRG.Done();
  }


  // clean up
  delete cmd_line;
  
  return 0;
}



