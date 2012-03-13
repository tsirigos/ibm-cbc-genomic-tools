
.SUFFIXES: .o .cpp
.cpp.o:
	$(CC) -c $(INCLUDE) $(CFLAGS) $<


# Macros

BINDIR = ../bin
CC = g++
CFLAGS = -g -O2 -Wall 
LIB = -lm -lgslcblas -lgsl 


# Explicit rules

gtools: genomic_regions genomic_overlaps genomic_scans genomic_apps vectors matrix permutation_test

clean:
	rm -f *.o
	rm -f $(BINDIR)/*.exe $(BINDIR)/genomic_overlaps $(BINDIR)/genomic_regions $(BINDIR)/genomic_scans $(BINDIR)/genomic_apps $(BINDIR)/matrix $(BINDIR)/permutation_test $(BINDIR)/vectors

genomic_regions: core.h core.o genomic_intervals.h genomic_intervals.o 
	$(CC) $(CFLAGS) $(INCLUDE) -o $(BINDIR)/genomic_regions genomic_regions.cpp core.o genomic_intervals.o $(LIB)
genomic_overlaps: core.h core.o genomic_intervals.h genomic_intervals.o 
	$(CC) $(CFLAGS) $(INCLUDE) -o $(BINDIR)/genomic_overlaps genomic_overlaps.cpp core.o genomic_intervals.o $(LIB)
genomic_scans: core.h core.o genomic_intervals.h genomic_intervals.o
	$(CC) $(CFLAGS) $(INCLUDE) -o $(BINDIR)/genomic_scans genomic_scans.cpp core.o genomic_intervals.o $(LIB)
genomic_apps: core.h core.o genomic_intervals.h genomic_intervals.o
	$(CC) $(CFLAGS) $(INCLUDE) -o $(BINDIR)/genomic_apps genomic_apps.cpp core.o genomic_intervals.o $(LIB)
vectors: core.h core.o
	$(CC) $(CFLAGS) $(INCLUDE) -o $(BINDIR)/vectors vectors.cpp core.o $(LIB)
matrix: core.h core.o
	$(CC) $(CFLAGS) $(INCLUDE) -o $(BINDIR)/matrix matrix.cpp core.o $(LIB)
permutation_test: core.h core.o
	$(CC) $(CFLAGS) $(INCLUDE) -o $(BINDIR)/permutation_test permutation_test.cpp core.o $(LIB)



# Implicit rules
core.o: core.h core.cpp
genomic_intervals.o: core.h core.cpp genomic_intervals.h genomic_intervals.cpp


