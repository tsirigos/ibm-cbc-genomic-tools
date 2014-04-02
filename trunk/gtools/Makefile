
.SUFFIXES: .o .cpp
.cpp.o:
	$(CC) -c $(INCLUDE) $(CFLAGS) $<


# Macros

BINDIR = ../bin
CC = g++
CFLAGS = -g -O2 -Wall 
INCLUDE = -Isamtools
HEADER = gzstream.h core.h genomic_intervals.h 
OBJ = gzstream.o core.o genomic_intervals.o 
LIB = -Lsamtools -lbam -lm -lz -lgslcblas -lgsl 


# Explicit rules


gtools: genomic_regions genomic_overlaps genomic_scans genomic_apps gtools_hic vectors matrix permutation_test

clean:
	rm -f *.o
	rm -rf $(BINDIR)/*.exe $(BINDIR)/*.dSYM $(BINDIR)/genomic_overlaps $(BINDIR)/genomic_regions $(BINDIR)/genomic_scans $(BINDIR)/genomic_apps $(BINDIR)/gtools_hic $(BINDIR)/matrix $(BINDIR)/permutation_test $(BINDIR)/vectors
	cd samtools; make clean

libbam.a: 
	cd samtools; make libbam.a; cd ..
genomic_regions: $(HEADER) $(OBJ) libbam.a
	$(CC) $(CFLAGS) $(INCLUDE) -o $(BINDIR)/genomic_regions genomic_regions.cpp $(OBJ) $(LIB)
genomic_overlaps: $(HEADER) $(OBJ) libbam.a
	$(CC) $(CFLAGS) $(INCLUDE) -o $(BINDIR)/genomic_overlaps genomic_overlaps.cpp $(OBJ) $(LIB)
genomic_scans: $(HEADER) $(OBJ) libbam.a
	$(CC) $(CFLAGS) $(INCLUDE) -o $(BINDIR)/genomic_scans genomic_scans.cpp $(OBJ) $(LIB)
genomic_apps: $(HEADER) $(OBJ) libbam.a
	$(CC) $(CFLAGS) $(INCLUDE) -o $(BINDIR)/genomic_apps genomic_apps.cpp $(OBJ) $(LIB)
gtools_hic: $(HEADER) $(OBJ) libbam.a
	$(CC) $(CFLAGS) $(INCLUDE) -o $(BINDIR)/gtools_hic genomic_apps.cpp $(OBJ) $(LIB)
vectors: $(HEADER) $(OBJ) libbam.a
	$(CC) $(CFLAGS) $(INCLUDE) -o $(BINDIR)/vectors vectors.cpp $(OBJ) $(LIB)
matrix: $(HEADER) $(OBJ) libbam.a
	$(CC) $(CFLAGS) $(INCLUDE) -o $(BINDIR)/matrix matrix.cpp $(OBJ) $(LIB)
permutation_test: $(HEADER) $(OBJ) libbam.a
	$(CC) $(CFLAGS) $(INCLUDE) -o $(BINDIR)/permutation_test permutation_test.cpp $(OBJ) $(LIB)



# Implicit rules
gzstream.o: gzstream.h gzstream.c
core.o: core.h core.cpp
genomic_intervals.o: gzstream.h gzstream.o core.h core.cpp genomic_intervals.h genomic_intervals.cpp


