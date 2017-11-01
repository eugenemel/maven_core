CC=c++
#CFLAGS=-O3 -std=c++11
CFLAGS=-O2 -std=c++11
MSTOOLKIT = ../MSToolkit
#CXXFLAGS += -O3 -Wall -Wextra -Wno-write-strings -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC -I. -I ../MSToolkit/include/ -L ../MSToolkit -lmstoolkitlite

all: formulaFitter peptide_ions digest mstoolkit

formulaFitter:
	$(CC) $(CFLAGS) -O3 -o formulaFitter formulaFitter.cpp -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC -I../MSToolkit/include -I../libmaven/ -I../pugixml/src  -L../build/lib -lmaven -lpugixml -lmstoolkitlite -lz

test_eigen: test_eigen.cpp
	$(CC) $(CFLAGS) -O3 -o test_eigen  test_eigen.cpp -I../

peptide_ions: peptide_ions.cpp ../libmaven/Peptide.cpp
	$(CC) $(CFLAGS) -O3 -o peptide_ions  ../libmaven/Peptide.cpp peptide_ions.cpp -I ../libmaven/

digest: digest.cpp ../libmaven/Peptide.cpp
	$(CC) $(CFLAGS) -O3 -o digest ../libmaven/Peptide.cpp digest.cpp -I ../libmaven/

mstoolkit: mstoolkit.cpp
	$(CC) $(CFLAGS) -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC  -o mstoolkit mstoolkit.cpp -I../MSToolkit/include -I../libmaven/ -I../pugixml/src  -L../build/lib -lmaven -lpugixml -lmstoolkitlite -lz
