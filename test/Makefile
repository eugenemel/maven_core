CC=c++

#OSX FLAGS
#CFLAGS=-g -Wall -Wno-sign-compare -std=c++0x
#LINUX FLAGS
CFLAGS=-O3 -std=c++11
#CFLAGS=-O3 -std=c++11 -stdlib=libstdc++
#CFLAGS=-O3 -stdlib=libstdc++ -Wno-c++11-extensions
#CFLAGS=-O3 -stdlib=libstdc++ -Wno-c++11-extensions
#CFLAGS=-O3 -Wno-c++11-extensions

DEFINES=-DLITTLE_ENDIAN

INCLUDES=-I../libmaven/ -I../pugixml/src/ -I../
LDFlAGS=-L/usr/lib/ -L ../build/lib/ -lmzorbi -lpugixml -lz -lsqlite3 # -lnetcdf -lcdfread

#INCLUDES=-I../libmzorbi/ -I../pugixml/src/ -I ./sqlite/
#OPENMP=-fopenmp


all: test_peptide

test_eigen: test_eigen.cpp
	$(CC) $(CFLAGS) -O3 -o test_eigen  test_eigen.cpp -I../

test_peptide: test_peptide.cpp ../libmaven/Peptide.cpp
	$(CC) $(CFLAGS) -O3 -o test_peptide  ../libmaven/Peptide.cpp test_peptide.cpp $(INCLUDES)
