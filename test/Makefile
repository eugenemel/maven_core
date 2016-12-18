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

INCLUDES=-I../libmzorbi/ -I../pugixml/src/ -I../
LDFlAGS=-L/usr/lib/ -L ../build/lib/ -lmzorbi -lpugixml -lz -lsqlite3 # -lnetcdf -lcdfread

#INCLUDES=-I../libmzorbi/ -I../pugixml/src/ -I ./sqlite/
#OPENMP=-fopenmp


all: buildExclusionList mzFileInfo test_frag test_eigen mzCSVConvert matchIsotpicPattern DO2Quant

test_eigen: test_eigen.cpp
	$(CC) -O3 -o test_eigen  test_eigen.cpp -I../

mzCSVConvert: mzCSVConvert.cpp
	$(CC) $(CFLAGS) -o mzCSVConvert mzCSVConvert.cpp $(INCLUDES) $(LDFlAGS)

clusterMS2: clusterMS2.cpp
	$(CC) $(CFLAGS) -o clusterMS2 clusterMS2.cpp $(INCLUDES) $(LDFlAGS)

test_frag: testFragMatch.cpp
	$(CC) $(CFLAGS) -o test_frag testFragMatch.cpp $(INCLUDES) $(LDFlAGS)

mzFileInfo: mzFileInfo.cpp
	$(CC) $(CFLAGS) -o mzFileInfo mzFileInfo.cpp $(INCLUDES) $(LDFlAGS)

DO2Quant: DO2Quant.cpp Peptide.hpp
	$(CC) $(CFLAGS) -o DO2Quant DO2Quant.cpp Peptide.cpp	$(INCLUDES) $(LDFlAGS)

LysQuant: LysQuant.cpp Peptide.hpp
	$(CC) $(CFLAGS) -o LysQuant LysQuant.cpp Peptide.cpp	$(INCLUDES) $(LDFlAGS)


mzImage: mzImage.cpp ../libmzorbi/parallelmassSlicer.cpp
	$(CC) $(CFLAGS) -o mzImage mzImage.cpp ../libmzorbi/parallelmassSlicer.cpp 	$(INCLUDES) $(LDFlAGS)


matchIsotpicPattern: matchIsotpicPattern.cpp
	$(CC) $(CFLAGS) -o matchIsotpicPattern  matchIsotpicPattern.cpp $(INCLUDES) $(LDFlAGS)


buildExclusionList: buildExclusionList.cpp
	$(CC) $(CFLAGS) -o buildExclusionList buildExclusionList.cpp $(INCLUDES) $(LDFlAGS)




