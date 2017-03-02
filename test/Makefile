CC=c++
CFLAGS=-O3 -std=c++11

all: peptide_ions digest

test_eigen: test_eigen.cpp
	$(CC) $(CFLAGS) -O3 -o test_eigen  test_eigen.cpp -I../

peptide_ions: peptide_ions.cpp ../libmaven/Peptide.cpp
	$(CC) $(CFLAGS) -O3 -o peptide_ions  ../libmaven/Peptide.cpp peptide_ions.cpp -I ../libmaven/

digest: digest.cpp ../libmaven/Peptide.cpp
	$(CC) $(CFLAGS) -O3 -o digest ../libmaven/Peptide.cpp digest.cpp -I ../libmaven/
