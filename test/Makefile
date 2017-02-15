CC=c++
CFLAGS=-O3 -std=c++11

all: test_peptide digest

test_eigen: test_eigen.cpp
	$(CC) $(CFLAGS) -O3 -o test_eigen  test_eigen.cpp -I../

test_peptide: test_peptide.cpp ../libmaven/Peptide.cpp
	$(CC) $(CFLAGS) -O3 -o test_peptide  ../libmaven/Peptide.cpp test_peptide.cpp $(INCLUDES)

digest: digest.cpp ../libmaven/Peptide.cpp
	$(CC) $(CFLAGS) -O3 -o digest ../libmaven/Peptide.cpp digest.cpp -I ../libmaven/
