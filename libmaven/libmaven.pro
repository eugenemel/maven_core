include(../libmaven.pri)
DESTDIR = $$OUTPUT_DIR/lib

TEMPLATE=lib
CONFIG += staticlib
TARGET = maven

#Issue 321: for lipidsummarizationutils.h and lipidsummarizationutils.cpp
QT += core

LIBS += -L. -lcdfread -lnetcdf -lmstoolkitlite

DEFINES += _LARGEFILE_SOURCE _FILE_OFFSET_BITS=64 GCC

INCLUDEPATH += ../pugixml/src/ ../libcdfread/ ../zlib/ ../MSToolkit/include

SOURCES=base64.cpp mzMassCalculator.cpp mzSample.cpp  mzUtils.cpp statistics.cpp mzFit.cpp mzAligner.cpp\
       PeakGroup.cpp EIC.cpp Scan.cpp Peak.cpp  \
       Compound.cpp \
       savgol.cpp \
       SavGolSmoother.cpp \
       parallelmassSlicer.cpp \
       PolyAligner.cpp \ 
       Fragment.cpp \
       BondBreaker.cpp \
       Peptide.cpp \
       sha1.cpp \
       ThreadSafeSmoother.cpp \
       directinfusionprocessor.cpp \
       lipidsummarizationutils.cpp


HEADERS += base64.h mzFit.h mzMassCalculator.h mzSample.h mzPatterns.h mzUtils.h  statistics.h SavGolSmoother.h PolyAligner.h Fragment.h parallelmassSlicer.h BondBreaker.h Peptide.h sha1.hpp \
    ThreadSafeSmoother.h \
    directinfusionprocessor.h \
    lipidsummarizationutils.h

message($$LIBS)
