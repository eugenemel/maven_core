include(../libmaven.pri)
DESTDIR = $$OUTPUT_DIR/lib

TEMPLATE=lib
CONFIG += staticlib
TARGET = mstoolkitlite
LIBS += -L. -lz -lpthread -lm -ldl

DEFINES += _LARGEFILE_SOURCE _FILE_OFFSET_BITS=64 GCC HAVE_EXPAT_CONFIG_H _NOSQLITE
QMAKE_CXXFLAGS_RELEASE   -= -std=c++11
QMAKE_CFLAGS_RELEASE   -= -std=c++11


INCLUDEPATH += ./include/

SOURCES = src/MSToolkit/MSObject.cpp   src/MSToolkit/MSReader.cpp   src/MSToolkit/Spectrum.cpp   src/MSToolkit/mzMLWriter.cpp\
src/expat-2.0.1/xmlparse.c      src/expat-2.0.1/xmlrole.c       src/expat-2.0.1/xmltok.c        src/expat-2.0.1/xmltok_impl.c   src/expat-2.0.1/xmltok_ns.c\
src/expat-2.0.1/xmlparse.cpp    src/expat-2.0.1/xmlrole.cpp     src/expat-2.0.1/xmltok.cpp      src/expat-2.0.1/xmltok_impl.cpp src/expat-2.0.1/xmltok_ns.cpp\
./src/mzParser/BasicChromatogram.cpp ./src/mzParser/PWIZface.cpp          ./src/mzParser/mzParser.cpp          ./src/mzParser/saxhandler.cpp\
./src/mzParser/BasicSpectrum.cpp     ./src/mzParser/RAMPface.cpp          ./src/mzParser/mzpMz5Config.cpp      ./src/mzParser/saxmzmlhandler.cpp\
./src/mzParser/Czran.cpp             ./src/mzParser/mz5handler.cpp        ./src/mzParser/mzpMz5Structs.cpp     ./src/mzParser/saxmzxmlhandler.cpp\
./src/mzParser/MSNumpress.cpp        ./src/mzParser/mzMLReader.cpp        ./src/mzParser/mzp_base64.cpp

HEADERS = include/MSNumpress.hpp include/Spectrum.h\
    include/expat.h include/latin1tab.h\
    include/trees.h include/xmltok_impl.h\
    include/MSObject.h include/ascii.h\
    include/expat_config.h\
    include/mzMLWriter.h include/utf8tab.h\
    include/MSReader.h include/asciitab.h\
    include/expat_external.h include/mzParser.h\
    include/winconfig.h include/MSToolkitTypes.h\
    include/crc32.h include/nametab.h\
    include/xmlrole.h include/iasciitab.h\
    include/internal.h include/xmltok.h


