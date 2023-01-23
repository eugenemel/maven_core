OUTPUT_DIR = $$(OUTPUT_DIR)
isEmpty(OUTPUT_DIR):OUTPUT_DIR=$$PWD/build

INSTALL_LIBDIR = $$(INSTALL_LIBDIR)
unix {
  !mac {
    isEmpty(INSTALL_LIBDIR):INSTALL_LIBDIR=lib
} }

INSTALL_PREFIX=$$(DESTDIR)$$INSTALL_PREFIX
DEFINES += INSTALL_LIBDIR=\\\"$$INSTALL_LIBDIR\\\"
DEFINES += _LARGEFILE_SOURCE _FILE_OFFSET_BITS=64 GCC

QMAKE_CXXFLAGS_RELEASE   += -O3 -std=c++11
QMAKE_CXXFLAGS_DEBUG   += -g -std=c++11
QMAKE_CXXFLAGS_WARN_ON = -Wall  -Wno-sign-conversion -Wno-sign-compare



QT += core
CONFIG += silent exceptions std++14
OBJECTS_DIR = tmp
MOC_DIR = tmp
UI_DIR   =  tmp
QMAKE_CC = gcc
QMAKE_CXX = g++


# Issue 600/601: decoding issues
DEFINES += ZLIB
LIBS += -lz

win32-g++:contains(QMAKE_HOST.arch, x86_64):{
    DEFINES -= CDFPARSER
    LIBS -= -lcdfread -lnetcdf
}

win32 {
    message("using win32 config")
    DEFINES += MINGW
    DEFINES += WIN32
}

mac {
    message("using mac config")
	QMAKE_CXXFLAGS_WARN_ON += -Wno-c++11-extensions
    CONFIG+=sdk_no_version_check

    DEFINES -= CDFPARSER
    LIBS -= -lcdfread -lnetcdf
}

unix {
    message("using unix config")
    DEFINES -= LITTLE_ENDIAN
 #  DEFINES += CDFPARSER
 #  LIBS += -lcdfread -lnetcdf
 #  LIBS += -lz -lcdfread -lnetcdf
}


INCLUDEPATH += $$PWD

LIBS += -L$$OUTPUT_DIR/lib -L$$OUTPUT_DIR/plugin
