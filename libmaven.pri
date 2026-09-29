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

    # Issue 853: on machines whose active Xcode SDK no longer ships
    # AGL.framework, Qt's own .prl files for internal QtGui/QtWidgets support
    # libraries still reference -framework AGL, breaking the final link step
    # of every Qt Widgets target. Fall back to the AGL.tbd stub still shipped
    # in an installed CommandLineTools SDK, if any (resolved dynamically so
    # this keeps working across CommandLineTools updates; a no-op if AGL
    # already links fine or no such SDK is present).
    AGL_FALLBACK_FRAMEWORKS = $$system(ls -d /Library/Developer/CommandLineTools/SDKs/MacOSX*.sdk/System/Library/Frameworks/AGL.framework 2>/dev/null | sort -V | tail -n 1 | xargs dirname 2>/dev/null)
    !isEmpty(AGL_FALLBACK_FRAMEWORKS) {
        message("Using AGL.framework fallback search path: $$AGL_FALLBACK_FRAMEWORKS")
        QMAKE_LFLAGS += -F$$AGL_FALLBACK_FRAMEWORKS
    }
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
