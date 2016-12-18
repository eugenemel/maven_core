include(../../libmaven.pri)
DESTDIR = $$OUTPUT_DIR/lib

TEMPLATE = lib
CONFIG += staticlib exceptions
TARGET = pugixml

QMAKE_CXXFLAGS_RELEASE += -O9


SOURCES=pugixml.cpp  pugixpath.cpp
HEADERS=pugixml.hpp  pugiconfig.hpp

contains(MEEGO_EDITION,harmattan) {
    target.path = /opt/src/lib
    INSTALLS += target
}
