include(../libmaven.pri)

TEMPLATE = app
TARGET = mzDeltas
DESTDIR = ../bin/
INCLUDEPATH += ../MSToolkit/include/ ../pugixml/src .. ../libmaven

SOURCES = mzDeltas.cpp

LIBS += -lmaven -lpugixml -lneural -lz -lmstoolkitlite
