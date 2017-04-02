TEMPLATE = subdirs
CONFIG += ordered qt thread 
SUBDIRS += libneural libcdfread pugixml/src libmaven

exists( peakdetector ) { SUBDIRS += peakdetector }
exists( maven ) { SUBDIRS += maven }
exists( mzLibraryBrowser )  { SUBDIRS += mzLibraryBrowser }
exists( MSToolkit)  { SUBDIRS += MSToolkit }

