TEMPLATE = subdirs
CONFIG += ordered qt thread 
SUBDIRS += libneural libcdfread pugixml/src 

exists( MSToolkit)  { SUBDIRS += MSToolkit }
exists( libmaven ) { SUBDIRS += libmaven }
#exists( peakdetector ) { SUBDIRS += peakdetector }
#exists( maven_peakdetector ) { SUBDIRS += maven_peakdetector }
exists( maven ) { SUBDIRS += maven }
#exists( mzLibraryBrowser )  { SUBDIRS += mzLibraryBrowser }

