TEMPLATE = subdirs
CONFIG += ordered qt thread
SUBDIRS += libneural libcdfread pugixml/src

exists( MSToolkit)  { SUBDIRS += MSToolkit }
exists( libmaven ) { SUBDIRS += libmaven }
