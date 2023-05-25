#include "isotopicenvelopeutils.h"

//private constructor: this is a singleton
IsotopeProcessorOptions::IsotopeProcessorOptions(){}

IsotopeProcessorOptions& IsotopeProcessorOptions::instance() {
    static IsotopeProcessorOptions options;
    return options;
}

void IsotopeProcessorOptions::setOptions(string config_file) {
    this->config_file = config_file;

    //TODO: import config file
}
