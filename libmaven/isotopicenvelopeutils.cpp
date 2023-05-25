#include "isotopicenvelopeutils.h"

//private constructor: this is a singleton
IsotopeProcessorOptions::IsotopeProcessorOptions(){}

IsotopeProcessorOptions& IsotopeProcessorOptions::instance() {
    static IsotopeProcessorOptions options;
    return options;
}
