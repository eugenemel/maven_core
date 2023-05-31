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

double IsotopicEnvelope::getIntensity() {
    if (totalIntensity < 0) {
        totalIntensity = std::accumulate(isotopeIntensity.begin(), isotopeIntensity.end(), 0.0);
    }
    return totalIntensity;
}

void IsotopicEnvelope::print() {
    if (!group) {
        cout << "Peak group is null." << endl;
        return;
    }

    if (!group->compound) {
        cout << "Compound is null." << endl;
        return;
    }

    if (!group->adduct) {
        cout << "Adduct is null." << endl;
        return;
    }

    cout << group->compound->name << " "
         << group->adduct->name << "@ rt="
         << group->medianRt() << ", totalIntensity="
         << getIntensity()
         << ":\n";

    for (unsigned int i = 0; i < isotopes.size(); i++) {
        Isotope isotope = isotopes.at(i);
        double intensity = isotopeIntensity.at(i);

        cout << "\t" << isotope.name << ": "
             << intensity << " "
             << "frac=" << (intensity/totalIntensity) << " "
             << "\n";
    }

    cout << endl;
}
