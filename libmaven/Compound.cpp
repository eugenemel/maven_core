#include "Fragment.h"
#include "mzSample.h"

Compound::Compound(string id, string name, string formula, int charge ) {
		this->id = id;
		this->name = name;
		this->formula = formula;
		this->charge = charge;
        this->exactMass =  mcalc->computeNeutralMass(formula);
		this->expectedRt = -1;
        this->logP = 0;
        this->isDecoy=false;

        cid=0;
        ionizationMode=0;
		precursorMz=0;
		productMz=0;
        collisionEnergy=0;
        virtualFragmentation=false;

}

MassCalculator* Compound::mcalc = new MassCalculator();

float Compound::ajustedMass(int charge) { 
	return Compound::mcalc->computeMass(formula,charge); 
}


FragmentationMatchScore Compound::scoreCompoundHit(Fragment* f, float productPpmTolr=20, bool searchProton=true) {
        FragmentationMatchScore s;
        Compound* cpd = this;

        if (!cpd or cpd->fragment_mzs.size() == 0) return s;

        Fragment t;
        t.precursorMz = cpd->precursorMz;
        t.mzs = cpd->fragment_mzs;
        t.intensity_array = cpd->fragment_intensity;
        t.annotations = cpd->fragment_iontype;

        if (searchProton)  { //special case, check for loss or gain of protons
            int N = t.mzs.size();
            for(int i=0; i<N;i++) {
                t.mzs.push_back( t.mzs[i] + PROTON);
                t.intensity_array.push_back( t.intensity_array[i] );
                t.mzs.push_back( t.mzs[i] - PROTON );
                t.intensity_array.push_back( t.intensity_array[i] );
            }
        }

        //theory fragmentation or library fragmentation = t
        //experimental data = f
        s = t.scoreMatch(f,productPpmTolr);

        return s;
}
