#include "mzSample.h"

//global options
int mzSample::filter_minIntensity = -1;
bool mzSample::filter_centroidScans = false;
int mzSample::filter_intensityQuantile = 0;
int mzSample::filter_polarity=0;
int mzSample::filter_mslevel=0;


mzSample::mzSample() {
    maxMz = maxRt = 0;
    minMz = minRt = 0;
    isBlank = false;
    isSelected=true;
    maxIntensity=0;
    minIntensity=0;
    totalIntensity=0;
    _normalizationConstant=1;
    _sampleOrder=-1;
    _C13Labeled=false;
    _N15Labeled=false;
    _S34Labeled=false; //Feng note: added to track S34 labeling state
    _D2Labeled=false; //Feng note: added to track D2 labeling state
    _setName =  "A";
    sampleId = -1;
    color[0]=color[1]=color[2]=0;
    color[3]=1.0;
    _averageFullScanTime = -1.0f; //initial value of -1 indicates that computation not yet performed from mzSample::getAverageFullScanTime()
}

mzSample::~mzSample() { 
        for(unsigned int i=0; i < scans.size(); i++ )
            if(scans[i]!=NULL) delete(scans[i]);
        scans.clear();
}

void mzSample::addScan(Scan*s ) {
        if (!s) return;

    //cerr << "addScan=" << "\tms="<< s->mslevel << "\tprecMz=" << s->precursorMz << "\trt=" << s->rt << endl;

	//skip scans that do not match mslevel
	if (mzSample::filter_mslevel and s->mslevel != mzSample::filter_mslevel ) {
		return;
	}
	//skip scans that do not match polarity 
	if (mzSample::filter_polarity and s->getPolarity() != mzSample::filter_polarity ) {
		return;
	}

        //unsigned int sizeBefore = s->intensity.size();
        if ( mzSample::filter_centroidScans == true ) {
            s->simpleCentroid();
        }

        //unsigned int sizeAfter1 = s->intensity.size();

        if ( mzSample::filter_intensityQuantile > 0) {
            s->quantileFilter(mzSample::filter_intensityQuantile);
        }
        //unsigned int sizeAfter2 = s->intensity.size();

        if ( mzSample::filter_minIntensity > 0) {
            s->intensityFilter(mzSample::filter_minIntensity);
        }
        //unsigned int sizeAfter3 = s->intensity.size();
        //cerr << "addScan " << sizeBefore <<  " " << sizeAfter1 << " " << sizeAfter2 << " " << sizeAfter3 << endl;

        scans.push_back(s);
        s->scannum=scans.size()-1;
}

string mzSample::cleanSampleName(string sampleName) {
        size_t pos =sampleName.find_last_of("/");
        if (pos != std::string::npos) {
                sampleName=sampleName.substr(pos+1, sampleName.length());
        }

        pos=sampleName.find_last_of("\\");
        if (pos != std::string::npos) {
                sampleName=sampleName.substr(pos+1, sampleName.length());
        }
        return sampleName;

}

void mzSample::loadSample(const char* filename) {
    loadSample(filename, true);
}

void mzSample::loadSample(const char* filename, bool isCorrectPrecursor) {

    string filenameString = string(filename);
    this->sampleName = cleanSampleName(filename);
    this->fileName = filenameString;

    if (mystrcasestr(filename,".mzCSV") != NULL ) {
        parseMzCSV(filename);
    } else if(mystrcasestr(filename,".mzdata") != NULL or mystrcasestr(filename,".mzdata.gz") != NULL ) {
        parseMzData(filename);
    } else if(mystrcasestr(filename,".mzxml") != NULL or mystrcasestr(filename,".mzxml.gz") != NULL ) {
        parseMzXML(filename);
    } else if(mystrcasestr(filename,".mzml") != NULL or mystrcasestr(filename,".mzml.gz") != NULL ) {
        parseMzML(filename);
    } else if(mystrcasestr(filename,".cdf") != NULL ) {
        parseCDF((char*) filename,1);
    } else {
        parseMzXML(filename);
    }

    //set min and max values for rt
    calculateMzRtRange();

    //recalculate precursor masses
    if (isCorrectPrecursor) {
        //cerr << "Recalculating Ms2 Precursor Masses" << endl;
        for(Scan* ms2scan: scans) {
            ms2scan->precursorMz=getMS1PrecursorMass(ms2scan,20);
        }
    }

    if (mystrcasestr(filename,"blan") != NULL) {
        this->isBlank = true;
        cerr << "Found Blank: " << filename << endl;
    }
}

void mzSample::loadMsToolsSample(const char* filename) {
using namespace MSToolkit;

    mzSample* sample = new mzSample();
    string filenameString = string(filename);
    this->sampleName = sample->cleanSampleName(filename);
    this->fileName = filenameString;

    MSToolkit::Spectrum spec;           // For holding spectrum.
    MSToolkit::MSReader mstReader;

   int iFirstScan=0;
   mstReader.readFile(filename, spec, iFirstScan);

   MSToolkit::MSSpectrumType filter = MSToolkit::MS1;
   mstReader.setFilter(filter);

   int iFileLastScan = mstReader.getLastScan();
   cerr << this->sampleName <<  " #scans=" << iFileLastScan << endl;

    for(int scanNum=0; scanNum<iFileLastScan; scanNum++) {

        //spec.clearPeaks();
        //spec.clearMZ();

        mstReader.readFile(NULL, spec,scanNum);

        //basic scan information
        Scan* scan = new Scan(this,
                scanNum,
                spec.getMsLevel(),
                spec.getRTime(), 0, 0);


        //precursor information
        if( spec.sizeMZ()) {
                scan->precursorMz = spec.getMZ(0);
                scan->precursorCharge = spec.getCharge();
        }

        //Activation method
        MSToolkit::MSActivation act = spec.getActivationMethod();
        string actMethod;
        switch(act){
            case MSToolkit::mstETD: actMethod="ETD"; break;
            case MSToolkit::mstETDSA: actMethod="ETDSA"; break;
            case MSToolkit::mstCID: actMethod="CID"; break;
            case MSToolkit::mstECD: actMethod="ECD"; break;
            case MSToolkit::mstPQD: actMethod = "PQD"; break;
            case MSToolkit::mstHCD: actMethod = "HCD"; break;
            case MSToolkit::mstNA: default: actMethod="UNKNOWN";
            break;
        }
        scan->activationMethod = actMethod;

        if (spec.size() ) { 
			for(int i=0; i < spec.size(); i++ ) {
				scan->intensity.push_back( spec[i].intensity );
				scan->mz.push_back( spec[i].mz );
            }
        }

		/*
        cerr << scan->scannum 
            << "\t" << scan->mslevel
            << "\t" << scan->totalIntensity()
            << "\t" << scan->precursorMz
            << "\t" << scan->precursorCharge
            << "\t" << scan->activationMethod
            << endl;
			*/

        this->addScan(scan);
    }

    //set min and max values for rt
   this->calculateMzRtRange();

    //recalculate precursor masses
    //cerr << "Recalculating Ms2 Precursor Masses" << endl;
    for(Scan* ms2scan: this->scans) {
        if(ms2scan->mslevel==2) {
            ms2scan->precursorMz=this->getMS1PrecursorMass(ms2scan,20);
        }
    }

    if (mystrcasestr(filename,"blan") != NULL) {
        this->isBlank = true;
        cerr << "Found Blank: " << filename << endl;
    }
}

void mzSample::parseMzCSV(const char* filename) {
		// file structure: scannum,rt,mz,intensity,mslevel,precursorMz,polarity,srmid
		int lineNum=0; 
		cerr << "Loading " << filename << endl;
		ifstream myfile(filename);
		if (! myfile.is_open() ) { cerr << "Can't open file " << filename; return; }

	    std::stringstream ss;
            std::string line;
            std::string polarity;

			int lastScanNum=-1;
			int scannum=0;
			float rt=0;
			float intensity=0;
			float mz=0;
			float precursorMz=0;
			int mslevel=0;

	   Scan* scan = NULL;
            int newscannum=0;

            while ( getline(myfile,line) ) {
                lineNum++;
                vector<string>fields;
                mzUtils::split(line,',', fields); 
                if (fields.size() >= 5 && lineNum > 1) {

                    ss.clear();

                    ss  << fields[0] <<  " " 
                        << fields[1] <<  " "
                        << fields[2] <<  " " 
                        << fields[3] <<  " "
                        << fields[4] << " "
                        << fields[5] << " "
                        << fields[6];

                    ss >> scannum >> rt >> mz >> intensity >> mslevel >> precursorMz >> polarity;

                    if ( scannum != lastScanNum ) {
                        newscannum++;
                        if (mslevel <= 0 ) mslevel=1;
                        int scanpolarity=0; 
                        if (polarity.empty() && fields.size()>7) polarity=fields[7];
                        if (!polarity.empty() && polarity[0] == '+' ) scanpolarity=1; 
                        if (!polarity.empty() && polarity[0] == '-' ) scanpolarity=-1; 
                        scan = new Scan(this,newscannum,mslevel,rt/60,precursorMz,scanpolarity);
                        if (mslevel > 1 ) scan->productMz=mz;

                        addScan(scan);
                        if (fields.size() > 7) scan->filterLine=fields[7];		//last field is srmId
                    }

                    scan->mz.push_back(mz);
                    scan->intensity.push_back(intensity);
                    lastScanNum = scannum;
                }
            }
}


void mzSample::writeMzCSV(const char* filename) {
       ofstream mzCSV;
       mzCSV.open(filename);
	   if (!mzCSV.is_open()) { cerr << "Unable to write to a file" << filename; return; }
		 
       mzCSV << "scannum,rt,mz,intensity,mslevel,precursorMz,polarity,srmid" << endl;
       for(unsigned int i=0;i<scans.size();i++ ) {
           Scan* scan = scans[i];
           for(unsigned int j=0; j< scan->nobs(); j++ ) {
               mzCSV << scan->scannum+1 << ","
               << scan->rt*60   << ","
               << scan->mz[j]   << ","
               << scan->intensity[j]   << ","
               << scan->mslevel   << ","
               << scan->precursorMz   << ","
               << (scan->getPolarity() > 0 ? "+" : "-" ) << ","
               << scan->filterLine << endl;
           }
       }
}

void mzSample::parseMzML(const char* filename) { 
    xml_document doc;

    const unsigned int parse_options = parse_minimal;

    bool loadok = doc.load_file(filename, parse_options);
    if (!loadok ) {
        cerr << "Failed to load " << filename << endl;
        return;
    }

    //Get a spectrumstore node
    xml_node chromatogramList = doc.first_child().first_element_by_path("mzML/run/chromatogramList");
    xml_node spectrumList = doc.first_child().first_element_by_path("mzML/run/spectrumList");

    if (!spectrumList.empty()) { parseMzMLSpectrumList(spectrumList); }
    else if (!chromatogramList.empty()) { parseMzMLChromatogromList(chromatogramList); }
}


void mzSample::parseMzMLChromatogromList(xml_node chromatogramList) {

    int scannum=0;
    for (xml_node chromatogram  = chromatogramList.child("chromatogram");
                                        chromatogram; chromatogram = chromatogram.next_sibling("chromatogram")) {

                        string chromatogramId = chromatogram.attribute("id").value();
                        vector<float> timeVector;
                        vector<float> intsVector;

                        xml_node binaryDataArrayList = chromatogram.child("binaryDataArrayList");
                        string precursorMzStr = chromatogram.first_element_by_path("precursor/isolationWindow/cvParam").attribute("value").value();
                        string productMzStr = chromatogram.first_element_by_path("product/isolationWindow/cvParam").attribute("value").value();

                        float precursorMz = -1.0f;
                        float productMz = -1.0f;

                        try {
                            precursorMz = string2float(precursorMzStr);
                            productMz = string2float(productMzStr);
                        } catch (...) {
                            //swallow, skip these scans
                        }

                        //Issue 347: retrieve scans
                        int scanpolarity = 0;
                        map<string,string>cvParams = mzML_cvParams(chromatogram);

                        if (cvParams.count("positive scan")) {
                            scanpolarity=1;
                        } else if(cvParams.count("negative scan")) {
                            scanpolarity=-1;
                        }

                        int mslevel=2;

        for( xml_node binaryDataArray= binaryDataArrayList.child("binaryDataArray");
                                binaryDataArray; binaryDataArray=binaryDataArray.next_sibling("binaryDataArray")) {
                                map<string,string>attr = mzML_cvParams(binaryDataArray);

                                int precision = 64;
                                if(attr.count("32-bit float")) precision=32;

								bool decompress = false;
								if(attr.count("zlib compression")) decompress=true;

                                string binaryDataStr = binaryDataArray.child("binary").child_value();
                                vector<float>binaryData = base64::decode_base64(binaryDataStr,precision/8,false,decompress);

                                if(attr.count("time array")) { timeVector = binaryData; }
                                if(attr.count("intensity array")) { intsVector = binaryData; }


        }

//                cerr << chromatogramId << endl;
//                cerr << timeVector.size()  << " ints=" << intsVector.size() << endl;
//                cerr << "pre: " << precursorMz  << " prod=" << productMz << endl;

                if (precursorMz > 0 and productMz > 0) {
                                for(unsigned int i=0; i < timeVector.size(); i++ ) {
                                                Scan* scan = new Scan(this,scannum++,mslevel,timeVector[i],precursorMz,scanpolarity);
                                                scan->productMz=productMz;
                                                scan->mz.push_back(productMz);
                                                scan->filterLine= chromatogramId;
                                                scan->intensity.push_back(intsVector[i]);
                                                addScan(scan);
                                }
                }
    }

        //renumber scans based on retention time
        std::sort(scans.begin(), scans.end(),Scan::compRt);
        for(unsigned int i=0; i< scans.size(); i++ ) {
                        scans[i]->scannum=i+1;
        }

}

void mzSample::parseMzMLSpectrumList(xml_node spectrumList) {

    //Iterate through spectrums
    int scannum=0;

    for (xml_node spectrum  = spectrumList.child("spectrum"); spectrum; spectrum = spectrum.next_sibling("spectrum")) {

        string spectrumId = spectrum.attribute("id").value();
        //cerr << "Processing: " << spectrumId << endl;

        if (spectrum.empty()) continue;
        map<string,string>cvParams = mzML_cvParams(spectrum);

         int mslevel=1;
            int scanpolarity=0;
        float rt=0;
        vector<float> mzVector;
        vector<float> intsVector;

        if(cvParams.count("ms level")) {
            string msLevelStr = cvParams["ms level"];
            mslevel=(int) string2float(msLevelStr);
        }

        if(cvParams.count("positive scan")) scanpolarity=1;
        else if(cvParams.count("negative scan")) scanpolarity=-1;
        else scanpolarity=0;

        xml_node scanNode = spectrum.first_element_by_path("scanList/scan");
        map<string,string>scanAttr = mzML_cvParams(scanNode);
        if(scanAttr.count("scan start time")) {
            string rtStr = scanAttr["scan start time"];
            rt = string2float(rtStr);
        }

        float lowerLimitMz = -1.0f;
        float upperLimitMz = -1.0f;

        xml_node scanWindowListNode = scanNode.first_element_by_path("scanWindowList/scanWindow");
        map<string,string> scanWindowAttr = mzML_cvParams(scanWindowListNode);
        if (scanWindowAttr.count("scan window lower limit")) {
            lowerLimitMz = string2float(scanWindowAttr["scan window lower limit"]);
        }
        if (scanWindowAttr.count("scan window upper limit")) {
            upperLimitMz = string2float(scanWindowAttr["scan window upper limit"]);
        }

        string filterString = "";
        if (scanAttr.count("filter string")){
            filterString = scanAttr["filter string"];
        }

        string precursorMzStr;
        float precursorMz = 0;
        float ms1PrecursorForMs3=0;

        float precursorIsolationWindow=0;
        float isolationWindowLowerOffset=0;
        float isolationWindowUpperOffset=0;

        //Issue 256: Support targeted MS3 experiments
        if (mslevel == 3) {
            xml_node precursorList = spectrum.first_element_by_path("precursorList");
            for (xml_node precursor  = precursorList.child("precursor"); precursor; precursor = precursor.next_sibling("precursor")) {

                bool isMs1Precursor = false;

                map<string, string> selectedIon = mzML_cvParams(precursor.first_element_by_path("selectedIonList/selectedIon"));
                precursorMzStr = selectedIon["selected ion m/z"];

                map<string,string>isolationWindow = mzML_cvParams(precursor.first_element_by_path("isolationWindow"));
                if ( precursorMzStr.length() == 0 && isolationWindow.find("isolation window target m/z") != isolationWindow.end()){
                    precursorMzStr = isolationWindow["isolation window target m/z"];
                }

                xml_node isolationWindowNode = precursor.first_element_by_path("isolationWindow");

                for (xml_node cv = isolationWindowNode.child("userParam"); cv; cv = cv.next_sibling("userParam")) {
                    string name = cv.attribute("name").value();
                    string value = cv.attribute("value").value();
                    if (name == "ms level" && value == "1") {
                        isMs1Precursor = true;
                        break;
                    }
                }

                if (isolationWindow.find("isolation window lower offset") != isolationWindow.end()) {
                    string precursorIsolationStrLower = isolationWindow["isolation window lower offset"];
                    isolationWindowLowerOffset = string2float(precursorIsolationStrLower);
                    if(isolationWindowLowerOffset>0) precursorIsolationWindow+=string2float(precursorIsolationStrLower);
                }
                if (isolationWindow.find("isolation window upper offset") != isolationWindow.end()) {
                    string precursorIsolationStrUpper = isolationWindow["isolation window upper offset"];
                    isolationWindowUpperOffset = string2float(precursorIsolationStrUpper);
                    if(isolationWindowUpperOffset>0) precursorIsolationWindow+=string2float(precursorIsolationStrUpper);
                }

                if (precursorIsolationWindow <= 0) precursorIsolationWindow = 1.0; //default to 1.0

                if(string2float(precursorMzStr)>0){
                    if (isMs1Precursor) {
                        ms1PrecursorForMs3=string2float(precursorMzStr);
                    } else { //ms2 precursor, or the true precursor
                        precursorMz=string2float(precursorMzStr);
                    }
                }
            }
        } else { //no precursor (mslevel == 1) or a single precursor (mslevel == 2). ms level higher than 3 is not supported

            map<string,string>selectedIon = mzML_cvParams(spectrum.first_element_by_path("precursorList/precursor/selectedIonList/selectedIon"));

            precursorMzStr = selectedIon["selected ion m/z"];

            map<string,string>isolationWindow = mzML_cvParams(spectrum.first_element_by_path("precursorList/precursor/isolationWindow"));

            if ( precursorMzStr.length() == 0 && isolationWindow.find("isolation window target m/z") != isolationWindow.end()){
                precursorMzStr = isolationWindow["isolation window target m/z"];
            }

            if (isolationWindow.find("isolation window upper offset") != isolationWindow.end()) {
                string precursorIsolationStrLower = isolationWindow["isolation window lower offset"];
                isolationWindowLowerOffset = string2float(precursorIsolationStrLower);
                if(isolationWindowLowerOffset>0) precursorIsolationWindow+=string2float(precursorIsolationStrLower);
            }

            if (isolationWindow.find("isolation window upper offset") != isolationWindow.end()) {
                string precursorIsolationStrUpper = isolationWindow["isolation window upper offset"];
                isolationWindowUpperOffset = string2float(precursorIsolationStrUpper);
                if(isolationWindowUpperOffset>0) precursorIsolationWindow+=string2float(precursorIsolationStrUpper);
            }

            if (precursorIsolationWindow <= 0) precursorIsolationWindow = 1.0; //default to 1.0

            if(string2float(precursorMzStr)>0){
                precursorMz=string2float(precursorMzStr);
            }
        }

        float injectionTime=0;
        string injectionTimeStr = scanAttr["ion injection time"];
        if(string2float(injectionTimeStr)>0) injectionTime=string2float(injectionTimeStr);

        string productMzStr = spectrum.first_element_by_path("product/isolationWindow/cvParam").attribute("value").value();
        float productMz = 0;   if(string2float(productMzStr)>0)   productMz=string2float(productMzStr);

        xml_node binaryDataArrayList = spectrum.child("binaryDataArrayList");
        if( ! binaryDataArrayList or binaryDataArrayList.empty()) continue;

        for(xml_node binaryDataArray= binaryDataArrayList.child("binaryDataArray");
                binaryDataArray; binaryDataArray=binaryDataArray.next_sibling("binaryDataArray")) {
                if( ! binaryDataArray or binaryDataArray.empty()) continue;

                map<string,string>attr = mzML_cvParams(binaryDataArray);

            int precision = 64;
            if(attr.count("32-bit float")) precision=32;

			bool decompress = false;
            if(attr.count("zlib compression")) decompress=true;

            string binaryDataStr = binaryDataArray.child("binary").child_value();
            if (!binaryDataStr.empty()) {
                vector<float>binaryData = base64::decode_base64(binaryDataStr,precision/8,false,decompress);
                if(attr.count("m/z array")) { mzVector = binaryData; }
                if(attr.count("intensity array")) { intsVector = binaryData; }
            }
        }

        //cerr << " scan=" << scannum << "\tms="<< mslevel << "\tprecMz=" << precursorMz << "\t rt=" << rt << "\tisolWin=" << precursorIsolationWindow << endl;
        Scan* scan = new Scan(this,scannum++,mslevel,rt,precursorMz,scanpolarity);
		scan->isolationWindow = precursorIsolationWindow;
        scan->productMz=productMz;
        scan->filterLine= spectrumId;
        scan->intensity = intsVector;
        scan->injectionTime = injectionTime;
        scan->mz= mzVector;
        scan->filterString = filterString;
        scan->lowerLimitMz = lowerLimitMz;
        scan->upperLimitMz = upperLimitMz;
        scan->ms1PrecursorForMs3 = ms1PrecursorForMs3;

        if (isolationWindowLowerOffset>0) scan->isolationWindowLowerOffset = isolationWindowLowerOffset;
        if (isolationWindowUpperOffset>0) scan->isolationWindowUpperOffset = isolationWindowUpperOffset;

        addScan(scan);
    }
 }

map<string,string> mzSample::mzML_cvParams(xml_node node) { 
		map<string,string>attr;
        if(!node || node.empty()) return attr;
    	for (xml_node cv = node.child("cvParam"); cv; cv = cv.next_sibling("cvParam")) {
				string name = cv.attribute("name").value();
				string value = cv.attribute("value").value();
				attr[name]=value;
				//cerr << name << "->" << value << endl;
		}
		return(attr);
}



void mzSample::parseMzData(const char* filename) { 
    xml_document doc;

    const unsigned int parse_options = parse_minimal;

    bool loadok = doc.load_file(filename, parse_options);
    if (!loadok ) {
        cerr << "Failed to load " << filename << endl;
        return;
    }

    //Get a spectrumstore node
    xml_node spectrumstore = doc.first_child().child("spectrumList");

	if(!spectrumstore) { 
		spectrumstore = doc.first_element_by_path("ExperimentCollection/Experiment/mzData/spectrumList");
	}

    //Iterate through spectrums
    int scannum=0;

    for (xml_node spectrum = spectrumstore.child("spectrum"); spectrum; spectrum = spectrum.next_sibling("spectrum")) {
        scannum++;
        float rt = 0;
        float precursorMz = 0;
		int   precursorCharge=0;
        float collisionEnergy=0;
        char scanpolarity = 0;	//default case

        xml_node spectrumInstrument = spectrum.first_element_by_path("spectrumDesc/spectrumSettings/spectrumInstrument");
        int mslevel =  spectrumInstrument.attribute("msLevel").as_int();
        //cerr << mslevel << " " << spectrum.attribute("msLevel").value() << endl;

        for( xml_node cvParam= spectrumInstrument.child("cvParam"); cvParam; cvParam= cvParam.next_sibling("cvParam")) {
            //cout << "cvParam=" << cvParam.attribute("name").value() << endl;
            //
            if (strncasecmp(cvParam.attribute("name").value(),"TimeInMinutes",10) == 0 ) {
                rt=cvParam.attribute("value").as_float();
                //cout << "rt=" << rt << endl;
            }

            if (strncasecmp(cvParam.attribute("name").value(),"time in seconds",10) == 0 ) {
                rt=cvParam.attribute("value").as_float()/60;
                //cout << "rt=" << rt << endl;
            }

            if (strncasecmp(cvParam.attribute("name").value(),"polarity",5) == 0 ) {
                if ( cvParam.attribute("value").value()[0] == 'p' || cvParam.attribute("value").value()[0] == 'P') {
                    scanpolarity = +1;
                } else {
                    scanpolarity = -1;
                }

            }
        }

		xml_node ionSelection = spectrum.first_element_by_path("spectrumDesc/precursorList/precursor/ionSelection");
		for( xml_node cvParam= ionSelection.child("cvParam"); cvParam; cvParam= cvParam.next_sibling("cvParam")) {
			if (strncasecmp(cvParam.attribute("name").value(),"MassToChargeRatio",10) == 0 ) {
				precursorMz=cvParam.attribute("value").as_float();
			}

			if (strncasecmp(cvParam.attribute("accession").value(),"MS:1000744",10) == 0 ) {
				precursorMz=cvParam.attribute("value").as_float();
			}

			if (strncasecmp(cvParam.attribute("accession").value(),"MS:1000041",10) == 0 ) {
				precursorCharge=cvParam.attribute("value").as_float();
			}

		}

        xml_node activation = spectrum.first_element_by_path("spectrumDesc/precursorList/precursor/activation");
		for( xml_node cvParam= activation.child("cvParam"); cvParam; cvParam= cvParam.next_sibling("cvParam")) {
			if (strncasecmp(cvParam.attribute("name").value(),"CollisionEnergy",10) == 0 ) {
				collisionEnergy=cvParam.attribute("value").as_float();
			}
		}


        //cout << spectrum.first_element_by_path("spectrumDesc/spectrumSettings/spectrumInstrument").child_value() << endl
        if (mslevel <= 0 ) mslevel=1;
        Scan* scan = new Scan(this,scannum,mslevel,rt,precursorMz,scanpolarity);
        scan->collisionEnergy=collisionEnergy;
		scan->precursorCharge=precursorCharge;
        addScan(scan);

        int precision1 = spectrum.child("intenArrayBinary").child("data").attribute("precision").as_int();
        string b64intensity = spectrum.child("intenArrayBinary").child("data").child_value();
        scan->intensity = base64::decode_base64(b64intensity,precision1/8,false,false);

        //cout << "mz" << endl;
        int precision2 = spectrum.child("mzArrayBinary").child("data").attribute("precision").as_int();
        string b64mz = spectrum.child("mzArrayBinary").child("data").child_value();
        scan->mz = base64::decode_base64(b64mz,precision2/8,false,false);

        //cout << "spectrum " << spectrum.attribute("title").value() << endl;
    }
}

void mzSample::parseMzXML(const char* filename) { 
	    xml_document doc;
            try {
                bool loadok = doc.load_file(filename,pugi::parse_minimal);

                if (!loadok ) {
                    cerr << "Failed to load " << filename << endl;
                    return;
                }

                //Get a spectrumstore node
                xml_node spectrumstore = doc.first_child().child("msRun");
                if (spectrumstore.empty()) {
                    xml_node scan = doc.first_child().child("scan");
                    if(!scan.empty()) { spectrumstore=doc.first_child(); }  
                    else { cerr << "parseMzXML: can't find <msRun> or <scan> section" << endl; return; }
                }

                xml_node msInstrument = spectrumstore.child("msInstrument");
                if (!msInstrument.empty()) {
                    xml_node msManufacturer = msInstrument.child("msManufacturer");
                    xml_node msModel = msInstrument.child("msModel");
                    xml_node msIonisation = msInstrument.child("msIonisation");
                    xml_node msMassAnalyzer = msInstrument.child("msMassAnalyzer");
                    xml_node msDetector = msInstrument.child("msDetector");
                    instrumentInfo[ "msManufacturer" ] =  msManufacturer.attribute("value").value();
                    instrumentInfo[ "msModel" ] =  msModel.attribute("value").value();
                    instrumentInfo[ "msIonisation" ] =  msIonisation.attribute("value").value();
                    instrumentInfo[ "msMassAnalyzer" ] =  msMassAnalyzer.attribute("value").value();
                    instrumentInfo[ "msDetector" ] =  msDetector.attribute("value").value();
                }


                //Iterate through spectrums
                int scannum=0;

                for (xml_node scan = spectrumstore.child("scan"); scan; scan = scan.next_sibling("scan")) {
                    scannum++;
                    if (strncasecmp(scan.name(),"scan",4) == 0) {
                        Scan* newScan = parseMzXMLScan(scan,scannum);
                        if(newScan) addScan(newScan);
                    }

                    for ( xml_node child = scan.first_child(); child; child = child.next_sibling()) {
                        scannum++;
                        if (strncasecmp(child.name(),"scan",4) == 0) {
                            Scan* newScan = parseMzXMLScan(child,scannum);
                            if(newScan) addScan(newScan);
                        }
                    }
                }
            } catch(char* err) {
                cerr << "Failed to load file: " <<  filename << " " << err << endl;
            }
}

Scan* mzSample::parseMzXMLScan(const xml_node& mzxml_scan_node, int scannum) {

    float rt = 0;
    float precursorMz = 0;
    float productMz=0;
    float collisionEnergy=0;
    int scanpolarity = 0;	//default case
    int msLevel=1;
    bool networkorder = false;
	int centroided=0;
    string filterLine;
    string scanType;

    for(xml_attribute attr = mzxml_scan_node.first_attribute(); attr; attr=attr.next_attribute()) {
        if (strncasecmp(attr.name(), "retentionTime",10) == 0 ) {
            if (strncasecmp(attr.value(),"PT",2) == 0 ) {
                rt=string2float( string(attr.value()+2));
            } else if (strncasecmp(attr.value(),"P",1) == 0 ) {
                rt=string2float( string(attr.value()+1));
            } else {
                rt=string2float( attr.value());
            }
            rt /=60;
        }

        if (strncasecmp(attr.name(),"polarity",5) == 0 ) {
            char p = attr.value()[0];
            switch(p) {
            case '+': scanpolarity=  1; break;
            case '-': scanpolarity= -1; break;
            default:
                //cerr << "Warning:: scan has unknown polarity type=" << scanpolarity << endl;
                scanpolarity=0;
                break;
            }
        }

        if (strncasecmp(attr.name(),"filterLine",9) == 0) filterLine = attr.value();
        if (strncasecmp(attr.name(),"mslevel",5) == 0 ) msLevel = string2integer(attr.value());
        if (strncasecmp(attr.name(),"basePeakMz",9) == 0 ) productMz = string2float(attr.value());
        if (strncasecmp(attr.name(),"collisionEnergy",12) == 0 ) collisionEnergy= string2float(attr.value());
        if (strncasecmp(attr.name(),"scanType",8) == 0 ) scanType = attr.value();
        if (strncasecmp(attr.name(),"centroided",10) == 0 ) centroided = string2integer(attr.value());

    }

    //work around .. get polarity from filterline
    if (scanpolarity == 0 && filterLine.size()>13 ) {
        if ( filterLine[12] == '+' ) {
            scanpolarity = 1;
        }  else {
            scanpolarity = -1;
        }
    }

    if (msLevel <= 0 or msLevel > 100) msLevel=1;

    //add new scan
    Scan* _scan = new Scan(this,scannum,msLevel,rt,precursorMz,scanpolarity);

    // cerr << "newScan=" << msLevel << " " << precursorMz << endl;
    if(centroided==1) {
        _scan->centroided=true;
    } else {
        _scan->centroided=false;
    }


    xml_node precursor = mzxml_scan_node.child("precursorMz");
    if (precursor) {
        _scan->precursorMz = string2float(mzxml_scan_node.child_value("precursorMz"));
		_scan->isolationWindow = 1.0; //default value of isolation window is 1 da.

        for(xml_attribute attr = precursor.first_attribute(); attr; attr=attr.next_attribute()) {
            if (strncasecmp(attr.name(),     "precursorIntensity",15) == 0) _scan->precursorIntensity = string2float(attr.value());
            else if (strncasecmp(attr.name(),"precursorCharge",15) == 0) _scan->precursorCharge = string2integer(attr.value());
            else if (strncasecmp(attr.name(),"activationMethod",15) == 0) _scan->activationMethod = string(attr.value());
            else if (strncasecmp(attr.name(),"precursorScanNum",15) == 0) _scan->precursorScanNum = string2integer(attr.value());
            else if (strncasecmp(attr.name(),"windowWideness",13) == 0) _scan->isolationWindow = string2float(attr.value());
         }
    /*
        cerr << "mzLevel=" << _scan->mslevel << endl;
        cerr << "precursorMz=" << _scan->precursorMz << endl;
        cerr << "precursorCharge=" << _scan->precursorCharge << endl;
        cerr << "activationMethod=" << _scan->activationMethod << endl;
        cerr << "precursorIntensity=" << _scan->precursorIntensity << endl;
    */

    }
    //addScan(_scan);

    xml_node peaks =  mzxml_scan_node.child("peaks");
    if ( ! peaks.empty() ) {
        string b64String(peaks.child_value());
        if ( b64String.empty()) return _scan;  //no m/z intensity values

		bool decompress = false; //decompress?
        if(strncasecmp(peaks.attribute("compressionType").value(),"zlib",4) == 0) decompress=true;

        if (!filterLine.empty()) _scan->filterLine=filterLine;
        if (!scanType.empty())   _scan->scanType= scanType;
        _scan->productMz=productMz;
        _scan->collisionEnergy=collisionEnergy;

        if (peaks.attribute("byteOrder").empty() || strncasecmp(peaks.attribute("byteOrder").value(),"network",5) == 0 ) {
            networkorder = true;
        }

        int precision=32;
        if (!peaks.attribute("precision").empty()) { precision = peaks.attribute("precision").as_int(); }


        vector<float> mzint = base64::decode_base64(b64String,precision/8,networkorder,decompress);
        int size = mzint.size()/2;

        _scan->mz.resize(size);
        _scan->intensity.resize(size);

       //cerr << "new scan=" << scannum << " msL=" << msLevel << " rt=" << rt << " precMz=" << precursorMz << " polar=" << scanpolarity << " prec=" << precision << " byteorder=" << networkorder << " comprezed=" << decompress << "\t" << _scan->totalIntensity() << endl;

       // cerr << "Network:" << networkorder << " precision" << precision <<  " size=" << size << endl;

        int j=0; int count=0;
        for(int i=0; i < size; i++ ) {
            float mzValue = mzint[j++];
            float intensityValue = mzint[j++];
            if (mzValue > 0 && intensityValue > 0 ) {
				//cerr << mzValue << " " << intensityValue << endl;
                _scan->mz[i]= mzValue;
                _scan->intensity[i] = intensityValue;
                count++;
            }
        }


        _scan->mz.resize(count);
        _scan->intensity.resize(count);

        if (filterLine.empty() && precursorMz > 0 ) {
            _scan->filterLine = scanType + ":" + float2string(_scan->precursorMz,4) + " [" + float2string(_scan->productMz,4) + "]";
            cerr << " addScan:" << _scan->filterLine << endl;
        } else if (scanType  == "MRM") {
            _scan->filterLine = scanType + ":" + float2string(_scan->precursorMz,4) + " [" + float2string(_scan->productMz,4) + "]";
        }
    }
    return _scan;
}

void mzSample::openStream(const char* filename) {
    string filenameString = string(filename);
    this->sampleName = cleanSampleName(filename);
    this->fileName = filenameString;

    _iostream.open(filename);
}


Scan* mzSample::randomAccessMzXMLScan(int seek_pos_start, int seek_pos_end) {

    int length = seek_pos_end-seek_pos_start;

    if(length <= 1) {
        cerr << "Empty mzXML buffer" << endl;
        return NULL;
    }

    //read scan text
    char *buffer = new char[length+1];
    _iostream.seekg(seek_pos_start,_iostream.beg); _iostream.read(buffer,length);
    buffer[length]='\0';

    //parse scan text
    Scan* scan=NULL;
    try {
        xml_document doc;
        const unsigned int parse_options = parse_minimal;
        bool loadok = doc.load(buffer, parse_options);
        if (loadok )  {
            xml_node scanXMLData = doc.child("scan");
            scan = parseMzXMLScan(scanXMLData,0);
        } else {
            cerr << "Faild to parse buffer len=" << length << "\t" << buffer[0] << " " << buffer[length-1] << " ok="<< loadok << endl;
        }
    } catch (char* err) {
        cerr << "error: randomAccessMzXMLScan:" << err << endl;
    }

    delete[] buffer;
    return scan;
}

void mzSample::summary() { 
		cerr << "Num of obs:" << this->scans.size() << endl;
		cerr << "Rt range:" << this->minRt  << " " << this->maxRt << endl;
		cerr << "Mz range:" << this->minMz  << " " << this->maxMz << endl;
}

void mzSample::calculateMzRtRange() { 

		if (scans.size() == 0 ) {
			cerr << "sample has no data" << endl;
			return;
		}

		minRt = scans[0]->rt;
		maxRt = scans[scans.size()-1]->rt;
		minMz =  FLT_MAX;
		maxMz =  0;
		minIntensity = FLT_MAX;
		maxIntensity = 0;
		totalIntensity = 0;
		int nobs = 0;

		for (unsigned int j=0; j < scans.size(); j++ ) {
			for (unsigned int i=0; i < scans[j]->mz.size(); i++ ) {
				totalIntensity +=  scans[j]->intensity[i];
				float mz = scans[j]->mz[i]; 
				if( mz < minMz && mz > 0   ) minMz = mz; //sanity check must be greater > 0
				if (mz > maxMz && mz < 1e9 ) maxMz = mz; //sanity check m/z over a billion
				if (scans[j]->intensity[i] < minIntensity ) minIntensity = scans[j]->intensity[i];
				if (scans[j]->intensity[i] > maxIntensity ) maxIntensity = scans[j]->intensity[i];
				nobs++;
			}
		}
		//sanity check
		if (minRt <= 0 ) minRt = 0;
		if (maxRt >= 1e4 ) maxRt = 1e4;
        //cerr << "calculateMzRtRange() rt=" << minRt << "-" << maxRt << " mz=" << minMz << "-" << maxMz << endl;
}

mzSlice mzSample::getMinMaxDimentions(const vector<mzSample*>& samples) {
    mzSlice d;
    d.rtmin=0;
    d.rtmax=0;
    d.mzmin=0;
    d.mzmax=0;

	if ( samples.size() > 0 ) {

			d.rtmin = samples[0]->minRt;
			d.rtmax = samples[0]->maxRt;
			d.mzmin = samples[0]->minMz;
			d.mzmax = samples[0]->maxMz;

			for(unsigned int i=1; i < samples.size(); i++) {
					if ( samples[i]->minRt < d.rtmin ) d.rtmin=samples[i]->minRt;
					if ( samples[i]->maxRt > d.rtmax ) d.rtmax=samples[i]->maxRt;
					if ( samples[i]->minMz < d.mzmin ) d.mzmin=samples[i]->minMz;
					if ( samples[i]->maxMz > d.mzmax ) d.mzmax=samples[i]->maxMz;
			}
	}

	cerr << "getMinMaxDimentions() " << d.rtmin << " " << d.rtmax << " " << d.mzmin << " " << d.mzmax << endl;

	return d;
}

float mzSample::getMaxRt(const vector<mzSample*>&samples) { 
		float maxRt=0;
		for(unsigned int i=0; i < samples.size(); i++ ) 
				if (samples[i]->maxRt > maxRt) 
						maxRt=samples[i]->maxRt; 

		return maxRt;
}

float mzSample::getAverageFullScanTime() {

    //Issue 371: Lazy implementation, if already computed, return
    if (_averageFullScanTime >= 0) return _averageFullScanTime;

    float s=0;
	int n=0;
    Scan* lscan = nullptr;
    Scan* tscan = nullptr;
	if ( scans.size() == 0 ) return 0;

	for(unsigned int i=1; i < scans.size(); i++ ) {
		if ( scans[i]->mslevel == 1 ) {
				tscan = scans[i];
				if ( lscan ) { s += tscan->rt-lscan->rt; n++; }
				lscan = tscan;
		}
	}
    if ( n > 0 ){
        _averageFullScanTime = s/n;
    } else {
        _averageFullScanTime = 0.0f;
    }

    return _averageFullScanTime;
}



void mzSample::enumerateSRMScans() {

    srmScans.clear();
    srmScansByMzs.clear();

    for(unsigned int i=0; i < scans.size(); i++ ) {

        if (scans[i]->filterLine.length()>0) {

            string stringKey = scans[i]->filterLine;
            if(srmScans.find(stringKey) == srmScans.end()) {
                srmScans.insert(make_pair(stringKey, vector<int>()));
            }
            srmScans[stringKey].push_back(static_cast<int>(i));

            pair<float, float> mzKey = make_pair(scans[i]->precursorMz, scans[i]->productMz);
            if (srmScansByMzs.find(mzKey) == srmScansByMzs.end()) {
                srmScansByMzs.insert(make_pair(mzKey, vector<int>()));
            }
            srmScansByMzs[mzKey].push_back(static_cast<int>(i));

        }
    }
    cerr << "mzSample::enumerateSRMScans(): Found " << srmScans.size() << " SRM scans" << endl;
}

Scan* mzSample::getScan(unsigned int scanNum) {
	if ( scanNum >= scans.size() ) scanNum = scans.size()-1;
	if ( scanNum < scans.size() ) {
		return(scans[scanNum]);
	} else {
		cerr << "Warning bad scan number " << scanNum << endl;
        return nullptr;
	}
}

EIC* mzSample::getEIC(float precursorMz, float collisionEnergy, float productMz, float amuQ1=0.1, float amuQ2=0.5) {
    EIC* e = new EIC();
    e->sampleName = sampleName;
    e->sample = this;
    e->totalIntensity=0;
    e->maxIntensity = 0;
    e->mzmin = 0;
    e->mzmax = 0;

    for(unsigned int i=0; i < scans.size(); i++ ) {
        Scan* scan = scans[i];
        if (scan->mslevel < 2) continue;
        if (precursorMz && abs(scan->precursorMz-precursorMz)>amuQ1 ) continue;
        if (productMz && abs(scan->productMz-productMz)>amuQ2) continue;
        //if (collisionEnergy && abs(scan->collisionEnergy-collisionEnergy) > 0.5) continue;

        float maxMz=0;
        float maxIntensity=0;
        for(unsigned int k=0; k < scan->nobs(); k++ ) {
            if (scan->intensity[k] > maxIntensity ) {
                maxIntensity=scan->intensity[k];
                maxMz = scan->mz[k];
            }
        }

        e->scannum.push_back(scan->scannum);
        e->rt.push_back( scan->rt );
        e->intensity.push_back(maxIntensity);
        e->mz.push_back(maxMz);
        e->totalIntensity += maxIntensity;
        if (maxIntensity>e->maxIntensity) e->maxIntensity = maxIntensity;
    }

    if ( e->rt.size() > 0 ) {
        e->rtmin = e->rt[0];
        e->rtmax = e->rt[ e->size()-1];
    }

    float scale = getNormalizationConstant();
    if(scale != 1.0) for (unsigned int j=0; j < e->size(); j++) { e->intensity[j] *= scale; }

    //if(e->size() == 0) cerr << "getEIC(Q1,CE,Q3): is empty" << precursorMz << " " << collisionEnergy << " " << productMz << endl;
    //std::cerr << "getEIC(Q1,CE,Q3): srm"  << precursorMz << " " << e->intensity.size() << endl;
    return e;
}



EIC* mzSample::getEIC(string srm) {

        EIC* e = new EIC();
	e->sampleName = sampleName;
	e->sample = this;
        e->totalIntensity=0;
        e->maxIntensity = 0;
	e->mzmin = 0;
	e->mzmax = 0;

	if (srmScans.size() == 0 ) enumerateSRMScans();

	if (srmScans.count(srm) > 0 ) {
            vector<int> srmscans = srmScans[srm];
            for (unsigned int i=0; i < srmscans.size(); i++ ) {
                Scan* scan = scans[srmscans[i]];
                float maxMz=0;
                float maxIntensity=0;
                for(unsigned int k=0; k < scan->nobs(); k++ ) {
                    if (scan->intensity[k] > maxIntensity ) {
                        maxIntensity=scan->intensity[k];
                        maxMz = scan->mz[k];
                    }
                }
                e->scannum.push_back(scan->scannum);
                e->rt.push_back( scan->rt );
                e->intensity.push_back(maxIntensity);
                e->mz.push_back(maxMz);
                e->totalIntensity += maxIntensity;

                if (maxIntensity>e->maxIntensity) e->maxIntensity = maxIntensity;
            }
	}

	if ( e->rt.size() > 0 ) { 
                e->rtmin = e->rt[0];
                e->rtmax = e->rt[ e->size()-1];
	}

	float scale = getNormalizationConstant();
	if(scale != 1.0) for (unsigned int j=0; j < e->size(); j++) { e->intensity[j] *= scale; }
	if(e->size() == 0) cerr << "getEIC(SRM STRING): is empty" << srm << endl;
    //std::cerr << "getEIC: srm"  << srm << " " << e->intensity.size() << endl;

	return e;
}

//Issue 347
//<precursor ion m/z, product ion m/z> for SRM scans
EIC* mzSample::getEIC(pair<float, float> mzKey) {

    EIC* e = new EIC();
    e->sampleName = sampleName;
    e->sample = this;
    e->totalIntensity=0;
    e->maxIntensity = 0;
    e->mzmin = 0;
    e->mzmax = 0;

    if (srmScansByMzs.size() == 0 ) enumerateSRMScans();

    if (srmScansByMzs.count(mzKey) > 0 ) {
            vector<int> srmscans = srmScansByMzs[mzKey];
            for (unsigned int i=0; i < srmscans.size(); i++ ) {
                Scan* scan = scans[static_cast<unsigned long>(srmscans[i])];
                float maxMz=0;
                float maxIntensity=0;
                for(unsigned int k=0; k < scan->nobs(); k++ ) {
                    if (scan->intensity[k] > maxIntensity ) {
                        maxIntensity=scan->intensity[k];
                        maxMz = scan->mz[k];
                    }
                }
                e->scannum.push_back(scan->scannum);
                e->rt.push_back( scan->rt );
                e->intensity.push_back(maxIntensity);
                e->mz.push_back(maxMz);
                e->totalIntensity += maxIntensity;

                if (maxIntensity>e->maxIntensity) e->maxIntensity = maxIntensity;
            }
    }

    if ( e->rt.size() > 0 ) {
                e->rtmin = e->rt[0];
                e->rtmax = e->rt[ e->size()-1];
    }

    float scale = getNormalizationConstant();
    if(scale != 1.0) for (unsigned int j=0; j < e->size(); j++) { e->intensity[j] *= scale; }
    if(e->size() == 0) cerr << "getEIC(pair<float, float> mzKey): is empty key=(" << mzKey.first << ", " << mzKey.second << ")" << endl;
    //std::cerr << "getEIC: srm"  << srm << " " << e->intensity.size() << endl;

return e;

}


EIC* mzSample::getEIC(float mzmin, float mzmax, float rtmin, float rtmax, int mslevel, string scanFilterString) {

    //ajust EIC retention time window to match sample retentention times
	if (rtmin < this->minRt ) rtmin = this->minRt;
	if (rtmax > this->maxRt ) rtmax = this->maxRt;
	if (mzmin < this->minMz ) mzmin = this->minMz;
	if (mzmax > this->maxMz ) mzmax = this->maxMz;

   //cerr << "getEIC()" << setprecision(7) << mzmin << " " << mzmax << " " << rtmin << " " << rtmax << endl;

	EIC* e = new EIC();
	e->sampleName = sampleName;
	e->sample = this;
	e->mzmin = mzmin;
	e->mzmax = mzmax;
	e->totalIntensity=0;
	e->maxIntensity = 0;

	int scanCount = scans.size();
	if (scanCount == 0) return e;

	if ( mzmin < minMz && mzmax < maxMz ) {
			cerr << "getEIC(): mzmin and mzmax are out of range" << endl;
			return e;
	}
    
	//binary search mz  domain iterator
	vector<float>::iterator mzItr;

        //binary search rt domain iterator
        Scan tmpScan(this,0,1,rtmin-0.1,0,-1);
        deque<Scan*>::iterator scanItr = lower_bound(scans.begin(), scans.end(),&tmpScan, Scan::compRt);
        if (scanItr >= scans.end() ) { return e; }

        //preallocated memory for arrays [ this should really be corrected for mslevel type ]
        int estimatedScans=scans.size();

        if (this->maxRt-this->minRt > 0 && (rtmax-rtmin)/(this->maxRt-this->minRt) <= 1 ) {
            estimatedScans=float (rtmax-rtmin)/(this->maxRt-this->minRt)*scans.size();
        }

	if (estimatedScans < 512 ) estimatedScans=512;
	if (estimatedScans > scans.size() ) estimatedScans=scans.size();

	e->scannum.reserve(estimatedScans);
	e->rt.reserve(estimatedScans);
	e->intensity.reserve(estimatedScans);
	e->mz.reserve(estimatedScans);

	int scanNum=scanItr-scans.begin()-1;
	for(; scanItr != scans.end(); ++scanItr) {
            Scan* scan = *(scanItr);

            scanNum++;
            if (scan->mslevel != mslevel) continue;
            if (scan->rt < rtmin) continue;
            if (scan->mz.size() == 0) continue;
            if (scan->rt > rtmax) break;

            //sim scan outside the range
            if (mslevel == 1 && (scan->mz.front() > mzmax || scan->mz.back() < mzmin)) continue;

            //Issue 222: respect filter string
            if (!scanFilterString.empty() && scan->filterString.find(scanFilterString) == string::npos) continue;

            float __maxMz=0;
            float __maxIntensity=0;

            //binary search
            mzItr = lower_bound(scan->mz.begin(), scan->mz.end(), mzmin);
            int lb = mzItr-scan->mz.begin();

            for(unsigned int k=lb; k < scan->nobs(); k++ ) {
                if (scan->mz[k] < mzmin) continue;
                if (scan->mz[k] > mzmax) break;
                if (scan->intensity[k] > __maxIntensity ) {
                    __maxIntensity=scan->intensity[k];
                    __maxMz = scan->mz[k];
                }
            }

            e->scannum.push_back(scanNum);
            e->rt.push_back(scan->rt);
            e->intensity.push_back(__maxIntensity);
            e->mz.push_back(__maxMz);
            e->totalIntensity += __maxIntensity;
            if (__maxIntensity>e->maxIntensity) e->maxIntensity = __maxIntensity;

	}

    //cerr << "estimatedScans=" << estimatedScans << " actul=" << e->scannum.size() << endl;
	if ( e->rt.size() > 0 ) { 
          e->rtmin = e->rt[0];
          e->rtmax = e->rt[ e->size()-1];
         //cerr << "getEIC()" << e->scannum[0] << " " << e->scannum[e->scannum.size()-1] << " " << scans.size() << endl;
	}

	//scale EIC by normalization constant
	float scale = getNormalizationConstant();
	if(scale != 1.0) for (unsigned int j=0; j < e->size(); j++) { e->intensity[j] *= scale; }
	
	//cerr << "getEIC: maxIntensity=" << e->maxIntensity << endl;
	//e->summary();

	if(e->size() == 0) cerr << "getEIC(mzrange,rtrange,mslevel): is empty: " << mzmin << " " << mzmax << " " << rtmin << " " << rtmax << endl;

	return(e);
}


EIC* mzSample::getTIC(float rtmin, float rtmax, int mslevel) { 

    //ajust EIC retention time window to match sample retentention times
    if (rtmin < this->minRt ) rtmin =this->minRt;
    if (rtmax > this->maxRt ) rtmax =this->maxRt;

    //cerr << "getEIC()" << setprecision(7) << mzmin << " " << mzmax << " " << rtmin << " " << rtmax << endl;

    EIC* e = new EIC();
    e->sampleName = sampleName;
    e->sample = this;
    e->mzmin = 0;
    e->mzmax = 0;
    e->totalIntensity=0;
    e->maxIntensity = 0;

    int scanCount = scans.size();
    if ( scanCount == 0 ) return e;

    for(int i=0;  i < scanCount; i++)
    {
        if (scans[i]->mslevel == mslevel) {
            Scan* scan = scans[i];
            float y = scan->totalIntensity();
            e->mz.push_back(0);
            e->scannum.push_back(i);
            e->rt.push_back(scan->rt);
            e->intensity.push_back(y);
            e->totalIntensity += y;
            if (y>e->maxIntensity) e->maxIntensity=y;
        }
    }
    if ( e->rt.size() > 0 ) {
        e->rtmin = e->rt[0];
        e->rtmax = e->rt[ e->size()-1];
    }
    return(e);
}

EIC* mzSample::getBIC(float rtmin, float rtmax, int mslevel) { 

    //ajust EIC retention time window to match sample retentention times
    if (rtmin < this->minRt ) rtmin =this->minRt;
    if (rtmax > this->maxRt ) rtmax =this->maxRt;

    //cerr << "getEIC()" << setprecision(7) << mzmin << " " << mzmax << " " << rtmin << " " << rtmax << endl;

    EIC* e = new EIC();
    e->sampleName = sampleName;
    e->sample = this;
    e->mzmin = 0;
    e->mzmax = 0;
    e->totalIntensity=0;
    e->maxIntensity = 0;

    int scanCount = scans.size();
    if ( scanCount == 0 ) return e;

    for(int i=0;  i < scanCount; i++)
    {
        if (scans[i]->mslevel == mslevel) {
            Scan* scan = scans[i];
			float maxMz=0;
			float maxIntensity=0;
    		for(unsigned int i=0;i<scan->intensity.size();i++)  {
					if(scan->intensity[i] >maxIntensity) { maxIntensity=scan->intensity[i]; maxMz=scan->mz[i] ; }
			}
            e->mz.push_back(maxMz);
            e->scannum.push_back(i);
            e->rt.push_back(scan->rt);
            e->intensity.push_back(maxIntensity);
            e->totalIntensity += maxIntensity;
            if (maxIntensity>e->maxIntensity) e->maxIntensity=maxIntensity;
        }
    }
    if ( e->rt.size() > 0 ) {
        e->rtmin = e->rt[0];
        e->rtmax = e->rt[ e->size()-1];
    }
    return(e);
}


//compute correlation between two mzs within some retention time window
float mzSample::correlation(float mz1,  float mz2, float ppm, float rt1, float rt2, bool debug) {
    
    float ppm1 = ppm*mz1/1e6;
    float ppm2 = ppm*mz2/1e6;
    int mslevel=1;
    EIC* e1 = mzSample::getEIC(mz1-ppm1, mz1+ppm1, rt1, rt2, mslevel);
    EIC* e2 = mzSample::getEIC(mz2-ppm2, mz2+ppm1, rt1, rt2, mslevel);

    double correlation = mzUtils::correlation(e1->intensity, e2->intensity);

    if (debug) {

        //note that if the vector is all 0s, this will still compute the correlation!
        cout << "mzSample::correlation() mz1=" << to_string(mz1) << ", mz2=" <<to_string(mz2) << endl;
        cout << "first: " << e1->intensity.size() <<", second: " << e2->intensity.size() << endl;

        for (unsigned int i = 0; i < e1->intensity.size(); i++) {
            cout << e1->intensity.at(i) << "\t" << e2->intensity.at(i) << endl;
        }

        cout << "mz_1 <- c(";
        for (unsigned int i = 0; i < e1->intensity.size(); i++) {
            if (i>0) {
                cout << ", ";
            }
            cout << e1->intensity.at(i);
        }
        cout << ")" << endl;

        cout << "mz_2 <- c(";
        for (unsigned int i = 0; i < e2->intensity.size(); i++) {
            if (i>0) {
                cout << ", ";
            }
            cout << e2->intensity.at(i);
        }
        cout << ")" << endl;

        cout << "correlation= " << correlation << endl;

        cout << "===========" << endl;
    }

    delete(e1);
    delete(e2);
    return(correlation);
}


int mzSample::parseCDF (char* filename, int is_verbose)
{
    #ifdef CDFPARSER
    int cdf=0;
    int errflag = 0;
    long nscans=0;
    long ninst=0;

    extern int ncopts;              /* from "netcdf.h" */
    ncopts = 0;

    static MS_Admin_Data            admin_data;
    static MS_Sample_Data           sample_data;
    static MS_Test_Data             test_data;
    static MS_Instrument_Data       inst_data;
    static MS_Raw_Data_Global       raw_global_data;
    static MS_Raw_Per_Scan          raw_data;
    double  mass_pt=0;
    double inty_pt=0;
    double inty=0;

    cdf = ms_open_read( filename );
    if ( -1 == cdf )
    {
        fprintf( stderr, "\nopen_cdf_ms: ms_open_read failed!" );
        return 0;
    }

    /* Initialize attribute data structures */

    ms_init_global( FALSE, &admin_data, &sample_data, &test_data, &raw_global_data );


    /* Read global information */

    if (MS_ERROR == ms_read_global( cdf, &admin_data, &sample_data, &test_data, &raw_global_data))
    {
        fprintf( stderr, "\nopen_cdf_ms: ms_read_global failed!" );
        ms_init_global( TRUE, &admin_data, &sample_data, &test_data, &raw_global_data);
        ms_close(cdf);
        return 0;
    }

    nscans = raw_global_data.nscans;

    switch (admin_data.experiment_type)
    {
    case 0:
        printf ("Centroid");
        break;
    case 1:
        printf ("Continuum");
        break;
    case 2:
        printf ("Library");
        break;
    default:
        printf ("Unknown: '%d'", admin_data.experiment_type);
        break;
    }

    printf ("\n\n-- Instrument Information --");
    ninst = admin_data.number_instrument_components;
    printf ("\nNumber_inst_comp\t%ld", ninst);


    printf ("\n\n-- Raw Data Information --");
    printf ("\nNumber of scans\t\t%ld", nscans);

    printf ("\nMass Range\t\t%.2f > %.2f",
            raw_global_data.mass_axis_global_min,
            raw_global_data.mass_axis_global_max);

    printf ("\nInty Range\t\t%.2f > %.2f",
            raw_global_data.intensity_axis_global_min,
            raw_global_data.intensity_axis_global_max);

    printf ("\nTime Range\t\t%.2f > %.2f",
            raw_global_data.time_axis_global_min,
            raw_global_data.time_axis_global_max);

    printf ("\nActual Run Time\t\t%.2f (%.2f min)",
            raw_global_data.run_time, raw_global_data.run_time/60.0);
    printf ("\nComments \t\t%s", raw_global_data.comments);


    if (errflag)                    /* if error occurred, clean up and leave */
    {
        ms_init_global( TRUE, &admin_data, &sample_data, &test_data, &raw_global_data );
        ms_close( cdf );
        return 0;
    }

    /* Check to see if scale factors and offsets are set to "NULL"
        values; if so, correct them for use below */

    if ((int)MS_NULL_FLT == (int)raw_global_data.mass_factor)
        raw_global_data.mass_factor = 1.0;

    if ((int)MS_NULL_FLT == (int)raw_global_data.time_factor)
        raw_global_data.time_factor = 1.0;

    if ((int)MS_NULL_FLT == (int)raw_global_data.intensity_factor)
        raw_global_data.intensity_factor = 1.0;

    if ((int)MS_NULL_FLT == (int)raw_global_data.intensity_offset)
        raw_global_data.intensity_offset = 0.0;

    if ((raw_global_data.mass_axis_global_min < 0) || (raw_global_data.mass_axis_global_max < 0))
    {
        /* this bug is frequently observed with files from HP/Agilent ChemStation */
        fprintf (stderr, "\n*** WARNING: Negative mass reported! Use '-v' for details.\n\n");
    }

    for (int scan = 0; scan < nscans; scan++)
    {
        ms_init_per_scan(FALSE, &raw_data, NULL);
        raw_data.scan_no = (long) scan;
        mass_pt = inty_pt = inty = 0.0;                     /* init */

        if (MS_ERROR == ms_read_per_scan(cdf, &raw_data, NULL))
        {               /* free allocated memory before leaving */
            fprintf(stderr, "\nreadchro: ms_read_per_scan failed (scan %d)!", scan);
            ms_init_per_scan(TRUE, &raw_data, NULL);
            return 0;
        }

        if (!raw_data.points) {       /* empty scan? */
            break;
        } else {                       /* there are data points */

            int polarity=0;
            if ( test_data.ionization_mode == polarity_plus) polarity=+1;
            else polarity = -1;


            Scan* myscan = new Scan(this,raw_data.actual_scan_no,
                                    test_data.scan_function - (int) resolution_proportional,
                                    raw_data.scan_acq_time/60,
                                    0,
                                    polarity);


            myscan->intensity.resize(raw_data.points);
            myscan->mz.resize(raw_data.points);

            if (admin_data.experiment_type == 0)
                myscan->centroided=true;
            else
                myscan->centroided=false;

            for (int i = 0; i < raw_data.points; i++)
            {
                switch( raw_global_data.mass_format )
                {
                case data_short:
                    mass_pt = (double) ((short *)raw_data.masses)[i];
                    break;

                case data_long:
                    mass_pt = (double) ((long *)raw_data.masses)[i];
                    break;

                case data_float:
                    mass_pt = (double) ((float *)raw_data.masses)[i];
                    break;

                case data_double:
                    mass_pt = ((double *)raw_data.masses)[i];
                    break;
                }

                mass_pt *= raw_global_data.mass_factor;

                switch( raw_global_data.intensity_format )
                {
                case data_short:
                    inty_pt = (double) ((short *)raw_data.intensities)[i];
                    break;

                case data_long:
                    inty_pt = (double) ((long *)raw_data.intensities)[i];
                    break;

                case data_float:
                    inty_pt = (double) ((float *)raw_data.intensities)[i];
                    break;

                case data_double:
                    inty_pt = (double) ((double *)raw_data.intensities)[i];
                    break;
                }

                inty_pt = inty_pt * raw_global_data.intensity_factor + raw_global_data.intensity_offset;
                //cerr << "mz/int" << mass_pt << " " << inty_pt << endl;
                myscan->intensity[i]=inty_pt;
                myscan->mz[i]=mass_pt;

               // if (raw_data.flags > 0) printf("\nWarning: There are flags in scan %ld (ignored).", scan);
            }       /* i loop */


            addScan(myscan);
        }


        ms_init_per_scan( TRUE, &raw_data, NULL );

    }   /* scan loop */
    #endif
    return 1;
}

Scan* mzSample::getAverageScan(float rtmin, float rtmax, int mslevel, int polarity, float sd ) {

    float rt=rtmin + (rtmax-rtmin)/2;
    int scanCount=0;
    int scannum=0;

    map<float,double> mz_intensity_map;
    map<float,double> mz_bin_map;
    map<float,int> mz_count;

    for(unsigned int s=0; s < scans.size(); s++) {
        if(scans[s]->getPolarity() != polarity
                || scans[s]->mslevel != mslevel
                || scans[s]->rt < rtmin
                || scans[s]->rt > rtmax) continue;

        Scan* scan = scans[s];
        scanCount++;
        for(unsigned int i=0; i < scan->mz.size(); i++) {
                float bin = FLOATROUND(scan->mz[i],sd);
                mz_intensity_map[bin] += ((double) scan->intensity[i]);
                mz_bin_map[bin] += ((double)(scan->intensity[i])*(scan->mz[i]));
                mz_count[bin]++;
        }
    }

    Scan* avgScan = new Scan(this,scannum,mslevel,rt/scanCount, 0, polarity);

    map<float,double>::iterator itr;
    for(itr = mz_intensity_map.begin(); itr != mz_intensity_map.end(); ++itr ) {
        float bin = (*itr).first;
        double totalIntensity=(*itr).second;
        double avgMz =  mz_bin_map[bin] / totalIntensity;
        avgScan->mz.push_back((float)avgMz);
        avgScan->intensity.push_back((float) totalIntensity / mz_count[bin]);
    }
    //cout << "getAverageScan() from:" << from << " to:" << to << " scanCount:" << scanCount << "scans. mzs=" << avgScan->nobs() << endl;
    return avgScan;
}


void mzSample::saveOriginalRetentionTimes() {
    if ( originalRetentionTimes.size() > 0 ) return;

    originalRetentionTimes=vector<float>(scans.size(),0);
    for(unsigned int ii=0; ii < scans.size(); ii++ ) {
        originalRetentionTimes[ii]=scans[ii]->rt;
    }
}

void mzSample::restoreOriginalRetentionTimes() {
    if ( originalRetentionTimes.size() == 0 ) return;

    for(unsigned int ii=0; ii < scans.size(); ii++ ) {
        scans[ii]->rt = originalRetentionTimes[ii];
    }
}


vector<float> mzSample::getIntensityDistribution(int mslevel) {

    vector<float> allintensities;
    for(unsigned int s=0; s < this->scans.size(); s++) {
        Scan* scan = this->scans[s];
        if (scan->mslevel != mslevel) continue;

        for(unsigned int i=0; i < scan->mz.size(); i++) {
            allintensities.push_back(scan->intensity[i]);
        }
    }

    return(quantileDistribution(allintensities));
}


void mzSample::applyPolynomialTransform() { 
	int poly_align_degree = polynomialAlignmentTransformation.size()-1;
	if (poly_align_degree <= 0) return;

	double* transform = &polynomialAlignmentTransformation.front();
	for(int i=0; i<scans.size(); i++ ) {
		float newrt = leasev(transform, poly_align_degree, scans[i]->rt);
		//cerr << "applyPolynomialTransform() " << scans[i]->rt << "\t" << newrt << endl;
		scans[i]->rt = newrt;
	}
}

double mzSample::getMS1PrecursorMass(Scan* ms2scan,float ppm) {

    if (ms2scan->precursorMz == 0 ) return 0;
    int scanNum = ms2scan->scannum;

    double adjPreMass=0;
    for(int i=scanNum; i>scanNum-50; i--) {
        Scan* scan = this->getScan(i);

        if (!scan or scan->mslevel > 1) continue;

        for( float mult=1; mult <= 5; mult++) {
            //find highest intensity posiion is last ms1 scan
            int pos = scan->findHighestIntensityPos(ms2scan->precursorMz,ppm*mult);
            if (pos >0) {
                adjPreMass= scan->mz[pos];
                break;
            }
        }

       //out << "HIT: " << " scan=" << scan->scannum << "  ms2=" << setprecision(10) << ms2scan->precursorMz << " ms1=" << scan->mz[pos] << endl;
        break;
    }

    if ( adjPreMass > 0 ) {
        //cerr << "adjPreMass : ms2p=" << ms2scan->precursorMz << " adjP" << adjPreMass << "  err="  << mzUtils::ppmDist((double) ms2scan->precursorMz,adjPreMass) << endl;
        return adjPreMass;
    } else {
        //cerr << "getMS1PrecurursorMass() CAN'T FIND PRECURSOR " << ms2scan->precursorMz << endl;
        return ms2scan->precursorMz;
    }
}

vector<Scan*> mzSample::getFragmentationEvents(mzSlice* slice) {
    vector<Scan*>matchedscans;
    for( unsigned int j=0; j < scans.size(); j++ ) {
            Scan* scan = scans[j];
            if (!scan or scan->mslevel <= 1 or scan->rt < slice->rtmin) continue; //skip ms1 events
            if (scan->rt > slice->rtmax) break;
            if (scan->precursorMz >= slice->mzmin and scan->precursorMz <= slice->mzmax) {
                matchedscans.push_back(scan);
            }
     }
    return matchedscans;
}

vector<Adduct*> Adduct::loadAdducts(string filename){
    vector<Adduct*> adducts;
    ifstream myfile(filename.c_str());

    if (! myfile.is_open()) return adducts;

    string line;
    while ( getline(myfile,line) ) {
        if (line.empty()) continue;
        if (line[0] == '#') continue;

        vector<string>fields;
        mzUtils::split(line,',', fields);

        //ionization
        if(fields.size() < 2 ) continue;
        string name=fields[0];
        int nmol=string2float(fields[1]);
        int charge=string2float(fields[2]);
        float mass=string2float(fields[3]);

        if ( name.empty() || nmol < 0 ) continue;
        Adduct* a = new Adduct(name,mass,charge,nmol);
        adducts.push_back(a);
    }
    myfile.close();
    return adducts;
}

string PeaksSearchParameters::encodeParams(){

    string encodedParams;

    //START SearchParameters

    //program level
    encodedParams = encodedParams + "searchVersion" + "=" + searchVersion + ";";

    //scan filter params
    encodedParams = encodedParams + "scanFilterMinFracIntensity" + "=" + to_string(scanFilterMinFracIntensity) + ";";
    encodedParams = encodedParams + "scanFilterMinSNRatio" + "=" + to_string(scanFilterMinSNRatio) + ";";
    encodedParams = encodedParams + "scanFilterMaxNumberOfFragments" + "=" + to_string(scanFilterMaxNumberOfFragments) + ";";
    encodedParams = encodedParams + "scanFilterBaseLinePercentile" + "=" + to_string(scanFilterBaseLinePercentile) + ";";
    encodedParams = encodedParams + "scanFilterIsRetainFragmentsAbovePrecursorMz" + "=" + to_string(scanFilterIsRetainFragmentsAbovePrecursorMz) + ";";
    encodedParams = encodedParams + "scanFilterPrecursorPurityPpm" + "=" + to_string(scanFilterPrecursorPurityPpm) + ";";
    encodedParams = encodedParams + "scanFilterMinIntensity" + "=" + to_string(scanFilterMinIntensity) + ";";

    //scan filter for MS1 scans
    encodedParams = encodedParams + "scanFilterMs1MinRt" + "=" + to_string(scanFilterMs1MinRt) + ";";
    encodedParams = encodedParams + "scanFilterMs1MaxRt" + "=" + to_string(scanFilterMs1MaxRt) + ";";

    //scan filter for MS2 scans
    encodedParams = encodedParams + "scanFilterMs2MinRt" + "=" + to_string(scanFilterMs2MinRt) + ";";
    encodedParams = encodedParams + "scanFilterMs2MaxRt" + "=" + to_string(scanFilterMs2MaxRt) + ";";

    //consensus spectrum params
    encodedParams = encodedParams + "consensusIsIntensityAvgByObserved" + "=" + to_string(consensusIsIntensityAvgByObserved) + ";";
    encodedParams = encodedParams + "consensusIsNormalizeTo10K" + "=" + to_string(consensusIsNormalizeTo10K) + ";";
    string consensusIntensityAgglomerationTypeStr = "UNSPECIFIED";
    if (consensusIntensityAgglomerationType == Fragment::ConsensusIntensityAgglomerationType::Mean) {
        consensusIntensityAgglomerationTypeStr = "MEAN";
    } else if (consensusIntensityAgglomerationType == Fragment::ConsensusIntensityAgglomerationType::Median) {
        consensusIntensityAgglomerationTypeStr = "MEDIAN";
    }
    encodedParams = encodedParams + "consensusIntensityAgglomerationType" + "=" + consensusIntensityAgglomerationTypeStr + ";";

    //ms1 consensus spectrum params
    encodedParams = encodedParams + "consensusMs1PpmTolr" + "=" + to_string(consensusMs1PpmTolr) + ";";
    encodedParams = encodedParams + "consensusMinNumMs1Scans" + "=" + to_string(consensusMinNumMs1Scans) + ";";
    encodedParams = encodedParams + "consensusMinFractionMs1Scans" + "=" + to_string(consensusMinFractionMs1Scans) + ";";

    //ms2 consensus spectrum params
    encodedParams = encodedParams + "consensusPpmTolr" + "=" + to_string(consensusPpmTolr) + ";";
    encodedParams = encodedParams + "consensusMinNumMs2Scans" + "=" + to_string(consensusMinNumMs2Scans) + ";";
    encodedParams = encodedParams + "consensusMinFractionMs2Scans" + "=" + to_string(consensusMinFractionMs2Scans) + ";";

    // ms1 matching
    encodedParams = encodedParams + "ms1PpmTolr" + "=" + to_string(ms1PpmTolr) + ";";

    // ms2 search
    encodedParams = encodedParams + "ms2MinNumMatches" + "=" + to_string(ms2MinNumMatches) + ";";
    encodedParams = encodedParams + "ms2MinNumDiagnosticMatches" + "=" + to_string(ms2MinNumDiagnosticMatches) + ";";
    encodedParams = encodedParams + "ms2MinNumUniqueMatches" + "=" + to_string(ms2MinNumUniqueMatches) + ";";
    encodedParams = encodedParams + "ms2PpmTolr" + "=" + to_string(ms2PpmTolr) + ";";
    encodedParams = encodedParams + "ms2MinIntensity" + "=" + to_string(ms2MinIntensity) + ";";

    // END SearchParameters
    // START PeaksSearchParameters

    //baseline
    encodedParams = encodedParams + "baselineSmoothingWindow" + "=" + to_string(baselineSmoothingWindow) + ";";
    encodedParams = encodedParams + "baselineDropTopX" + "=" + to_string(baselineDropTopX) + ";";

    //eic
    encodedParams = encodedParams + "eicSmoothingWindow" + "=" + to_string(eicSmoothingWindow) + ";";
    encodedParams = encodedParams + "eicEicSmoothingAlgorithm" + "=" + eicEicSmoothingAlgorithm + ";";
    encodedParams = encodedParams + "eicMaxPeakGroupRtDiff" + "=" + to_string(eicMaxPeakGroupRtDiff) + ";";
    encodedParams = encodedParams + "eicNoNoiseObs" + "=" + to_string(eicNoNoiseObs) + ";";

    //quality
    encodedParams = encodedParams + "qualitySignalBaselineRatio" + "=" + to_string(qualitySignalBaselineRatio) + ";";
    encodedParams = encodedParams + "qualitySignalBlankRatio" + "=" + to_string(qualitySignalBlankRatio) + ";";
    encodedParams = encodedParams + "qualityMinPeakGroupIntensity" + "=" + to_string(qualityMinPeakGroupIntensity) + ";";
    encodedParams = encodedParams + "qualityMinPeakWidth" + "=" + to_string(qualityMinPeakWidth) + ";";
    encodedParams = encodedParams + "qualityMinGoodPeakPerGroup" + "=" + to_string(qualityMinGoodPeakPerGroup) + ";";
    encodedParams = encodedParams + "qualityMinPeakQuality" + "=" + to_string(qualityMinPeakQuality) + ";";
    encodedParams = encodedParams + "qualityClassifierModelName" + "=" + qualityClassifierModelName + ";";

    //isotopes
    encodedParams = encodedParams + "isotopesIsRequireMonoisotopicPeaks" + "=" + to_string(isotopesIsRequireMonoisotopicPeaks) + ";";
    encodedParams = encodedParams + "isotopesExtractIsotopicPeaks" + "=" + to_string(isotopesExtractIsotopicPeaks) + ";";
    encodedParams = encodedParams + "isotopesMzTolerance" + "=" + to_string(isotopesMzTolerance) + ";";

    //ms1
    encodedParams = encodedParams + "ms1IsMatchRtFlag" + "=" + to_string(ms1IsMatchRtFlag) + ";";
    encodedParams = encodedParams + "ms1RtTolr" + "=" + to_string(ms1RtTolr) + ";";
    encodedParams = encodedParams + "ms1MassSliceMergeNumScans" + "=" + to_string(ms1MassSliceMergeNumScans) + ";";
    encodedParams = encodedParams + "ms1MassSliceMergePpm" + "=" + to_string(ms1MassSliceMergePpm) + ";";

    //ms2 matching
    encodedParams = encodedParams + "ms2IsMatchMs2" + "=" + to_string(ms2IsMatchMs2) + ";";
    encodedParams = encodedParams + "ms2ScoringAlgorithm" + "=" + ms2ScoringAlgorithm + ";";
    encodedParams = encodedParams + "ms2MinScore" + "=" + to_string(ms2MinScore) + ";";

    //matching options
    encodedParams = encodedParams + "matchingIsRequireAdductPrecursorMatch" + "=" + to_string(matchingIsRequireAdductPrecursorMatch) + ";";
    encodedParams = encodedParams + "matchingIsRetainUnknowns" + "=" + to_string(matchingIsRetainUnknowns) + ";";
    encodedParams = encodedParams + "matchingLibraries" + "=" + matchingLibraries + ";";

    string matchingPolicyStr;
    if (matchingPolicy == PeakGroupCompoundMatchingPolicy::SINGLE_TOP_HIT) {
        matchingPolicyStr = "SINGLE_TOP_HIT";
    } else if (matchingPolicy == PeakGroupCompoundMatchingPolicy::ALL_MATCHES) {
        matchingPolicyStr = "ALL_MATCHES";
    } else if (matchingPolicy == PeakGroupCompoundMatchingPolicy::TOP_SCORE_HITS) {
        matchingPolicyStr = "TOP_SCORE_HITS";
    }
    encodedParams = encodedParams + "matchingPolicy" + "=" + matchingPolicyStr + ";";

    // END PeaksSearchParameters

    return encodedParams;
}

shared_ptr<PeaksSearchParameters> PeaksSearchParameters::decode(string encodedParams){
    shared_ptr<PeaksSearchParameters> peaksSearchParameters = shared_ptr<PeaksSearchParameters>(new PeaksSearchParameters());

    unordered_map<string, string> decodedMap = mzUtils::decodeParameterMap(encodedParams); //use semicolon (default)

    //START SearchParameters

    //program level
    if (decodedMap.find("searchVersion") != decodedMap.end()) {
        peaksSearchParameters->searchVersion = decodedMap["searchVersion"];
    }

    //scan filter params
    if (decodedMap.find("scanFilterMinFracIntensity") != decodedMap.end()){
        peaksSearchParameters->scanFilterMinFracIntensity = stof(decodedMap["scanFilterMinFracIntensity"]);
    }
    if (decodedMap.find("scanFilterMinSNRatio") != decodedMap.end()){
        peaksSearchParameters->scanFilterMinSNRatio = stof(decodedMap["scanFilterMinSNRatio"]);
    }
    if (decodedMap.find("scanFilterMaxNumberOfFragments") != decodedMap.end()) {
        peaksSearchParameters->scanFilterMaxNumberOfFragments = stoi(decodedMap["scanFilterMaxNumberOfFragments"]);
    }
    if (decodedMap.find("scanFilterBaseLinePercentile") != decodedMap.end()) {
        peaksSearchParameters->scanFilterBaseLinePercentile = stoi(decodedMap["scanFilterBaseLinePercentile"]);
    }
    if (decodedMap.find("scanFilterIsRetainFragmentsAbovePrecursorMz") != decodedMap.end()) {
        peaksSearchParameters->scanFilterIsRetainFragmentsAbovePrecursorMz = decodedMap["scanFilterIsRetainFragmentsAbovePrecursorMz"] == "1";
    }
    if (decodedMap.find("scanFilterPrecursorPurityPpm") != decodedMap.end()){
        peaksSearchParameters->scanFilterPrecursorPurityPpm = stof(decodedMap["scanFilterPrecursorPurityPpm"]);
    }
    if (decodedMap.find("scanFilterMinIntensity") != decodedMap.end()){
        peaksSearchParameters->scanFilterMinIntensity = stof(decodedMap["scanFilterMinIntensity"]);
    }

    //scan filter for MS1 scans
    if (decodedMap.find("scanFilterMs1MinRt") != decodedMap.end()) {
        peaksSearchParameters->scanFilterMs1MinRt = stof(decodedMap["scanFilterMs1MinRt"]);
    }
    if (decodedMap.find("scanFilterMs1MaxRt") != decodedMap.end()) {
        peaksSearchParameters->scanFilterMs1MaxRt = stof(decodedMap["scanFilterMs1MaxRt"]);
    }

    //scan filter for MS2 scans
    if (decodedMap.find("scanFilterMs2MinRt") != decodedMap.end()) {
        peaksSearchParameters->scanFilterMs2MinRt = stof(decodedMap["scanFilterMs2MinRt"]);
    }
    if (decodedMap.find("scanFilterMs2MaxRt") != decodedMap.end()) {
        peaksSearchParameters->scanFilterMs2MaxRt = stof(decodedMap["scanFilterMs2MaxRt"]);
    }

    //consensus spectrum params (all ms levels)

    if (decodedMap.find("consensusIsIntensityAvgByObserved") != decodedMap.end()){
        peaksSearchParameters->consensusIsIntensityAvgByObserved = decodedMap["consensusIsIntensityAvgByObserved"] == "1";
    }
    if (decodedMap.find("consensusIsNormalizeTo10K") != decodedMap.end()){
        peaksSearchParameters->consensusIsNormalizeTo10K = decodedMap["consensusIsNormalizeTo10K"] == "1";
    }
    if (decodedMap.find("consensusIntensityAgglomerationType") != decodedMap.end()) {
        string consensusIntensityAgglomerationTypeStr = decodedMap["consensusIntensityAgglomerationType"];
        if (consensusIntensityAgglomerationTypeStr == "MEAN") {
            peaksSearchParameters->consensusIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Mean;
        } else if (consensusIntensityAgglomerationTypeStr == "MEDIAN") {
            peaksSearchParameters->consensusIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Median;
        }
    }

    //ms1 consensus spectrum params
    if (decodedMap.find("consensusMs1PpmTolr") != decodedMap.end()){
        peaksSearchParameters->consensusMs1PpmTolr = stof(decodedMap["consensusMs1PpmTolr"]);
    }
    if (decodedMap.find("consensusMinNumMs1Scans") != decodedMap.end()){
        peaksSearchParameters->consensusMinNumMs1Scans = stoi(decodedMap["consensusMinNumMs1Scans"]);
    }
    if (decodedMap.find("consensusMinFractionMs1Scans") != decodedMap.end()){
        peaksSearchParameters->consensusMinFractionMs1Scans = stof(decodedMap["consensusMinFractionMs1Scans"]);
    }

    //ms2 consensus spectrum params
    if (decodedMap.find("consensusPpmTolr") != decodedMap.end()){
        peaksSearchParameters->consensusPpmTolr = stof(decodedMap["consensusPpmTolr"]);
    }
    if (decodedMap.find("consensusMinNumMs2Scans") != decodedMap.end()){
        peaksSearchParameters->consensusMinNumMs2Scans = stoi(decodedMap["consensusMinNumMs2Scans"]);
    }
    if (decodedMap.find("consensusMinFractionMs2Scans") != decodedMap.end()){
        peaksSearchParameters->consensusMinFractionMs2Scans = stof(decodedMap["consensusMinFractionMs2Scans"]);
    }

    // ms1 matching
    if (decodedMap.find("ms1PpmTolr") != decodedMap.end()) {
        peaksSearchParameters->ms1PpmTolr = stof(decodedMap["ms1PpmTolr"]);
    }

    // ms2 search
    if (decodedMap.find("ms2MinNumMatches") != decodedMap.end()) {
        peaksSearchParameters->ms2MinNumMatches = stoi(decodedMap["ms2MinNumMatches"]);
    }
    if (decodedMap.find("ms2MinNumDiagnosticMatches") != decodedMap.end()) {
        peaksSearchParameters->ms2MinNumDiagnosticMatches = stoi(decodedMap["ms2MinNumDiagnosticMatches"]);
    }
    if (decodedMap.find("ms2MinNumUniqueMatches") != decodedMap.end()) {
        peaksSearchParameters->ms2MinNumUniqueMatches = stoi(decodedMap["ms2MinNumUniqueMatches"]);
    }
    if (decodedMap.find("ms2PpmTolr") != decodedMap.end()) {
        peaksSearchParameters->ms2PpmTolr = stof(decodedMap["ms2PpmTolr"]);
    }
    if (decodedMap.find("ms2MinIntensity") != decodedMap.end()) {
        peaksSearchParameters->ms2MinIntensity = stof(decodedMap["ms2MinIntensity"]);
    }

    // END SearchParameters
    // START PeaksSearchParameters

    //baseline
    if (decodedMap.find("baselineSmoothingWindow") != decodedMap.end()) {
        peaksSearchParameters->baselineSmoothingWindow = stof(decodedMap["baselineSmoothingWindow"]);
    }
    if (decodedMap.find("baselineDropTopX") != decodedMap.end()) {
        peaksSearchParameters->baselineDropTopX = stoi(decodedMap["baselineDropTopX"]);
    }

    //eic
    if (decodedMap.find("eicSmoothingWindow") != decodedMap.end()){
        peaksSearchParameters->eicSmoothingWindow = stof(decodedMap["eicSmoothingWindow"]);
    }
    if (decodedMap.find("eicEicSmoothingAlgorithm") != decodedMap.end()){
        peaksSearchParameters->eicEicSmoothingAlgorithm = decodedMap["eicEicSmoothingAlgorithm"];
    }
    if (decodedMap.find("eicMaxPeakGroupRtDiff") != decodedMap.end()){
        peaksSearchParameters->eicMaxPeakGroupRtDiff = stof(decodedMap["eicMaxPeakGroupRtDiff"]);
    }
    if (decodedMap.find("eicNoNoiseObs") != decodedMap.end()){
        peaksSearchParameters->eicNoNoiseObs = stoi(decodedMap["eicNoNoiseObs"]);
    }

    //quality
    if (decodedMap.find("qualitySignalBaselineRatio") != decodedMap.end()){
        peaksSearchParameters->qualitySignalBaselineRatio = stof(decodedMap["qualitySignalBaselineRatio"]);
    }
    if (decodedMap.find("qualitySignalBlankRatio") != decodedMap.end()){
        peaksSearchParameters->qualitySignalBlankRatio = stof(decodedMap["qualitySignalBlankRatio"]);
    }
    if (decodedMap.find("qualityMinPeakGroupIntensity") != decodedMap.end()){
        peaksSearchParameters->qualityMinPeakGroupIntensity = stof(decodedMap["qualityMinPeakGroupIntensity"]);
    }
    if (decodedMap.find("qualityMinPeakWidth") != decodedMap.end()){
        peaksSearchParameters->qualityMinPeakWidth = stoi(decodedMap["qualityMinPeakWidth"]);
    }
    if (decodedMap.find("qualityMinGoodPeakPerGroup") != decodedMap.end()){
        peaksSearchParameters->qualityMinGoodPeakPerGroup = stof(decodedMap["qualityMinGoodPeakPerGroup"]);
    }
    if (decodedMap.find("qualityMinPeakQuality") != decodedMap.end()){
        peaksSearchParameters->qualityMinPeakQuality = stof(decodedMap["qualityMinPeakQuality"]);
    }
    if (decodedMap.find("qualityClassifierModelName") != decodedMap.end()){
        peaksSearchParameters->qualityClassifierModelName = decodedMap["qualityClassifierModelName"];
    }

    //isotopes
    if (decodedMap.find("isotopesIsRequireMonoisotopicPeaks") != decodedMap.end()){
        peaksSearchParameters->isotopesIsRequireMonoisotopicPeaks = decodedMap["isotopesIsRequireMonoisotopicPeaks"] == "1";
    }
    if (decodedMap.find("isotopesExtractIsotopicPeaks") != decodedMap.end()) {
        peaksSearchParameters->isotopesExtractIsotopicPeaks = decodedMap["isotopesExtractIsotopicPeaks"] == "1";
    }
    if (decodedMap.find("isotopesMzTolerance") != decodedMap.end()){
        peaksSearchParameters->isotopesMzTolerance = stof(decodedMap["isotopesMzTolerance"]);
    }

    //ms1
    if (decodedMap.find("ms1IsMatchRtFlag") != decodedMap.end()){
        peaksSearchParameters->ms1IsMatchRtFlag = decodedMap["ms1IsMatchRtFlag"] == "1";
    }
    if (decodedMap.find("ms1RtTolr") != decodedMap.end()) {
        peaksSearchParameters->ms1RtTolr = stof(decodedMap["ms1RtTolr"]);
    }
    if (decodedMap.find("ms1MassSliceMergeNumScans") != decodedMap.end()) {
        peaksSearchParameters->ms1MassSliceMergeNumScans = stof(decodedMap["ms1MassSliceMergeNumScans"]);
    }
    if (decodedMap.find("ms1MassSliceMergePpm") != decodedMap.end()) {
        peaksSearchParameters->ms1MassSliceMergePpm = stof(decodedMap["ms1MassSliceMergePpm"]);
    }

    //ms2 matching
    if (decodedMap.find("ms2IsMatchMs2") != decodedMap.end()){
        peaksSearchParameters->ms2IsMatchMs2 = decodedMap["ms2IsMatchMs2"] == "1";
    }
    if (decodedMap.find("ms2ScoringAlgorithm") != decodedMap.end()){
        peaksSearchParameters->ms2ScoringAlgorithm = decodedMap["ms2ScoringAlgorithm"];
    }
    if (decodedMap.find("ms2MinScore") != decodedMap.end()){
        peaksSearchParameters->ms2MinScore = stof(decodedMap["ms2MinScore"]);
    }

    //matching options
    if (decodedMap.find("matchingIsRequireAdductPrecursorMatch") != decodedMap.end()){
        peaksSearchParameters->matchingIsRequireAdductPrecursorMatch = decodedMap["matchingIsRequireAdductPrecursorMatch"] == "1";
    }
    if (decodedMap.find("matchingIsRetainUnknowns") != decodedMap.end()){
        peaksSearchParameters->matchingIsRetainUnknowns = decodedMap["matchingIsRetainUnknowns"] == "1";
    }
    if (decodedMap.find("matchingLibraries") != decodedMap.end()) {
        peaksSearchParameters->matchingLibraries = decodedMap["matchingLibraries"];
    }
    if (decodedMap.find("matchingPolicy") != decodedMap.end()) {
        string matchingPolicyStr = decodedMap["matchingPolicy"];
        if (matchingPolicyStr == "ALL_MATCHES") {
            peaksSearchParameters->matchingPolicy = PeakGroupCompoundMatchingPolicy::ALL_MATCHES;
        } else if (matchingPolicyStr == "SINGLE_TOP_HIT") {
            peaksSearchParameters->matchingPolicy = PeakGroupCompoundMatchingPolicy::SINGLE_TOP_HIT;
        } else if (matchingPolicyStr == "TOP_SCORE_HITS") {
            peaksSearchParameters->matchingPolicy = PeakGroupCompoundMatchingPolicy::TOP_SCORE_HITS;
        }
    }

    // END PeaksSearchParameters

    return peaksSearchParameters;
}


set<mzSample*> SRMTransition::getSamples(){
       set<mzSample*> samples{};
       for (auto p : mzSlices) {
           samples.insert(p.first);
       }
       return samples;
   }

string IsotopeParameters::encodeParams() {

    string encodedParams;

    encodedParams = encodedParams + "searchVersion" + "=" + searchVersion + ";";

    encodedParams = encodedParams + "ppm" + "=" + to_string(ppm) + ";";
    encodedParams = encodedParams + "maxIsotopeScanDiff" + "=" + to_string(maxIsotopeScanDiff) + ";";
    encodedParams = encodedParams + "maxNaturalAbundanceErr" + "=" + to_string(maxNaturalAbundanceErr) + ";";
    encodedParams = encodedParams + "minIsotopicCorrelation" + "=" + to_string(minIsotopicCorrelation) + ";";
    encodedParams = encodedParams + "isC13Labeled" + "=" + to_string(isC13Labeled) + ";";
    encodedParams = encodedParams + "isN15Labeled" + "=" + to_string(isN15Labeled) + ";";
    encodedParams = encodedParams + "isS34Labeled" + "=" + to_string(isS34Labeled) + ";";
    encodedParams = encodedParams + "isD2Labeled" + "=" + to_string(isD2Labeled) + ";";

    encodedParams = encodedParams + "eic_smoothingAlgorithm" + "=" + to_string(eic_smoothingAlgorithm) + ";";
    encodedParams = encodedParams + "eic_smoothingWindow" + "=" + to_string(eic_smoothingWindow) + ";";

    encodedParams = encodedParams + "isIgnoreNaturalAbundance" + "=" + to_string(isIgnoreNaturalAbundance) + ";";
    encodedParams = encodedParams + "isExtractNIsotopes" + "=" + to_string(isExtractNIsotopes) + ";";
    encodedParams = encodedParams + "maxIsotopesToExtract" + "=" + to_string(maxIsotopesToExtract) + ";";
    encodedParams = encodedParams + "avgScanTime" + "=" + to_string(avgScanTime) + ";";

    if (adduct) {
        encodedParams = encodedParams + "adductName" + "=" + adduct->name + ";"; // use with Adduct* parameter
    } else {
        encodedParams = encodedParams + "adductName" + "=" + adductName + ";"; // use with Adduct* parameter
    }

    if (clsf) {
        encodedParams = encodedParams + "clsfFile" + "=" + clsf->getModelFilename() + ";"; // use with Classifier* parameter
    } else {
        encodedParams = encodedParams + "clsfFile" + "=" + clsfFile + ";";     // use with Classifier* parameter
    }

    return encodedParams;
}

IsotopeParameters IsotopeParameters::decode(string encodedParams) {

    IsotopeParameters isotopeParameters;

    unordered_map<string, string> decodedMap = mzUtils::decodeParameterMap(encodedParams);

    if (decodedMap.find("searchVersion") != decodedMap.end()) {
        isotopeParameters.searchVersion = decodedMap["searchVersion"];
    }
    if (decodedMap.find("ppm") != decodedMap.end()) {
        isotopeParameters.ppm = stof(decodedMap["ppm"]);
    }
    if (decodedMap.find("maxIsotopeScanDiff") != decodedMap.end()) {
        isotopeParameters.maxIsotopeScanDiff = stod(decodedMap["maxIsotopeScanDiff"]);
    }
    if (decodedMap.find("maxNaturalAbundanceErr") != decodedMap.end()) {
        isotopeParameters.maxNaturalAbundanceErr = stod(decodedMap["maxNaturalAbundanceErr"]);
    }
    if (decodedMap.find("minIsotopicCorrelation") != decodedMap.end()) {
        isotopeParameters.minIsotopicCorrelation = stod(decodedMap["minIsotopicCorrelation"]);
    }
    if (decodedMap.find("isC13Labeled") != decodedMap.end()) {
        isotopeParameters.isC13Labeled = decodedMap["isC13Labeled"]=="1";
    }
    if (decodedMap.find("isN15Labeled") != decodedMap.end()) {
        isotopeParameters.isN15Labeled = decodedMap["isN15Labeled"]=="1";
    }
    if (decodedMap.find("isS34Labeled") != decodedMap.end()) {
        isotopeParameters.isS34Labeled = decodedMap["isS34Labeled"]=="1";
    }
    if (decodedMap.find("isD2Labeled") != decodedMap.end()) {
        isotopeParameters.isD2Labeled = decodedMap["isD2Labeled"]=="1";
    }
    if (decodedMap.find("eic_smoothingAlgorithm") != decodedMap.end()) {
        isotopeParameters.eic_smoothingAlgorithm = static_cast<EIC::SmootherType>(stoi(decodedMap["eic_smoothingAlgorithm"]));
    }
    if (decodedMap.find("eic_smoothingWindow") != decodedMap.end()) {
        isotopeParameters.eic_smoothingWindow = stof(decodedMap["eic_smoothingWindow"]);
    }
    if (decodedMap.find("isIgnoreNaturalAbundance") != decodedMap.end()) {
        isotopeParameters.isIgnoreNaturalAbundance = decodedMap["isIgnoreNaturalAbundance"]=="1";
    }
    if (decodedMap.find("isExtractNIsotopes") != decodedMap.end()) {
        isotopeParameters.isExtractNIsotopes = decodedMap["isExtractNIsotopes"]=="1";
    }
    if (decodedMap.find("maxIsotopesToExtract") != decodedMap.end()) {
        isotopeParameters.maxIsotopesToExtract = stoi(decodedMap["maxIsotopesToExtract"]);
    }
    if (decodedMap.find("avgScanTime") != decodedMap.end()) {
        isotopeParameters.avgScanTime = stof(decodedMap["avgScanTime"]);
    }
    if (decodedMap.find("adductName") != decodedMap.end()) {
        isotopeParameters.adductName = decodedMap["adductName"];
    }
    if (decodedMap.find("clsfFile") != decodedMap.end()) {
        isotopeParameters.clsfFile = decodedMap["clsfFile"];
    }

    isotopeParameters.isotopeParametersType = IsotopeParametersType::SAVED;

    return isotopeParameters;
}
