// Issue 856: Standalone comparison driver for ExperimentAnchorPoints::compute().
//
// Loads the real X0497-M043_stds .mzML samples and the real anchor points
// file used by the mzkitchen run under investigation, runs the *unmodified*
// ExperimentAnchorPoints::compute() exactly as peakdetector.cpp's
// anchorPointsBasedAlignment() does, and reports, per anchor-point compound
// and sample, how far compute()'s rt_update lands from the RT a human
// manually verified in manual_bookmarks.mzrollDB.
//
// See $RUN_DIR/research/diagnosis.md for the full write-up this driver's
// output is evidence for.

#include "mzSample.h"
#include <sqlite3.h>
#include <algorithm>
#include <cmath>
#include <dirent.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

using namespace std;

namespace {

struct LibraryAnchor {
    double mz;
    double rt;
    string name;
};

string basename(const string& path) {
    size_t pos = path.find_last_of('/');
    return pos == string::npos ? path : path.substr(pos + 1);
}

string trim(string s) {
    while (!s.empty() && (s.back() == '\r' || s.back() == '\n' || s.back() == ' ')) s.pop_back();
    size_t start = 0;
    while (start < s.size() && s[start] == ' ') start++;
    return s.substr(start);
}

// Parses the real .apts anchor points file format used by this run:
//   "<mz>,<rt> #<compound name>"
// (ExperimentAnchorPoints::computeAnchorPointSetFromFile() only reads the
// mz/rt fields via a ',' split; the trailing "#name" comment is not
// retained by the class itself, so this driver re-parses it independently
// in order to label the output.)
vector<LibraryAnchor> parseAptsFile(const string& path) {
    vector<LibraryAnchor> result;
    ifstream in(path);
    if (!in.is_open()) {
        cerr << "Could not open anchor points file: " << path << endl;
        return result;
    }
    string line;
    while (getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        size_t commaPos = line.find(',');
        if (commaPos == string::npos) continue;
        string mzStr = line.substr(0, commaPos);
        string rest = line.substr(commaPos + 1);
        size_t hashPos = rest.find('#');
        string rtStr = hashPos == string::npos ? rest : rest.substr(0, hashPos);
        string name = hashPos == string::npos ? "" : trim(rest.substr(hashPos + 1));

        LibraryAnchor a;
        try {
            a.mz = stod(mzStr);
            a.rt = stod(rtStr);
        } catch (...) {
            continue;
        }
        a.name = name;
        result.push_back(a);
    }
    return result;
}

vector<string> findMzMLFiles(const string& dir) {
    vector<string> files;
    DIR* d = opendir(dir.c_str());
    if (!d) {
        cerr << "Could not open sample directory: " << dir << endl;
        return files;
    }
    struct dirent* entry;
    while ((entry = readdir(d)) != nullptr) {
        string name = entry->d_name;
        if (name.size() > 5 && name.substr(name.size() - 5) == ".mzML") {
            files.push_back(dir + "/" + name);
        }
    }
    closedir(d);
    sort(files.begin(), files.end());
    return files;
}

struct SqliteHandle {
    sqlite3* db = nullptr;
    explicit SqliteHandle(const string& path) {
        if (sqlite3_open_v2(path.c_str(), &db, SQLITE_OPEN_READONLY, nullptr) != SQLITE_OK) {
            cerr << "Could not open sqlite db (read-only): " << path << " (" << sqlite3_errmsg(db) << ")" << endl;
            sqlite3_close(db);
            db = nullptr;
        }
    }
    ~SqliteHandle() { if (db) sqlite3_close(db); }
};

// Ground truth extracted from manual_bookmarks.mzrollDB: for a given
// compound name and sample (matched by mzML basename), the RT a human
// verified in the MAVEN GUI, plus the peak intensity (useful for judging
// signal quality of the sample the RT was picked from).
struct ManualBookmarks {
    map<string, int> groupIdByCompound;                 // peakgroups.compoundName -> groupId
    map<string, int> sampleIdByFilename;                // basename(samples.filename) -> sampleId
    map<pair<int, int>, pair<double, double>> peakByGroupSample; // (groupId,sampleId) -> (rt, peakIntensity)

    bool lookup(const string& compoundName, const string& sampleFilename, double& rtOut, double& intensityOut) const {
        auto gIt = groupIdByCompound.find(compoundName);
        if (gIt == groupIdByCompound.end()) return false;
        auto sIt = sampleIdByFilename.find(basename(sampleFilename));
        if (sIt == sampleIdByFilename.end()) return false;
        auto pIt = peakByGroupSample.find({gIt->second, sIt->second});
        if (pIt == peakByGroupSample.end()) return false;
        rtOut = pIt->second.first;
        intensityOut = pIt->second.second;
        return true;
    }
};

ManualBookmarks loadManualBookmarks(const string& dbPath) {
    ManualBookmarks mb;
    SqliteHandle handle(dbPath);
    if (!handle.db) return mb;

    sqlite3_stmt* stmt = nullptr;
    if (sqlite3_prepare_v2(handle.db, "select sampleId, filename from samples", -1, &stmt, nullptr) == SQLITE_OK) {
        while (sqlite3_step(stmt) == SQLITE_ROW) {
            int sid = sqlite3_column_int(stmt, 0);
            const unsigned char* fname = sqlite3_column_text(stmt, 1);
            if (fname) mb.sampleIdByFilename[basename(reinterpret_cast<const char*>(fname))] = sid;
        }
    }
    sqlite3_finalize(stmt);

    stmt = nullptr;
    if (sqlite3_prepare_v2(handle.db, "select groupId, compoundName from peakgroups", -1, &stmt, nullptr) == SQLITE_OK) {
        while (sqlite3_step(stmt) == SQLITE_ROW) {
            int gid = sqlite3_column_int(stmt, 0);
            const unsigned char* name = sqlite3_column_text(stmt, 1);
            if (name) mb.groupIdByCompound[reinterpret_cast<const char*>(name)] = gid;
        }
    }
    sqlite3_finalize(stmt);

    stmt = nullptr;
    if (sqlite3_prepare_v2(handle.db, "select groupId, sampleId, rt, peakIntensity from peaks", -1, &stmt, nullptr) == SQLITE_OK) {
        while (sqlite3_step(stmt) == SQLITE_ROW) {
            int gid = sqlite3_column_int(stmt, 0);
            int sid = sqlite3_column_int(stmt, 1);
            double rt = sqlite3_column_double(stmt, 2);
            double intensity = sqlite3_column_double(stmt, 3);
            mb.peakByGroupSample[{gid, sid}] = {rt, intensity};
        }
    }
    sqlite3_finalize(stmt);

    return mb;
}

string findLibraryName(const vector<LibraryAnchor>& library, double mz, double rt) {
    for (auto& lib : library) {
        if (fabs(lib.mz - mz) < 1e-3 && fabs(lib.rt - rt) < 1e-3) return lib.name;
    }
    return "";
}

} // namespace

int main(int argc, char** argv) {

    // Defaults: the real X0497-M043_stds dataset and the real parameters
    // this dataset's mzkitchen run actually used for its executed module,
    // "peakdetector_mzkitchen_search" (see config.json in the run's output
    // folder). Override via flags for use against other datasets.
    string sampleDir = "/Users/phillipseitzer/Projects/mass_spec-issues/Issue2007-RTAlignment-X0497-DB/X0497-M043_stds";
    string anchorPointsFile = "/Users/phillipseitzer/Projects/mass_spec-issues/Issue2007-RTAlignment-X0497-DB/Pos_Data_Variability_DBRT_v2.apts";
    string manualBookmarksDb = "/Users/phillipseitzer/Projects/mass_spec-issues/Issue2007-RTAlignment-X0497-DB/manual_bookmarks.mzrollDB";
    float precursorPPM = 20;
    float maxRtWindow = 5;
    int eicSmoothingWindow = 5;
    float minPeakIntensity = 1e4; // peakdetector.cpp:88 default (standardsAlignment_minPeakIntensity is not set in config.json)
    bool isRtAlignmentAnchorPointReference = false; // config.json's peakdetector_mzkitchen_search: rtAlignmentAnchorPointReference=0
    bool debug = false;

    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        if (arg == "--samples" && i + 1 < argc) sampleDir = argv[++i];
        else if (arg == "--apts" && i + 1 < argc) anchorPointsFile = argv[++i];
        else if (arg == "--manual-db" && i + 1 < argc) manualBookmarksDb = argv[++i];
        else if (arg == "--ppm" && i + 1 < argc) precursorPPM = static_cast<float>(atof(argv[++i]));
        else if (arg == "--rt-window" && i + 1 < argc) maxRtWindow = static_cast<float>(atof(argv[++i]));
        else if (arg == "--smoothing" && i + 1 < argc) eicSmoothingWindow = atoi(argv[++i]);
        else if (arg == "--min-intensity" && i + 1 < argc) minPeakIntensity = static_cast<float>(atof(argv[++i]));
        else if (arg == "--anchor-point-reference") isRtAlignmentAnchorPointReference = true;
        else if (arg == "--debug") debug = true;
        else if (arg == "--help") {
            cout << "usage: rt_alignment_compare [--samples DIR] [--apts FILE] [--manual-db FILE]\n"
                 << "                             [--ppm N] [--rt-window N] [--smoothing N]\n"
                 << "                             [--min-intensity N] [--anchor-point-reference] [--debug]\n";
            return 0;
        }
    }

    cout << "=== rt_alignment_compare ===" << endl;
    cout << "samples dir:          " << sampleDir << endl;
    cout << "anchor points file:   " << anchorPointsFile << endl;
    cout << "manual bookmarks db:  " << manualBookmarksDb << endl;
    cout << "precursorPPM:         " << precursorPPM << endl;
    cout << "maxRtWindow:          " << maxRtWindow << endl;
    cout << "eic_smoothingWindow:  " << eicSmoothingWindow << endl;
    cout << "minPeakIntensity:     " << minPeakIntensity << endl;
    cout << "anchorPointReference: " << (isRtAlignmentAnchorPointReference ? "true" : "false") << endl;
    cout << endl;

    // 1. Load samples, exactly as peakdetector.cpp's loadSamples() does when
    // no sampleIdMappingFile is given: natural-name sort, sequential sampleId.
    // (peakdetector.cpp:1521-1539 -- this run's config.json leaves
    // sampleIdMappingFile empty, so isAnchorPointSample is left at its
    // mzSample.cpp:31 default of `true` for every sample.)
    vector<string> mzmlFiles = findMzMLFiles(sampleDir);
    if (mzmlFiles.empty()) {
        cerr << "No .mzML files found in " << sampleDir << endl;
        return 1;
    }

    vector<mzSample*> samples;
    for (auto& f : mzmlFiles) {
        cout << "Loading " << f << " ..." << endl;
        mzSample* s = new mzSample();
        s->loadSample(f.c_str());
        if (s->scans.empty()) {
            cerr << "  WARNING: 0 scans loaded, skipping" << endl;
            delete s;
            continue;
        }
        samples.push_back(s);
    }
    if (samples.empty()) {
        cerr << "No samples with scan data were loaded. Aborting." << endl;
        return 1;
    }

    sort(samples.begin(), samples.end(), [](mzSample* a, mzSample* b) { return a->sampleName < b->sampleName; });
    for (size_t i = 0; i < samples.size(); i++) {
        samples[i]->setSampleId(static_cast<int>(i));
        samples[i]->setSampleOrder(static_cast<int>(i));
    }

    cout << endl << "Loaded " << samples.size() << " samples (sampleId assigned by natural name order):" << endl;
    for (auto* s : samples) {
        cout << "  [" << s->sampleId << "] " << s->sampleName
             << " isAnchorPointSample=" << (s->isAnchorPointSample ? "true" : "false") << endl;
    }
    cout << endl;

    vector<LibraryAnchor> library = parseAptsFile(anchorPointsFile);
    cout << "Parsed " << library.size() << " named anchor points from " << anchorPointsFile << endl << endl;

    // 2. Run the real, unmodified ExperimentAnchorPoints::compute(), exactly
    // as peakdetector.cpp:2298-2314's anchorPointsBasedAlignment() does.
    ExperimentAnchorPoints experiment(
        samples,
        anchorPointsFile,
        precursorPPM,
        maxRtWindow,
        eicSmoothingWindow,
        minPeakIntensity,
        isRtAlignmentAnchorPointReference);

    experiment.compute(debug, true); // (debug, isClean) -- matches peakdetector.cpp:2311

    if (!experiment.referenceSample) {
        cerr << "compute() did not select a reference sample. Aborting." << endl;
        return 1;
    }
    cout << "Reference sample chosen by determineReferenceSample(): "
         << experiment.referenceSample->sampleName << endl << endl;

    // 3. Load manual_bookmarks.mzrollDB ground truth.
    ManualBookmarks manual = loadManualBookmarks(manualBookmarksDb);
    if (manual.groupIdByCompound.empty()) {
        cerr << "Could not load compound identities from " << manualBookmarksDb << ". Aborting." << endl;
        return 1;
    }

    // 4. Per anchor point (compound) x sample: observed RT, rt_update
    // (the reference sample's own RT for this anchor point -- this is
    // exactly what ends up in rt_update_key.rt_update), the manually
    // verified RT for that same compound in that same physical sample, and
    // the absolute deviation between rt_update and the manual RT.
    cout << fixed << setprecision(4);
    cout << left << setw(28) << "compound"
         << setw(20) << "sample"
         << right << setw(10) << "observed"
         << setw(12) << "rt_update"
         << setw(12) << "manual_rt"
         << setw(14) << "|deviation|"
         << "  from_eic" << endl;

    vector<double> deviations;
    int rowsWithoutManualMatch = 0;
    int rowsSkippedNoCompoundIdentity = 0;

    for (auto& aps : experiment.anchorPointSets) {
        if (!aps.isValid) continue;

        string compoundName = findLibraryName(library, aps.referenceMz, aps.referenceRt);
        if (compoundName.empty()) {
            // This is the synthetic "last RT in file" anchor point
            // (AnchorPointSet::lastRt(), mzAligner.cpp:765) -- it has no
            // compound identity in manual_bookmarks, so it is out of scope
            // for the R1/R4 compound-identity comparison (documented in
            // diagnosis.md).
            rowsSkippedNoCompoundIdentity++;
            continue;
        }

        auto refIt = aps.sampleToPoints.find(experiment.referenceSample);
        if (refIt == aps.sampleToPoints.end()) continue;
        double rtUpdate = refIt->second->rt; // what rt_update_key.rt_update will be for every sample, for this anchor point

        for (auto* sample : samples) {
            auto it = aps.sampleToPoints.find(sample);
            if (it == aps.sampleToPoints.end()) continue;

            double observedRt = it->second->rt;
            bool fromEIC = it->second->isRtFromEIC;

            double manualRt = 0.0, manualIntensity = 0.0;
            bool hasManual = manual.lookup(compoundName, sample->sampleName, manualRt, manualIntensity);

            cout << left << setw(28) << compoundName.substr(0, 27)
                 << setw(20) << sample->sampleName.substr(0, 19)
                 << right << setw(10) << observedRt
                 << setw(12) << rtUpdate
                 << setw(12);
            if (hasManual) {
                cout << manualRt;
                double deviation = fabs(rtUpdate - manualRt);
                cout << setw(14) << deviation;
                deviations.push_back(deviation);
            } else {
                cout << "n/a";
                cout << setw(14) << "n/a";
                rowsWithoutManualMatch++;
            }
            cout << "  " << (fromEIC ? "yes" : "interp") << endl;
        }
    }

    cout << endl;
    cout << rowsSkippedNoCompoundIdentity << " anchor-point-set row(s) skipped: no compound identity (e.g. the synthetic last-RT anchor)." << endl;
    cout << rowsWithoutManualMatch << " row(s) had no matching manual_bookmarks peak (compound+sample not found there)." << endl;

    if (!deviations.empty()) {
        sort(deviations.begin(), deviations.end());
        double sum = 0;
        for (double d : deviations) sum += d;
        double mean = sum / deviations.size();
        double median = deviations.size() % 2 == 0
            ? (deviations[deviations.size() / 2 - 1] + deviations[deviations.size() / 2]) / 2.0
            : deviations[deviations.size() / 2];

        cout << endl << "compute(): n=" << deviations.size()
             << " mean(|rt_update - manual_rt|)=" << mean
             << " median(|rt_update - manual_rt|)=" << median
             << " max=" << deviations.back() << endl;
    } else {
        cerr << "No comparable rows -- could not compute mean/median deviation." << endl;
        return 1;
    }

    return 0;
}
