#include "lipidsummarizationutils.h"
#include <QRegExp>
#include <QStringList>
#include <regex>

using namespace std;

/**
 * @brief LipidSummarizationUtils::getNameComponents
 * @param lipidName
 * --> must be a Level 4 (or above) lipid name (will not work on Level 3)
 * @return pair<string, vector<string>> where pair.first = class and pair.second = chains (in order)
 */
pair<string, vector<string>> LipidSummarizationUtils::getNameComponents(string lipidName){

//    QString lipidNameAsQString = QString(lipidName.c_str());

//    QRegExp rx("(\\,|\\(|\\)|\\/)");
//    QStringList qStringComponents = lipidNameAsQString.split(rx);

//    vector<string> chains;

//    if (qStringComponents.size() < 2) {
//        return make_pair("", chains); //no class information and no chains
//    }

//    for (unsigned int i = 1; i < qStringComponents.size(); i++) {
//        if (qStringComponents.at(i).toStdString().find(":") != string::npos){
//            chains.push_back(qStringComponents.at(i).toStdString());
//        }
//    }

//    pair<string, vector<string>> components = make_pair(qStringComponents.at(0).toStdString(), chains);


    regex rx("(\\,|\\(|\\)|\\/)");

    sregex_token_iterator iter(lipidName.begin(), lipidName.end(), rx, -1);

    vector<string> chains;
    string lipidClass("");

    sregex_token_iterator end;
    for(; iter != end; ++iter){
        string bit = *(iter);

        //first match is lipid class
        if (lipidClass == "") {
            lipidClass = bit;
        }

        //semicolon indicates acyl chain
        if (bit.find(":") != string::npos) {
            chains.push_back(bit);
        }
    }

    pair<string, vector<string>> components = make_pair(lipidClass, chains);
    return components;
}

/**
 * @brief LipidSummarizationUtils::getAcylChainLengthSummary
 * @param lipidName
 * @return Level 3 summary (acyl chain lengths, but not information about attachment to head group)
 */
string LipidSummarizationUtils::getAcylChainLengthSummary(string lipidName){
    string lipidNameNoPositionLevel = LipidSummarizationUtils::getSummary(LipidSummarizationUtils::getNameComponents(lipidName), 3);
    return lipidNameNoPositionLevel == "" ? lipidName : lipidNameNoPositionLevel;
}

/**
 * @brief LipidSummarizationUtils::getAcylChainCompositionSummary
 * @param lipidName
 * @return  Level 2 summary (# of C-C and C=C bonds)
 */
string LipidSummarizationUtils::getAcylChainCompositionSummary(string lipidName){
    string lipidNameChainSumLevel = LipidSummarizationUtils::getSummary(LipidSummarizationUtils::getNameComponents(lipidName), 2);
    return lipidNameChainSumLevel == "" ? lipidName : lipidNameChainSumLevel;
}

/**
 * @brief LipidSummarizationUtils::getLipidClassSummary
 * @param lipidName
 * @return Level 1 summary (lipid class name)
 */
string LipidSummarizationUtils::getLipidClassSummary(string lipidName){
    string lipidClass = LipidSummarizationUtils::getSummary(LipidSummarizationUtils::getNameComponents(lipidName), 1);
    return lipidClass == "" ? lipidName : lipidClass;
}

/**
 * @brief LipidSummarizationUtils::getSummarized
 * summarize the lipid name based on the information that is available in the name,
 * to the level of detail desired.
 *
 * Common input is Summary level 4, which includes the acyl chain composition and lengths,
 * but not the location(s) of the stereochemistry of the double bond(s).
 *
 * @param lipidNameComponents
 * @param summaryLevel
 *  Summary level 3: known chain legnths, head group attachment unknown
 *  Summary level 2: known composition of acyl chains (# of C-C and C=C bonds), chain lengths unknown
 *  Summary level 1: lipid class only (all chain information unknown)
 * @return
 */
string LipidSummarizationUtils::getSummary(pair<string, vector<string>> lipidNameComponents, int summaryLevel) {
    string lipidNameSummarized = string("");

    //Currently only summary levels 1-3 are supported
    if (summaryLevel < 1 || summaryLevel > 3) {
        return lipidNameSummarized;
    }

    //lipid class information is necessary for summary levels 1-3
    if (lipidNameComponents.first != ""){

        lipidNameSummarized.append(lipidNameComponents.first);

        if (summaryLevel == 1) { //lipid class only
            return lipidNameSummarized;
        }

        //chain information is necessary for summary levels 2 and 3
        if (lipidNameComponents.second.size() > 0) {

            lipidNameSummarized.append("(");

            if (summaryLevel == 3) {

                string linkageType = "";
                vector<string> chains;

                for (auto chain : lipidNameComponents.second){

                    QRegExp rx("-");
                    QString chainAsQString(chain.c_str());

                    QStringList qStringComponents = chainAsQString.split(rx);

                    if (qStringComponents.size() == 1) {
                        chains.push_back(chain);
                    } else if (qStringComponents.size() == 2){
                        linkageType = qStringComponents.at(0).toStdString().append("-");
                        chains.push_back(qStringComponents.at(1).toStdString());
                    }
                }

                if (linkageType != "") {
                    lipidNameSummarized.append(linkageType);
                }

                //WARNING: this sorting assumes that the vector contains no empty strings
                sort(chains.begin(), chains.end(),
                     [](const string& lhs, const string& rhs){

                    char lhsFirst = lhs.at(0);
                    char rhsFirst = rhs.at(0);

                    if (isdigit(lhsFirst) && isdigit(rhsFirst)) {
                        return lhs.compare(rhs) < 0; //TODO: look for digits up to colon, parse as number
                    } else if (isdigit(lhsFirst) && !isdigit(rhsFirst)) {
                        return lhs > rhs;
                    } else if (!isdigit(lhsFirst) && isdigit(rhsFirst)) {
                        return lhs < rhs;
                    } else {
                        if (lhsFirst == 't' && rhsFirst != 't') {
                            return lhs < rhs;
                        } else if (lhsFirst != 't' && rhsFirst == 't') {
                            return lhs > rhs;
                        } else {
                            if (lhsFirst == 'd' && rhsFirst != 'd') {
                                return lhs < rhs;
                            } else if (lhsFirst != 'd' && rhsFirst == 'd') {
                                return lhs > rhs;
                            } else {
                                if (lhsFirst == 'm' && rhsFirst != 'm') {
                                    return lhs < rhs;
                                } else if (lhsFirst != 'm' && rhsFirst == 'm') {
                                    return lhs > rhs;
                                } else {
                                    return lhs.compare(rhs) < 0; //TODO: look for digits up to colon, parse as number
                                }
                            }
                        }
                    }
                });

                for (unsigned int i = 0; i < chains.size(); i++){
                    if (i > 0) {
                        lipidNameSummarized.append("_");
                    }
                    lipidNameSummarized.append(chains.at(i));
                }

                lipidNameSummarized.append(")");

            } else if (summaryLevel == 2) {

                int alkaneSum = 0;
                int alkeneSum = 0;
                string linkageType = "";
                int numHydroxyl = 0;

                for (auto chain : lipidNameComponents.second){

                    vector<string> chainBits;

                    QRegExp rx(":");
                    QString chainQString(chain.c_str());

                    QStringList chainBitsRaw = chainQString.split(rx);

                    for (auto chainBit : chainBitsRaw) {

                        if (chainBit.startsWith("m", Qt::CaseSensitive)){
                            numHydroxyl++;
                            string chainBitStd = chainBit.toStdString();
                            chainBitStd = chainBitStd.substr(1, chainBitStd.size()-1);
                            chainBit = QString(chainBitStd.c_str());

                        } else if (chainBit.startsWith("d", Qt::CaseSensitive)) {
                            numHydroxyl = numHydroxyl + 2;
                            string chainBitStd = chainBit.toStdString();
                            chainBitStd = chainBitStd.substr(1, chainBit.size()-1);
                            chainBit = QString(chainBitStd.c_str());

                        } else if (chainBit.startsWith("t", Qt::CaseSensitive)) {
                            numHydroxyl = numHydroxyl + 3;
                            string chainBitStd = chainBit.toStdString();
                            chainBitStd = chainBitStd.substr(1, chainBit.size()-1);
                            chainBit = QString(chainBitStd.c_str());

                        }

                        //Issue 124: for p- and o- linked lipids
                        QRegExp rx("-");
                        QStringList chainBitPieces = chainBit.split(rx);

                        if (chainBitPieces.size() == 1) {
                            chainBits.push_back(chainBit.toStdString());
                        } else if (chainBitPieces.size() == 2){
                            linkageType = chainBitPieces.at(0).toStdString().append("-");
                            chainBits.push_back(chainBitPieces.at(1).toStdString());
                        }
                    }

                    //all chains must be sensible, otherwise do not try to to infer anything
                    if (chainBits.size() != 2) {
                      return string("");
                    }

                    alkaneSum += stoi(chainBits.at(0));
                    alkeneSum += stoi(chainBits.at(1));
                }

                if (linkageType != "") {
                    lipidNameSummarized.append(linkageType);
                }

                lipidNameSummarized.append(to_string(alkaneSum));
                lipidNameSummarized.append(":");
                lipidNameSummarized.append(to_string(alkeneSum));

                if (numHydroxyl > 0) {
                    lipidNameSummarized.append(",");
                    lipidNameSummarized.append(to_string(numHydroxyl));
                    lipidNameSummarized.append("-OH");
                }

                lipidNameSummarized.append(")");
            }

        }
    }

    return lipidNameSummarized;
}
