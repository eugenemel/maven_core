#include "lipidsummarizationutils.h"
#include <QRegExp>
#include <QStringList>

using namespace std;

/**
 * @brief LipidUtils::getNameComponents
 * @param lipidName
 * --> must be a Level 4 (or above) lipid name (will not work on Level 3)
 * @return pair<string, vector<string>> where pair.first = class and pair.second = chains (in order)
 */
pair<string, vector<string>> LipidSummarizationUtils::getNameComponents(string lipidName){

    QString lipidNameAsQString = QString(lipidName.c_str());

    QRegExp rx("(\\, |\\( |\\) |\\/)");
    QStringList qStringComponents = lipidNameAsQString.split(rx, QString::SkipEmptyParts);

    vector<string> chains;

    if (qStringComponents.size() < 1) {
        return make_pair("", chains); //no class information and no chains
    }

    for (unsigned int i = 1; i < qStringComponents.size(); i++) {
        if (qStringComponents.at(i).toStdString().find(":") != string::npos){
            chains.push_back(qStringComponents.at(i).toStdString());
        }
    }

    pair<string, vector<string>> components = make_pair(qStringComponents.at(0).toStdString(), chains);
    return components;
}

string LipidSummarizationUtils::getAcylChainLengthSummary(string lipidName){
    return "TODO";
}

string LipidSummarizationUtils::getAcylChainCompositionSummary(string lipidName){
    return "TODO";
}

string LipidSummarizationUtils::getLipidClassSummary(string lipidName){
    return "TODO";
}

string LipidSummarizationUtils::getSummary(pair<string, vector<string>> lipidNameComponents, int summaryLevel) {
    return "TODO";
}
