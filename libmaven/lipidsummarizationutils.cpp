#include "lipidsummarizationutils.h"
#include <regex>

using namespace std;

/**
 * @brief LipidSummarizationUtils::getNameComponents
 * @param lipidName
 * --> must be a Level 4 (or above) lipid name (will not work on Level 3)
 * @return pair<string, vector<string>> where pair.first = class and pair.second = chains (in order)
 *
 * the lipid class is assumed to be everything from the start of the string until the first open paranthesis (
 * the rest of the string are the chains.
 * If the last character is a ")", it is removed.
 *
 * The chains are split based on the forward slash / character.
 */
LipidNameComponents LipidSummarizationUtils::getNameComponents(string lipidName){

    int initialLevel = 4;

    LipidNameComponents lipidNameComponents;

    bool isFoundChains = false;
    string::size_type chainStart = lipidName.size();
    for (string::size_type i = 0; i < lipidName.size(); i++){
        if (lipidName[i] == '(') {
            chainStart = i;
            isFoundChains = true;
            break;
        }
    }

    if (!isFoundChains) {
        return lipidNameComponents;
    }

    string lipidClass = lipidName.substr(0, chainStart);

    string::size_type lastPos = lipidName.size()-1;
    if (lipidName[lastPos] == ')'){
        lastPos--;
    }

    string chainsSubstring = lipidName.substr(chainStart+1, (lastPos-chainStart));

    regex rx("/");

    sregex_token_iterator iter(chainsSubstring.begin(), chainsSubstring.end(), rx, -1);

    vector<string> chains;

    sregex_token_iterator end;
    for(; iter != end; ++iter){
        string chain = *(iter);
        chains.push_back(chain);

        //Issue 321: additional oxygenation information
        if (
                chain.find(";OH") != string::npos ||
                chain.find(";(OH") != string::npos || //handles any number of (OH)
                chain.find(";Ep") != string::npos ||
                chain.find(";OO") != string::npos || //handles ;OO and ;OOH
                chain.find(";oxo") != string::npos ||
                chain.find(";COOH") != string::npos
             ) {
            initialLevel = 5;
        }
    }

    lipidNameComponents.lipidClass = lipidClass;
    lipidNameComponents.chains = chains;
    lipidNameComponents.initialLevel = initialLevel;

    return lipidNameComponents;
}

/**
 * @brief LipidSummarizationUtils::getSnPositionSummary
 * @param lipidName
 * @return Level 4 summary (sn position, with oxygenation information summarized)
 */
string LipidSummarizationUtils::getSnPositionSummary(std::string lipidName){

    LipidNameComponents lipidNameComponents = LipidSummarizationUtils::getNameComponents(lipidName);

    if (lipidNameComponents.initialLevel == 4) {
        return lipidName;
    }

    return getSummaryLevel5ToLevel4(lipidNameComponents);
}

/**
 * @brief LipidSummarizationUtils::getAcylChainLengthSummary
 * @param lipidName
 * @return Level 3 summary (acyl chain lengths, but not information about attachment to head group)
 */
string LipidSummarizationUtils::getAcylChainLengthSummary(string lipidName){

    LipidNameComponents lipidNameComponents = LipidSummarizationUtils::getNameComponents(lipidName);

    //Issue 321: If the initial lipid is at level 5, first downgrade to level 4.
    if (lipidNameComponents.initialLevel == 5) {
        lipidName = LipidSummarizationUtils::getSummaryLevel5ToLevel4(lipidNameComponents);
        lipidNameComponents = LipidSummarizationUtils::getNameComponents(lipidName);
    }

    string summarizedLipidName = LipidSummarizationUtils::getSummary(LipidSummarizationUtils::getNameComponents(lipidName), 3);
    return summarizedLipidName == "" ? lipidName : summarizedLipidName;
}

/**
 * @brief LipidSummarizationUtils::getAcylChainCompositionSummary
 * @param lipidName
 * @return  Level 2 summary (# of C-C and C=C bonds)
 */
string LipidSummarizationUtils::getAcylChainCompositionSummary(string lipidName){

    LipidNameComponents lipidNameComponents = LipidSummarizationUtils::getNameComponents(lipidName);

    //Issue 321: If the initial lipid is at level 5, first downgrade to level 4.
    if (lipidNameComponents.initialLevel == 5) {
        lipidName = LipidSummarizationUtils::getSummaryLevel5ToLevel4(lipidNameComponents);
        lipidNameComponents = LipidSummarizationUtils::getNameComponents(lipidName);
    }

    string summarizedLipidName = LipidSummarizationUtils::getSummary(LipidSummarizationUtils::getNameComponents(lipidName), 2);
    return summarizedLipidName == "" ? lipidName : summarizedLipidName;
}

/**
 * @brief LipidSummarizationUtils::getLipidClassSummary
 * @param lipidName
 * @return Level 1 summary (lipid class name)
 */
string LipidSummarizationUtils::getLipidClassSummary(string lipidName){

    LipidNameComponents lipidNameComponents = LipidSummarizationUtils::getNameComponents(lipidName);

    return lipidNameComponents.lipidClass == "" ? lipidName : lipidNameComponents.lipidClass;
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
string LipidSummarizationUtils::getSummary(LipidNameComponents lipidNameComponents, int summaryLevel) {

    //Currently only summary levels 1-3 are supported
    if (summaryLevel < 1 || summaryLevel > 3) {
        return "";
    }

    //lipid class information is necessary for summary levels 1-3
    if (lipidNameComponents.lipidClass != ""){

        if (summaryLevel == 1) { //lipid class only
            return lipidNameComponents.lipidClass;
        }

        //chain information is necessary for summary levels 2 and 3
        if (lipidNameComponents.chains.size() > 0) {

            if (summaryLevel == 3) {

                return (getSummaryLevel4ToLevel3(lipidNameComponents));

            } else if (summaryLevel == 2) {

                return (getSummaryLevel4ToLevel2(lipidNameComponents));
            }

        }
    }

    //fallback to empty string, indicates that summarization was not possible
    return "";
}

string LipidSummarizationUtils::getSummaryLevel5ToLevel4(LipidNameComponents lipidNameComponents){

   string lipidNameSummarized("");
   lipidNameSummarized.append(lipidNameComponents.lipidClass);
   lipidNameSummarized.append("(");

   vector<string> updatedChains = vector<string>(lipidNameComponents.chains.size());

   regex rx(";");
   for (unsigned int i = 0; i < updatedChains.size(); i++) {

       string chain = lipidNameComponents.chains[i];

       sregex_token_iterator iter(chain.begin(), chain.end(), rx, -1);
       vector<string> oxBits{};

       sregex_token_iterator end;
       for (; iter != end; ++iter){
           string component = *(iter);
           oxBits.push_back(component);
       }

       if (oxBits.size() == 1) {
           //no oxidation of any kind

            updatedChains[i] = oxBits[0];
       } else {

           int numAddOx = 0;

           int numSingle = 0;
           int numDouble = 0;

           for (string::size_type pos = 0; pos < oxBits[0].size(); pos++) {
               char c = oxBits[0][pos];
               if (c == ':') {
                   numSingle = stoi(oxBits[0].substr(0, pos));
                   numDouble = stoi(oxBits[0].substr(pos+1, (oxBits[0].size()-pos)));
                   break;
               }
           }

           for (unsigned int j = 1; j < oxBits.size(); j++) {

               string oxBit = oxBits[j];

               //single hydroxyl
               if (oxBit == "OH") {
                   numAddOx++;
               }

               //multiple hydroxyl
               auto multipleHydroxylPos = oxBit.find("(OH)");
               if (multipleHydroxylPos != string::npos) {
                   int numHydroxyl = stoi(oxBit.substr(multipleHydroxylPos+4, oxBit.size()-multipleHydroxylPos));
                   numAddOx += numHydroxyl;
               }

               //epoxide
               if (oxBit == "Ep") {
                   numAddOx++;
                   numDouble++;
               }

               //hydroperoxide
               if (oxBit == "OOH") {
                   numAddOx = numAddOx + 2;
               }

               //peroxide
               if (oxBit == "OO") {
                   numAddOx = numAddOx + 2;
                   numDouble++;
               }

               //aldehyde
               if (oxBit == "oxo") {
                   numAddOx++;
                   numDouble++;
               }

               //carboxylic acid
               if (oxBit == "COOH") {
                   numAddOx = numAddOx + 2;
               }

           }

           string updatedChain = to_string(numSingle).append(":").append(to_string(numDouble));
           if (numAddOx > 0) {
               updatedChain = updatedChain.append(";O");
               if (numAddOx > 1) {
                   updatedChain = updatedChain.append(to_string(numAddOx));
               }
           }

           updatedChains[i] = updatedChain;
       }
   }

   for (unsigned int i = 0; i < updatedChains.size(); i++) {
       if (i > 0) {
           lipidNameSummarized.append("/");
       }
       lipidNameSummarized.append(updatedChains[i]);
   }
   lipidNameSummarized.append(")");

   return lipidNameSummarized;
}

string LipidSummarizationUtils::getSummaryLevel4ToLevel3(LipidNameComponents lipidNameComponents) {

    string lipidNameSummarized("");
    lipidNameSummarized.append(lipidNameComponents.lipidClass);
    lipidNameSummarized.append("(");

    string linkageType = "";
    vector<string> chains;

    for (auto chain : lipidNameComponents.chains){

        regex rx("-");
        sregex_token_iterator iter(chain.begin(), chain.end(), rx, -1);

        vector<string> chainComponents{};

        sregex_token_iterator end;

        for (; iter != end; ++iter){
            string component = *(iter);
            chainComponents.push_back(component);
        }

        if (chainComponents.size() == 1) {
            chains.push_back(chain);
        } else if (chainComponents.size() == 2){
            linkageType = chainComponents[0].append("-");
            chains.push_back(chainComponents[1]);
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

    return lipidNameSummarized;
}

string LipidSummarizationUtils::getSummaryLevel4ToLevel2(LipidNameComponents lipidNameComponents) {

    string lipidNameSummarized("");
    lipidNameSummarized.append(lipidNameComponents.lipidClass);
    lipidNameSummarized.append("(");

    int alkaneSum = 0;
    int alkeneSum = 0;
    string linkageType = "";
    int numOxygenations = 0;

    for (auto chain : lipidNameComponents.chains){

        vector<string> chainBits;

        regex rx(":");
        sregex_token_iterator iter(chain.begin(), chain.end(), rx, -1);
        vector<string> chainBitsRaw{};

        sregex_token_iterator end;

        for (; iter != end; ++iter){
            string component = *(iter);
            chainBitsRaw.push_back(component);
        }

        for (auto chainBit : chainBitsRaw) {

            string chainBitOHCorrected = chainBit;

            //Issue 321: Support old approach

            if (chainBit[0] == 'm'){

                numOxygenations++;
                chainBitOHCorrected = chainBit.substr(1);

            } else if (chainBit[0] == 'd') {

                numOxygenations = numOxygenations + 2;
                chainBitOHCorrected = chainBit.substr(1);

            } else if (chainBit[0] == 't') {

                numOxygenations = numOxygenations + 3;
                chainBitOHCorrected = chainBit.substr(1);
            }

            //Issue 321: oxygenations in format of x:y;On e.g. 16:1;O2 (16 carbons, 1 C=C, 2 oxygenations)

            int numOxSemicolonFormat = 0;
            for (string::size_type i = 0; i < chainBitOHCorrected.size(); i++) {
                if (chainBitOHCorrected[i] == ';') {
                    if (i < chainBitOHCorrected.size()-2) {
                        string::size_type end = chainBitOHCorrected.size()-i-1;
                        numOxSemicolonFormat = stoi(chainBitOHCorrected.substr(i+2, end));
                    } else {
                        numOxSemicolonFormat = 1;
                    }
                    break;
                }
            }

            numOxygenations += numOxSemicolonFormat;


            //Issue 124: for plasmalogens (e.g.PE(p-16:0/18:2) or PE(o-16:0/18:2))
            regex rx("-");
            sregex_token_iterator iter(chainBitOHCorrected.begin(), chainBitOHCorrected.end(), rx, -1);

            vector<string> chainBitPieces{};

            sregex_token_iterator end;

            for (; iter != end; ++iter){
                string component = *(iter);
                chainBitPieces.push_back(component);
            }

            if (chainBitPieces.size() == 1) {
                chainBits.push_back(chainBitOHCorrected);
            } else if (chainBitPieces.size() == 2) {
                linkageType = chainBitPieces[0].append("-");
                chainBits.push_back(chainBitPieces[1]);
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

    //Issue 321: Switch to lipidmaps 2020 syntax
    if (numOxygenations > 0) {
        lipidNameSummarized.append(";O");
        if (numOxygenations > 1) {
            lipidNameSummarized.append(to_string(numOxygenations));
        }
    }

    lipidNameSummarized.append(")");

    return lipidNameSummarized;
}
