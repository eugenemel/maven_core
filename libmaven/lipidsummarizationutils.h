#ifndef LIPIDSUMMARIZATIONUTILS_H
#define LIPIDSUMMARIZATIONUTILS_H

#include <string>
#include <vector>

class LipidSummarizationUtils {

private:
    LipidSummarizationUtils() {} // prevent instantiation

public:
    static std::pair<std::string, std::vector<std::string>> getNameComponents(std::string lipidName);

    //Level 3 (Molecular Species Level)
    static std::string getAcylChainLengthSummary(std::string lipidName);

    //Level 2 (Species Level)
    static std::string getAcylChainCompositionSummary(std::string lipidName);

    //Level 1 (Lipid Class)
    static std::string getLipidClassSummary(std::string lipidName);

    //LipidMaps 2020: "Full structure level ---> -1 (avoid giving level number for future levels)
    //this level doesn't have a proper number because it isn't summarized.
    static std::string getFullStructureLevelSummaryAttributeKey(){return "FULL_STRUCTURE_SUMMARY";}

    //LipidMaps 2020: "Structure defined level" --> level 5
    static std::string getStructureDefinedLevelSummaryAttributeKey(){return "STRUCTURE_DEFINED_SUMMARY";}

    //LipidMaps 2020: "sn-Position level" --> level 4
    static std::string getSnPositionLevelSummaryAttributeKey(){return "SN_POSITION_SUMMARY";}

    //LipidMaps 2020: "Molecular Species Level" ---> level 3
    static std::string getAcylChainLengthSummaryAttributeKey(){return "ACYL_CHAIN_LENGTH_SUMMARY";}

    //LipidMaps 2020: "Species Level" ---> level 2
    static std::string getAcylChainCompositionSummaryAttributeKey() {return "ACYL_CHAIN_COMPOSITION_SUMMARY";}

    //Lipid Class ---> level 1
    static std::string getLipidClassSummaryKey() {return "LIPID_CLASS";}

private:
    static std::string getSummary(std::pair<std::string, std::vector<std::string>> lipidNameComponents, int summaryLevel);

    static std::string getSummaryLevel4ToLevel3(std::pair<std::string, std::vector<std::string>> lipidNameComponents);
    static std::string getSummaryLevel4ToLevel2(std::pair<std::string, std::vector<std::string>> lipidNameComponents);

};

#endif // LIPIDSUMMARIZATIONUTILS_H
