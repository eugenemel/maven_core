#ifndef LIPIDSUMMARIZATIONUTILS_H
#define LIPIDSUMMARIZATIONUTILS_H

#include <string>
#include <vector>

class LipidSummarizationUtils {

private:
    LipidSummarizationUtils() {} // prevent instantiation

public:
    static std::pair<std::string, std::vector<std::string>> getNameComponents(std::string lipidName);

    static std::string getAcylChainLengthSummary(std::string lipidName);
    static std::string getAcylChainCompositionSummary(std::string lipidName);
    static std::string getLipidClassSummary(std::string lipidName);


    static std::string getAcylChainLengthSummaryAttributeKey(){return "ACYL_CHAIN_LENGTH_SUMMARY";}
    static std::string getAcylChainCompositionSummaryAttributeKey() {return "ACYL_CHAIN_COMPOSITION_SUMMARY";}
    static std::string getLipidClassSummaryKey() {return "LIPID_CLASS";}

private:
    static std::string getSummary(std::pair<std::string, std::vector<std::string>> lipidNameComponents, int summaryLevel);

};

#endif // LIPIDSUMMARIZATIONUTILS_H
