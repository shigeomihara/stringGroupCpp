/**********************************************************************************
 *        PCSData.hh  2024.1.23--
 *********************************************************************************/
#if !defined(__PCS_DATA_HH)
#define __PCS_DATA_HH

#include <vector>
#include "globalConstants.hh"
#include "stringGroup.hh"

struct SimPCSData{
std::vector<double> G;
    std::vector<double> T;
    std::vector<double> I;
    std::vector<double> V;
    std::vector<std::string> Time;
};

class PCSData{
private:
    GlobalConstants &m_gc;
public:
    PCSData(GlobalConstants &gc);
    std::vector<double> G;
    std::vector<double> T;
    std::vector<double> I;
    std::vector<double> V;
    std::vector<std::string> Time;
    
    void readFile();
    void readFile(std::string fileName);
    
    void averageGTIV();
    void makeSimPCSDataPmax(SimPCSData &sd);
    void addFault(StringGroup &sg, const double G, const double T);
    
    static std::vector<std::string> splitString(std::string str, char c);
};

#endif
