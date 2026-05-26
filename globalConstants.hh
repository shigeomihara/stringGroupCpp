/**************************************************************************
 *        globalConstants.hh
 ************************************************************************/
#if !defined(__GLOBALCONSTANTS_HH)
#define __GLOBALCONSTANTS_HH
#include <string>
#include "searchHistoryAM.hh"

class GlobalConstants{
public:
    // ND240HA
    float IscSTC {8.61};   // A
    float VocSTC {37.05};   // V
    float ImpSTC {8.16};   // A
    float VmpSTC {29.42};   // V
    const float alphaSc {0.00038};   // /degC
    const float betaOc {-0.00329};   // /degC
    const float gammaMp {-0.00440};   // /degC

    // ND240HA
    const float IscSTC_ND240HA {8.61};   // A
    const float VocSTC_ND240HA {37.05};   // V
    const float ImpSTC_ND240HA {8.16};   // A
    const float VmpSTC_ND240HA {29.42};   // V

    // ND265HD
    const float IscSTC_ND265HD {8.70};   // A   ND265HD
    const float VocSTC_ND265HD {37.91};   // V
    const float ImpSTC_ND265HD {8.38};   // A
    const float VmpSTC_ND265HD {31.63};   // V
    
    const int Ns {60};
    const int Ns_sub {20};   // the number of cells in a substring
    
    const int S {84};   // the number of strings in a stringGroup
    const int L {14*3};   // the number of substrings in a string

    static const int K {4};   // ストリング群のパラメータ推定をする測定データ時間数

    // search history of adaptive model
    SearchHistoryAM sh;
    
    bool haveGPU {false};  // GPUを使わないときは、これをfalseにする

    // Constructors
    GlobalConstants(bool haveGPU);
    GlobalConstants(bool haveGPU, std::string kataban);
    
    float getVocSTCss(){ return VocSTC*Ns_sub/Ns; }
    float getVmpSTCss(){ return VmpSTC*Ns_sub/Ns; }
};

#endif
