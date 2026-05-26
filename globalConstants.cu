// -*- C++ -*-
/**************************************************************************
 *        globalConstants.cu
 ************************************************************************/
#include <string>
#include "globalConstants.hh"

/************ コンストラクタ   ************/
GlobalConstants::GlobalConstants(bool haveGPU): sh{haveGPU}{
    this->haveGPU = haveGPU;
}

GlobalConstants::GlobalConstants(bool haveGPU, std::string kataban): sh{haveGPU}{
    this->haveGPU = haveGPU;
    
    if(kataban == "ND240HA"){
	IscSTC = IscSTC_ND240HA;
	VocSTC = VocSTC_ND240HA;
	ImpSTC = ImpSTC_ND240HA;
	VmpSTC = VmpSTC_ND240HA;
    }else if(kataban == "ND265HD"){
	IscSTC = IscSTC_ND265HD;
	VocSTC = VocSTC_ND265HD;
	ImpSTC = ImpSTC_ND265HD;
	VmpSTC = VmpSTC_ND265HD;
    }
}
