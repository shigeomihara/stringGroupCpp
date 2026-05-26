// -*- C++ -*-
/************************************************************************************
 *                 Substring.cu
 ***********************************************************************************/
// #include "globalConstants.hh"
#include "substring.hh"
// #include "adaptiveModel.hh"
// #include "bisection.hh"
// #include <cfloat>  // DBL_MAX
#include <stdexcept>
#include <iostream>
#include <string>

/** Constructor     ****************************************************************
 *  出力：m_gc, m_amInit, m_G, m_T, m_Tkel, m_kT, m_etabVt, m_Ib0, m_lnIb0
 *  関数呼び出し：AdaptiveModelのコンストラクター    *******************************/
Substring::Substring(GlobalConstants &gcc, double GG, double TT): m_gc{gcc}, m_amInit{gcc,GG,TT},
							       m_G{GG}, m_T{TT}{
    m_Tkel = TT+273.15;
    m_kT = 1.380649E-23*m_Tkel;
    m_etabVt = 1.1*m_kT/1.602176634E-19;
    m_Ib0 = m_gc.IscSTC*(1.0+m_gc.alphaSc*(TT-25.0))/( exp(0.6/(m_etabVt))-1.0 );
    m_lnIb0 = log(m_Ib0);

    m_ImaxMax = m_amInit.getImax();  // Mar.4,2026
    m_VmaxInit = m_amInit.getVmax();
    m_VmaxSum = m_VmaxInit*m_gc.Ns_sub;  // Mar.4,2026
    m_ImaxMin = m_amInit.getImax();  // Mar.4,2026
    m_nVmin = m_gc.Ns_sub*m_amInit.getVmin();
    m_VbMin = -m_etabVt*std::log(m_IbMax/m_Ib0+1.0);

    m_VbypassDiodeOn = -m_etabVt*log(1.0E-3/m_Ib0 +1);  // 1mA以上ならオンとみなす
}

/** Constructor    ******************************************************************
 *  出力：m_gc, m_amInit, m_G, m_T, m_Tkel, m_kT, m_etabVt, m_Ib0, m_lnIb0
 *  関数呼び出し：AdaptiveModelのコンストラクター    *******************************/
Substring::Substring(GlobalConstants gcc, double GG, double TT,
		     double IscSTC, double VocSTC, double ImpSTC, double VmpSTC)
    : m_gc{gcc}, m_amInit{IscSTC, VocSTC, ImpSTC, VmpSTC,
		      gcc.alphaSc, gcc.betaOc, gcc.gammaMp, GG, TT, gcc.Ns_sub, gcc},
      m_G{GG}, m_T{TT}{
    m_Tkel = TT+273.15;
    m_kT = 1.380649E-23*m_Tkel;
    m_etabVt = 1.1*m_kT/1.602176634E-19;
    m_Ib0 = m_gc.IscSTC*(1.0+m_gc.alphaSc*(TT-25.0))/( exp(0.6/(m_etabVt))-1.0 );
    m_lnIb0 = log(m_Ib0);

    double Emin = m_amInit.calc();
    
    // printf("For datasheet values, Rs=%g, eta=%g, Rh=%g, I0=%g, Iph=%g, Emin=%g Substring::Substring\n",
    // 	   m_amInit.getRs(),m_amInit.getEta(), m_amInit.getRh(), m_amInit.getI0(), m_amInit.getIph(), Emin);//////////t
}

/** Calculate substring current I (A) for voltage V (V)    **************************
 *  入力： m_Ib0, m_etabVt
 *  出力： m_bypassDiodeOn
 *  関数呼び出し： AdaptiveModel::getI     *****************************************/
double Substring::getI(double V){
    double Ic;
    if(m_nVmin <= V && V <= m_VmaxSum){
	double Ic1 {m_amInit.getImin()};
	double Ic2 {m_ImaxMax};
	m_Vfunc2 = V;
	Bisection<Substring> bis {this};
	Ic = bis.search2(Ic1, Ic2);
    }else if(V < m_nVmin) Ic = m_ImaxMin;
    else Ic = m_amInit.getImin();
    
    double Ib;
    if(V < m_VbMin) Ib = m_IbMax;
    else Ib = m_Ib0*( exp(-V/m_etabVt)-1.0 );
    
    if(Ib > m_Ib0) m_bypassDiodeOn = true;
    else  m_bypassDiodeOn = false;
    return Ic+Ib;
}

/***************************************************************/
double Substring::func2(double Ic) const{
    double Vcells {m_amInit.getV(Ic)*(m_gc.Ns_sub-m_vamNumSum)};
    for(int i=0; i<m_vam.size(); i++) Vcells += m_vam[i].getV(Ic)*m_vamNum[i];
    return Vcells - m_Vfunc2;
}

/** Calculate substring voltage V (V) for current I (A)    **************************
 *  入力： m_Ib0, m_etabVt, m_lnIb0
 *  出力： m_bypassDiodeOn, m_Ifunc
 *  関数呼び出し： AdaptiveModel::getI, AdaptiveModel::getV, Bisection::search, func ***/
double Substring::getV(double I){
    m_bypassDiodeOn = false;
    
    // double x1 {cuda::std::numeric_limits<double>::min()};
    double x1 {std::numeric_limits<double>::min()};
    double x2;
    if(I>=0.0) x2 = I/m_Ib0+1.0;
    else x2 = 1.0;

    m_Ifunc = I+m_Ib0;
    
    Bisection<Substring> bis {this};
    double x = bis.search(x1, x2);

    double Ix {m_Ifunc-x*m_Ib0};
    double V {m_amInit.getV(Ix)*(m_gc.Ns_sub-m_vamNumSum)};
    for(int i=0; i<m_vam.size(); i++) V += m_vam[i].getV(Ix)*m_vamNum[i];
    
    if(V <= m_VbypassDiodeOn) m_bypassDiodeOn = true;
    return V;
}

/** getVで使われる２分法用のメソッド  **********************************************
 *  入力： m_Ib0, m_etabVt, m_lnIb0, m_Ifunc
 *  関数呼び出し： AdaptiveModel::getV    *****************************************/
double Substring::func(double x) const{
    double Ix {m_Ifunc-x*m_Ib0};
    double V {m_amInit.getV(Ix)*(m_gc.Ns_sub-m_vamNumSum)};
    for(int i=0; i<m_vam.size(); i++) V += m_vam[i].getV(Ix)*m_vamNum[i];
    return V + m_etabVt*std::log(x);
}

/** m_aG[]、m_aT[] と m_aEtabVt[], m_aIb0[], m_aLnIb0[] をセットする  ****************
 *  AdaptiveModel m_amInitのパラメータも対応してセットする
 *  出力： m_aG[], m_aT[], m_aEtabVt[], m_aIb0[], m_aLnIb0[]
 *  関数呼び出し： AdaptiveModel::calcPCS
 *  戻り値： AdaptiveModel::calcのE値     *******************************************/
double Substring::setVars(double GG[], double TT[]){
    for(int k=0; k<m_gc.K; k++){
	m_aG[k] = GG[k];
	m_aT[k] = TT[k];
    }
    
    // 配列m_aEtabVt, m_aIb0, m_aLnIb0をセットする
    for(int k=0; k<m_gc.K; k++){
	double Tkel_k = TT[k]+273.15;
	double kT_k = 1.380649E-23*Tkel_k;
	m_aEtabVt[k] = 1.1*kT_k/1.602176634E-19;

	// m_aIb0[k] = getIb0(k);
	m_aIb0[k] = m_gc.IscSTC*(1.0+m_gc.alphaSc*(m_aT[k]-25.0))/( exp(0.6/m_aEtabVt[k])-1.0 );
	m_aLnIb0[k] = log(m_aIb0[k]);
    }

    return m_amInit.calcPCS(GG, TT);
}

/** m_aG[]、m_aT[] と m_aEtabVt[], m_aIb0[], m_aLnIb0[] をセットする  ****************
 *  StringGroup::getE用、   速度重視
 *  AdaptiveModel m_amInitのパラメータも対応してセットする
 *  出力： m_aG[], m_aT[], m_aEtabVt[], m_aIb0[], m_aLnIb0[]
 *  関数呼び出し： AdaptiveModel::calcPCS
 *  戻り値： AdaptiveModel::calcのE値     *******************************************/
double Substring::setVarsGetE(){
    // for(int k=0; k<m_gc.K; k++){
    // 	m_aG[k] = GG[k];
    // 	m_aT[k] = TT[k];
    // }
    
    // // 配列m_aEtabVt, m_aIb0, m_aLnIb0をセットする
    // for(int k=0; k<m_gc.K; k++){
    // 	double Tkel_k = TT[k]+273.15;
    // 	double kT_k = 1.380649E-23*Tkel_k;
    // 	m_aEtabVt[k] = 1.1*kT_k/1.602176634E-19;

    // 	// m_aIb0[k] = getIb0(k);
    // 	m_aIb0[k] = m_gc.IscSTC*(1.0+m_gc.alphaSc*(m_aT[k]-25.0))/( exp(0.6/m_aEtabVt[k])-1.0 );
    // 	m_aLnIb0[k] = log(m_aIb0[k]);
    // }

    // return m_amInit.calcPCS(GG, TT);
    return m_amInit.calcPCSGetE();
}

/***********************************************************************/
void Substring::setGT(double G, double T){
    m_G = G;
    m_T = T;
    m_amInit.setGT(G, T);
    for(int i=0; i<m_vam.size(); i++) m_vam[i].setGT(G, T);

    m_Tkel = m_T+273.15;
    m_kT = 1.380649E-23*m_Tkel;
    m_etabVt = 1.1*m_kT/1.602176634E-19;
    m_Ib0 = m_gc.IscSTC*(1.0+m_gc.alphaSc*(m_T-25.0))/( exp(0.6/(m_etabVt))-1.0 );
    m_lnIb0 = log(m_Ib0);

    m_ImaxMax = m_amInit.getImax();
    m_ImaxMin = m_amInit.getImax();
    m_VmaxInit = m_amInit.getVmax();
    m_VmaxSum = m_VmaxInit*m_gc.Ns_sub;
    for(int i=0; i<m_vam.size(); i++){
	double Imax {m_vam[i].getImax()};
	if(Imax > m_ImaxMax) m_ImaxMax = Imax;
	if(Imax < m_ImaxMin) m_ImaxMin = Imax;
	
	double Vmax {m_vam[i].getVmax()};
	m_VmaxSum -= m_VmaxInit*m_vamNum[i];
	m_VmaxSum += Vmax*m_vamNum[i];
    }
    
    // m_nVmin = m_gc.Ns_sub*m_amInit.getVmin();
    m_VbMin = -m_etabVt*std::log(m_IbMax/m_Ib0+1.0);
    m_VbypassDiodeOn = -m_etabVt*log(1.0E-3/m_Ib0 +1);  // 1mA以上ならオンとみなす    
}

/****** Suppose refreshment after set5params of adaptiveModel ******/
void Substring::refreshMaxMin(){
    m_amInit.refreshImaxVmax();
    for(int i=0; i<m_vam.size(); i++) m_vam[i].refreshImaxVmax();
    
    m_ImaxMax = m_amInit.getImax();
    m_ImaxMin = m_amInit.getImax();
    m_VmaxInit = m_amInit.getVmax();
    m_VmaxSum = m_VmaxInit*m_gc.Ns_sub;
    for(int i=0; i<m_vam.size(); i++){
	double Imax {m_vam[i].getImax()};
	if(Imax > m_ImaxMax) m_ImaxMax = Imax;
	if(Imax < m_ImaxMin) m_ImaxMin = Imax;
	
	double Vmax {m_vam[i].getVmax()};
	m_VmaxSum -= m_VmaxInit*m_vamNum[i];
	m_VmaxSum += Vmax*m_vamNum[i];
    }
}

/** setVarsのあと、Preparation of getV(Ik), getI(Vk) for k   ************************
 *  入力： m_amInit, m_aIb0[k], m_aEtabVt[k], m_aLnIb0[k]
 *  出力： m_BypassDiodeOn[k]      *************************************************/
void Substring::getVIjunbi(int k){
    m_aBypassDiodeOn[k] = false;
    m_amInit.getVIjunbi(k);

    m_Ib0 = m_aIb0[k];
    m_etabVt = m_aEtabVt[k];
    m_lnIb0 = m_aLnIb0[k];
}

/** After setting for getV(Ik) of k   ***********************************************
 *  入力： m_bypassDiodeOn
 *  出力： m_BypassDiodeOn[k]    ****************************************************/
void Substring::getVIafter(int k){
    m_aBypassDiodeOn[k] = m_bypassDiodeOn;
}

/********************************************************************************/
void Substring::addAM(AdaptiveModel &am, int amNum){
    m_vam.push_back(am);
    m_vamNum.push_back(amNum);
    // printf("ImaxMax=%g, ImaxMin=%g, VmaxSum=%g\n", m_ImaxMax, m_ImaxMin, m_VmaxSum);/////////
    double Imax {am.getImax()};
    if(Imax > m_ImaxMax) m_ImaxMax = Imax;
    if(Imax < m_ImaxMin) m_ImaxMin = Imax;
    double Vmax {am.getVmax()};
    m_VmaxSum -= m_VmaxInit*amNum;
    m_VmaxSum += Vmax*amNum;
    m_vamNumSum += amNum;
    // printf("ImaxMax=%g, ImaxMin=%g, VmaxSum=%g\n", m_ImaxMax, m_ImaxMin, m_VmaxSum);/////////
}

/********************************************************************************/
void Substring::addAM(double X[4], int amNum){
    AdaptiveModel am {static_cast<double>(X[0]), static_cast<double>(X[1]),
		      static_cast<double>(X[2]*X[0]),
		      static_cast<double>(X[3]*X[1]),
		      m_gc.alphaSc, m_gc.betaOc, m_gc.gammaMp, m_G, m_T, 1, m_gc};
    m_vam.push_back(am);
    m_vamNum.push_back(amNum);
    // printf("ImaxMax=%g, ImaxMin=%g, VmaxSum=%g\n", m_ImaxMax, m_ImaxMin, m_VmaxSum);/////////
    double Imax {am.getImax()};
    if(Imax > m_ImaxMax) m_ImaxMax = Imax;
    if(Imax < m_ImaxMin) m_ImaxMin = Imax;
    double Vmax {am.getVmax()};
    m_VmaxSum -= m_VmaxInit*amNum;
    m_VmaxSum += Vmax*amNum;
    m_vamNumSum += amNum;
    // printf("ImaxMax=%g, ImaxMin=%g, VmaxSum=%g\n", m_ImaxMax, m_ImaxMin, m_VmaxSum);/////////
}

/** AdaptiveModel m_amInitのSTCパラメータをセットする   ********************************
 *  関数呼び出し： AdaptiveModel::setSTCparams    *********************************/
void Substring::setAM_STCparams(double X[4]){
    // if(X[0] != m_Xpre[0] || X[1] != m_Xpre[1] || X[2] != m_Xpre[2] || X[3] != m_Xpre[3]){
	m_amInit.setSTCparams(X);
    // 	for(int i=0; i<4; i++) m_Xpre[i] = X[i];
    // }
}

/**     *************************************************************************/
void Substring::print() const{
    m_amInit.print();
    for(int k=0; k<GlobalConstants::K; k++){
	printf("Substring: k=%d, aG=%g, aT=%g, aIb0=%g, aEtabVt=%g, aLnIb0=%g\n",
	       k, m_aG[k], m_aT[k], m_aIb0[k], m_aEtabVt[k], m_aLnIb0[k]);
    }
}

/**     *************************************************************************/
FiveParamsPlus Substring::get5paramsPlus() const{
    FiveParamsPlus fpp;

    fpp.fiveParams = m_amInit.get5params();
    for(int k=0; k<GlobalConstants::K; k++){
	fpp.aBypassDiodeOn[k] = m_aBypassDiodeOn[k];
    }

    return fpp;
}

/**     *************************************************************************/
void FiveParamsPlus::print() const{
    for(int k=0; k<GlobalConstants::K; k++){
	printf("k=%d: Rs=%g, eta=%g, Rh=%g, I0=%g, Iph=%g, bypassDiodeOn=%d\n",
	       k, fiveParams.aRs[k],
	       fiveParams.aEta[k], fiveParams.aRh[k],
	       fiveParams.aI0[k], fiveParams.aIph[k], aBypassDiodeOn[k]);
    }
}
