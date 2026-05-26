/*******************************************************************
 *            Substring.hh
 ******************************************************************/
#if !defined(__SUBSTRING_HH)
#define __SUBSTRING_HH

#include <iostream>
#if !defined(UBUNTU)
#include <cuda/std/limits>  // numeric_limits
#endif
#include <cmath>  // std::log
#include <vector>
#include "globalConstants.hh"
#include "adaptiveModel.hh"
#include "bisection.hh"

struct FiveParamsPlus{
    FiveParams fiveParams;
    bool aBypassDiodeOn[GlobalConstants::K];

    void print() const;
};

class Substring{
private:
    GlobalConstants &m_gc;
    AdaptiveModel m_amInit;
    std::vector<int> m_vamNum;
    int m_vamNumSum {0};
    
    // コンストラクタでセットされる
    double m_G;   // W/m^2
    double m_T;   // degC
    double m_Tkel;
    double m_kT;
    double m_etabVt;
    double m_Ib0;
    double m_lnIb0;

    bool m_bypassDiodeOn {false};
    double m_VbypassDiodeOn;
    
    double m_Ifunc;   // func の bisection search 用
    double m_Vfunc2;   // func2 の bisection search 用

    double m_aG[GlobalConstants::K], m_aT[GlobalConstants::K];
    double m_aEtabVt[GlobalConstants::K], m_aIb0[GlobalConstants::K],
	m_aLnIb0[GlobalConstants::K];
    
    bool m_aBypassDiodeOn[GlobalConstants::K];

    double m_Xpre[4];

    double m_VmaxSum;
    double m_VmaxInit;
    double m_ImaxMax;
    double m_ImaxMin;
    double m_nVmin;
    const double m_IbMax {1000.0};
    double m_VbMin;
public:
    std::vector<AdaptiveModel> m_vam;
    // Constructor
    Substring(GlobalConstants &gcc, double GG, double TT);
    
    Substring(GlobalConstants gcc, double GG, double TT,
	      double IscSTC, double VocSTC, double ImpSTC, double VmpSTC);

    double getI(double V);
    double getV(double I);

    double func(double x) const;  // getVで使われる２分法用
    double func2(double x) const;  // getIで使われる２分法用

    bool getBypassDiodeOn() const{ return m_bypassDiodeOn; }
    
    double getIb0() const{ return m_Ib0; }
    double getEtabVt() const{ return m_etabVt; }
    double getIsc() const{ return m_amInit.getIsc(); }

    void setGT(double G, double T);
    void refreshMaxMin();

    double setVars(double G[], double T[]);
    double setVarsGetE();
    
    void getVIjunbi(int k);
    void getVIafter(int k);
    
    double getIb0(int k) const{ return m_aIb0[k]; };
    double getEtabVt(int k) const{ return m_aEtabVt[k]; }
    double getIsc(int k) const{ return m_amInit.getIsc(k); }
    
    bool getBypassDiodeOn(int k) const{ return m_aBypassDiodeOn[k]; }

    void addAM(AdaptiveModel &am, int amNum);
    void addAM(double X[4], int amNum);
    void setAM_STCparams(double X[4]);

    FiveParams get5params() const{ return m_amInit.get5params(); }
    FiveParamsPlus get5paramsPlus() const;

    void print() const;
};

#endif
