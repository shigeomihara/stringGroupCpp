/**********************************************************************
 *       stringModel.hh
 *********************************************************************/
#if !defined(__STRINGMODEL_HH)
#define __STRINGMODEL_HH

// #include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include "globalConstants.hh"
#include "substring.hh"

class StringModel{
private:
    GlobalConstants &m_gc;

    double m_G;   // W/m^2
    double m_T;   // degC

    Substring m_ssInit;
    std::vector<int> m_vssNum;
    int m_vssNumSum {0};
    
    double m_etabVt;
    double m_Ib0;  // バイパスダイオードの反転飽和電流

    double m_IscSub, m_Vmin, m_Vmax;
    double m_mIb099;
    double m_IscSub2;
    
    double m_Vfunc;   // for bisection search

    // 時系列としての測定G,T, PCS下の全サブストリングで一様
    // double aG[GlobalConstants::K], aT[GlobalConstants::K];
    double m_aVmin[GlobalConstants::K], m_aVmax[GlobalConstants::K],
	m_aIb0[GlobalConstants::K], m_aIscSub[GlobalConstants::K],
	m_aEtabVt[GlobalConstants::K]; 
    double m_aIk[GlobalConstants::K];  // 時点kに対するI, i.e. I_k
public:
    std::vector<Substring> m_vss;
    // Constructor
    StringModel(GlobalConstants &gcc, double GG, double TT);
    
    double getIb0() const{ return m_Ib0; };
    double getVmin() const{ return m_Vmin; };
    double getVmax() const{ return m_Vmax; };
    
    double getV(double I);
    double getI(double V);

    void getI(int l, double Vk[GlobalConstants::K], double Ik[GlobalConstants::K]);

    void getVjunbi(int k);
    void getIjunbi(int k);
    
    void getVjunbiVss(int k, int l);
    void getIjunbi(int k, int l);
    
    double func(double x);  // getIで使われる２分法用

    void addSS();
    void addSS(Substring &ss, int l);
    void setGT(double G, double T);
    void refreshMaxMin();

    double setVars(double GG[], double TT[]);
    double setVarsVss(double GG[], double TT[], int l);
    double setVarsVssGetE(int l);
    
    void setSS_STCparams(double X[GlobalConstants::K], int l);

    void set_aIk(double Ik[GlobalConstants::K]);
    void set_aIk_for(double V[GlobalConstants::K]);
    double get_aIk(int k) const{ return m_aIk[k]; }

    FiveParams get5params(int l) const;
    FiveParamsPlus get5paramsPlus(int l) const;

    void print(int l) const;
};

#endif
