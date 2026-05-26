/***************************************************************************
 *         StringGroup.hh  2024.1.21--
 ***************************************************************************/
#if !defined(__STRING_GROUP_HH)
#define __STRING_GROUP_HH

#include "stringModel.hh"
// #include <string>
#include "nr3.h"
#include "mins.h"

class StringGroup{
private:
    GlobalConstants &m_gc;
    
    StringModel m_smInit;
    std::vector<int> m_vsmNum;
    int m_vsmNumSum {0};

    // smInitのもの
    double m_G;   // W/m^2
    double m_T;   // degC

    double m_Ib0;
    double m_Ib0sg;

    double m_Vmin;
    double m_Vmax;
    
    double m_Imin;
    double m_Imax;
    
    // for bisection search
    double m_Ifunc;

    double ELPSOsearch(double X[GlobalConstants::K], int s, int l);
    double getE(double X[GlobalConstants::K], int s, int l);
    
    double m_aG[GlobalConstants::K], m_aT[GlobalConstants::K];
    double m_aIsg[GlobalConstants::K], m_aVsg[GlobalConstants::K];

    double m_bunbo;
public:
    std::vector<StringModel> m_vsm; // パラメータが変更されたストリング 0,1,2,...
    // Constructor
    StringGroup(GlobalConstants &gcc, double GG, double TT);

    double getI(double V);
    double getV(double I);

    void getI(double Ik[GlobalConstants::K]);

    double getVmin() const{ return m_Vmin; };
    double getVmax() const{ return m_Vmax; };
    double getImin() const{ return m_Imin; };
    double getImax() const{ return m_Imax; };

    // for bisection search
    double func(double V);

    void addSM(StringModel &sm, int s);
    void setGT(double G, double T);
    void refreshMaxMin();

    // for getPmax
    double operator()(const double V);   // returns -P at V=V
    double getPmax(double &Vmp, double &Imp);

    void paramEstim(double G[], double T[], double I[], double V[], std::string Time[]);
};

// class XX {
//     int a = 3;
// public:
//     XX(){ return a; }

// };


#endif
