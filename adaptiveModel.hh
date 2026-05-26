/*********************************************************************
 *      adaptiveModel.hh
 ********************************************************************/
#if !defined(_ADAPTIVE_MODEL_HH)
#define _ADAPTIVE_MODEL_HH

#if !defined(UBUNTU)
#include <cuda_runtime.h>
#endif
#include "globalConstants.hh"
#include <vector>

struct FiveParams{
    double aRs[GlobalConstants::K], aEta[GlobalConstants::K], aRh[GlobalConstants::K],
	aI0[GlobalConstants::K], aIph[GlobalConstants::K];

    void print() const;
};

class AdaptiveModel{
private:
    GlobalConstants &m_gc;
    
#if !defined(UBUNTU)
    // const unsigned int numBlocksRs {4};
    // const unsigned int numBlocksEta {4};
    // const unsigned int numBlocksRs {8};
    // const unsigned int numBlocksEta {8};
    const unsigned int m_numBlocksRs {32};
    const unsigned int m_numBlocksEta {32};
    const dim3 m_numBlocks {m_numBlocksRs, m_numBlocksEta};
    // const dim3 m_threadsPerBlock {32, 32};   // 32*32=10^10=1024
    const dim3 m_threadsPerBlock {16, 16};   // 32*32=10^10=1024x
    // const dim3 m_threadsPerBlock {4, 4};   // 32*32=10^10=1024

    // const int nRs {256};  // 各GPUスレッドが調べるRsの個数
    // const int nEta {256};  // 各GPUスレッドが調べるEtaの個数
    // const int nRs {128};  // 各GPUスレッドが調べるRsの個数
    // const int nEta {128};  // 各GPUスレッドが調べるEtaの個数
    // const int m_nRs {64};  // 各GPUスレッドが調べるRsの個数
    // const int m_nEta {64};  // 各GPUスレッドが調べるEtaの個数
    const int m_nRs {32};  // 各GPUスレッドが調べるRsの個数
    const int m_nEta {32};  // 各GPUスレッドが調べるEtaの個数
    // const int m_nRs {16};  // 各GPUスレッドが調べるRsの個数
    // const int m_nEta {16};  // 各GPUスレッドが調べるEtaの個数
    // const int nRs {8};  // 各GPUスレッドが調べるRsの個数
    // const int nEta {8};  // 各GPUスレッドが調べるEtaの個数

    const int m_numThreadsRs = m_numBlocksRs*m_threadsPerBlock.x;
    const int m_numThreadsEta = m_numBlocksEta*m_threadsPerBlock.y;
    const int m_numThreads = m_numThreadsRs*m_numThreadsEta;
    const int m_NRs = m_numThreadsRs*m_nRs;
    const int m_NEta = m_numThreadsEta*m_nEta;
#endif
#ifdef UBUNTU
    const int m_NRs = 5000;
    const int m_NEta = 5000;
#endif    
    
    // Rs, etaの最初の探索範囲
    const double m_Rs0init = 0.001, m_Rs1init = 10.0, m_eta0init = 1.0, m_eta1init = 2.5;
    
    const double m_deltaOc = 0.0539, m_deltaMp = 0.0265;  // ND240-HA
    const double m_etad = -5.7E-4;  // /K

    // *** Notice *** If you modify the parameters higher than here,
    // you must remove the file 'searchHistoryAM.dat'.
    
    // このクラスのメインの状態変数、コンストラクタで設定される
    double m_IscSTC, m_VocSTC, m_ImpSTC, m_VmpSTC, m_alphaSc, m_betaOc, m_gammaMp, m_G, m_T;
    int m_Ns;

    // 関数setVarsで設定される
    double m_Rs0=m_Rs0init, m_Rs1=m_Rs1init, m_eta0=m_eta0init, m_eta1=m_eta1init;

    double m_GpGSTC;
    double m_TmTSTC;
    double m_lnGpGSTC;
    double m_PmaxSTC;
  
    double m_Tkel;
    double m_NsVt;
    double m_W01;
    double m_W02;

    double m_Voc, m_Isc, m_Vmp, m_Imp;  // G,Tに対する値
    // 関数setVarsで設定される（ここまで）

    // 関数calcA0B0で設定される
    double m_A0, m_B0, m_eta00, m_absA0, m_absB0;  // eta00はA0, B0の計算に使われたetaの値、1.0, 1,1, 1,2, 1,3のいずれか

    // Five parameters of single diode model
    double m_Rs, m_eta, m_Rh, m_I0, m_Iph;

    // コンストラクタ、calcPCS()から呼び出される
    void setVars();
    void calcA0B0();
    void setVarsGetE();
  
    double calcRh(double Rs, double eta) const;
    double calcI0(double Rs, double eta, double Rh) const;
    double calcIph(double Rs, double eta, double Rh) const;

    bool sameCalc();  // Exist the same calc in m_vsr already?

    // PCSデータによるパラメータ推定用、calcPCSでセットされる
    double m_aG[GlobalConstants::K], m_aT[GlobalConstants::K];
    double m_aNsVt[GlobalConstants::K];
    double m_aIsc[GlobalConstants::K], m_aVoc[GlobalConstants::K],
	m_aImp[GlobalConstants::K], m_aVmp[GlobalConstants::K]; 
    double m_aRs[GlobalConstants::K], m_aEta[GlobalConstants::K], m_aRh[GlobalConstants::K],
	m_aI0[GlobalConstants::K], m_aIph[GlobalConstants::K];
    double m_aLnGpGSTC[GlobalConstants::K];
    double m_aW01[GlobalConstants::K], m_aW02[GlobalConstants::K];
    
    double SearchNoGPU(double kRs, double dEta);

    const double m_Vmin {-100.0};   // the least voltage of this cell (V)
    const double m_Imin {-1000.0};   // the least current of this cell (A)
    double m_Imax;
    double m_Vmax;

public:
    // Constructors G (W/m^2), T (degC)
    AdaptiveModel(double IscSTC, double VocSTC, double ImpSTC, double VmpSTC, double alphaSc,
		  double betaOc, double gammaMp, double G, double T, int Ns, GlobalConstants gc);
    AdaptiveModel(GlobalConstants &gc, double G, double T);
    AdaptiveModel(GlobalConstants gc);  // Constructor doing nothing

    double calc();
    // double calcRh(double Rs, double eta);//////////////////////t

    double calcNoGPU();

    double getRs() const{ return m_Rs; }
    double getEta() const{ return m_eta; }
    double getRh() const{ return m_Rh; }
    double getI0() const{ return m_I0; }
    double getIph() const{ return m_Iph; }

    double getV(double I) const;
    double getI(double V) const;

    double getIsc() const{ return m_Isc; }

    // PCSデータによるパラメータ推定用
    double calcPCS(double G[], double T[]);
    double calcPCSGetE();
    
    double getIsc(int k) const{ return m_aIsc[k]; }
    
    void getVIjunbi(int k);

    void set5params(double Rss, double Etaa, double Rhh, double I00, double Iphh);
    void set5params(double Rss, double Etaa, double Rhh, double I00, double Iphh, double NsVtt, int k);
    void setGT(double G, double T);
    void refreshImaxVmax();

    void setSTCparams(double X[4]);

    FiveParams get5params() const;

    void print() const;
    
    double getVmax() const{ return m_Vmax; }
    double getImax() const{ return m_Imax; }
    double getVmin() const{ return m_Vmin; }
    double getImin() const{ return m_Imin; }
};

// Kernel
#if !defined(UBUNTU)
__global__ void Search(double *d_Rs0nth, double *d_eta0nth,
		       double *d_RsMin,  double *d_RhMin, double *d_etaMin,double *d_Emin,
		       double Isc, double Voc, double Imp, double Vmp,
		       double alphaSc, double betaOc, double gammaMp,
		       double NsVt, double I01d, double I02d,
		       double etad, double TInv,
		       double kRs, double dEta, int nRs, int nEta,
		       double absA0, double absB0);
#endif
#endif
