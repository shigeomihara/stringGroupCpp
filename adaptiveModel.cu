// -*- C++ -*-
/******************************************************************************
 *      adaptiveModel.cu
 *****************************************************************************/
#include "adaptiveModel.hh"
#include "LambertWhara.hh"
#include "globalConstants.hh"
#include <iostream>
#include <cmath>
#include <cstdio>
#include <limits>
#if !defined(UBUNTU)
#include <cuda_runtime.h>
#endif
#include <iomanip>  // calcPCSparams
#include <stdexcept>
#include <string>

/** Constructor G (W/m^2), T (degC)    **********************************************
 *  出力：IscSTC, VocSTC, ImpSTC, VmpSTC, alphaSc, betaOc, gammaMp, G, T, Ns
 *  関数呼び出し：setVars, calcA0B0    *********************************************/
AdaptiveModel::AdaptiveModel(double IscSTCC, double VocSTCC, double ImpSTCC, double VmpSTCC,
			     double alphaScc, double betaOcc, double gammaMpp,
			     double GG, double TT, int Nss, GlobalConstants gc)
    :m_IscSTC{IscSTCC}, m_VocSTC{VocSTCC}, m_ImpSTC{ImpSTCC}, m_VmpSTC{VmpSTCC},
     m_alphaSc{alphaScc}, m_betaOc{betaOcc}, m_gammaMp{gammaMpp}, m_G{GG}, m_T{TT},
     m_Ns{Nss}, m_gc{gc}{

    setVars();
    calcA0B0();
    calc();
    m_Imax = getI(m_Vmin);
    m_Vmax = getV(m_Imin);
}

/** GlobalConstantsを与えるコンストラクタ  ******************************************
 *  出力：IscSTC, VocSTC, ImpSTC, VmpSTC, alphaSc, betaOc, gammaMp, G, T, Ns
 *  関数呼び出し：setVars, calcA0B0     ********************************************/
AdaptiveModel::AdaptiveModel(GlobalConstants &gc, double GG, double TT)
    :m_IscSTC{gc.IscSTC}, m_VocSTC{gc.VocSTC/gc.Ns},
     m_ImpSTC{gc.ImpSTC}, m_VmpSTC{gc.VmpSTC/gc.Ns},
     m_alphaSc{gc.alphaSc}, m_betaOc{gc.betaOc}, m_gammaMp{gc.gammaMp},
     m_G{GG}, m_T{TT}, m_Ns{1}, m_gc{gc}{  // modified m_Ns{gc.Ns_sub}, Mar.4,2026

    setVars();
    calcA0B0();
    calc();    // added Mar.4, 2026
    m_Imax = getI(m_Vmin);
    m_Vmax = getV(m_Imin);
    // testVars(33.4);////////////////////////
}

/********* Empty Constructor, 9 Apr.2026  *************************/
AdaptiveModel::AdaptiveModel(GlobalConstants gc)
    :m_IscSTC{gc.IscSTC}, m_VocSTC{gc.VocSTC/gc.Ns},
     m_ImpSTC{gc.ImpSTC}, m_VmpSTC{gc.VmpSTC/gc.Ns},
     m_alphaSc{gc.alphaSc}, m_betaOc{gc.betaOc}, m_gammaMp{gc.gammaMp},
     m_Ns{1}, m_gc{gc}{  // modified m_Ns{gc.Ns_sub}, Mar.4,2026
}

/** クラスの変数たちを設定する。コンストラクタから呼び出される。  *******************
 *  入力：Rs0init, Rs1init, eta0init, eta1init, G, T, VocSTC, IscSTC, VmpSTC, ImpSTC, Ns,
 *        etad, alphaSC, betaOc, gammaMp
 *  出力：左辺  *********************************************************************/
void AdaptiveModel::setVars(){
    // printf("AdaptiveModel::setVars()$$ NRs=%d, NEta=%d\n", NRs, NEta);//////////////t
    
    m_Rs0 = m_Rs0init;
    m_Rs1 = m_Rs1init;
    m_eta0 = m_eta0init;
    m_eta1 = m_eta1init;

    m_GpGSTC = m_G/1000.0;
    m_lnGpGSTC = std::log(m_GpGSTC);
    m_TmTSTC = m_T - 25.0;
    m_PmaxSTC = m_ImpSTC*m_VmpSTC;
    m_Tkel = m_T+273.15;
    double kT = 1.380649E-23*m_Tkel;
    m_NsVt = m_Ns*kT/1.602176634E-19;

    double Eg = 1.8516E-19 -1.12488E-22*m_Tkel*m_Tkel/(m_Tkel+1108.0);
    double Egd = -1.12488E-22*(m_Tkel*m_Tkel+2.0*m_Tkel*1108.0)/( (m_Tkel+1108.0)*(m_Tkel+1108.0) );
    m_W01 = 3.0/m_Tkel -Egd/kT +Eg/(kT*m_Tkel);
    m_W02 = (Eg/kT -3.0*std::log(m_Tkel))*m_etad;

    m_Isc = m_IscSTC*m_GpGSTC*(1.0+m_alphaSc*m_TmTSTC);
    m_Voc = m_VocSTC/(1.0-m_deltaOc*m_lnGpGSTC)*(1.0+m_betaOc*m_TmTSTC);
    m_Vmp = m_VmpSTC*( 1.0/(1.0-m_deltaOc*m_lnGpGSTC) +m_betaOc*m_TmTSTC );
    double Pmax = m_PmaxSTC*m_GpGSTC/(1.0-m_deltaMp*m_lnGpGSTC)*(1.0+m_gammaMp*m_TmTSTC);
    m_Imp = Pmax/m_Vmp;
}    

/** A0, ,absA0, B0, absB0, eta00を計算する。コンストラクタから実行される  2024.1.8--
 *  入力：Voc, Isc, Vmp, Imp, NsVt, alphaSc, betaOc, gammaMp, Tkel, etad, W01, W02
 *  出力：A0, ,absA0, B0, absB0, eta00  *********************************************/
void AdaptiveModel::calcA0B0(){
    double eta = 1.0;
    
    for(int count=0; count<3; count++){
	double x = m_Voc/(eta*m_NsVt);
	double z = m_Vmp/(eta*m_NsVt);

	double emx = std::exp(-x);
	double ezmx = std::exp(z-x);
	
	double Rh = (m_Vmp+(m_Voc-m_Vmp)*emx-m_Voc*ezmx) /(m_Isc-m_Imp+m_Imp*emx-m_Isc*ezmx);
	
#ifdef UBUNTU
	if(std::isnan(Rh)){   // Rh<0は許している
	    eta += 0.1;
	    continue;
	}
#else
	if(isnan(Rh)){   // Rh<0は許している
	    eta += 0.1;
	    continue;
	}
#endif
	
	double I01 = m_Isc-m_Voc/Rh;
	
	double one_pXm1 = emx/(1.0-emx);
	double ZpXm1 = ezmx/(1.0-emx);

	double I0 = I01*one_pXm1;

	this->m_A0 = -z*I01*ZpXm1+m_Imp-m_Vmp/Rh;
	this->m_absA0 = std::abs(m_A0);
	

	if(m_absA0<std::sqrt(std::numeric_limits<double>::min()) || m_absA0>std::sqrt(std::numeric_limits<double>::max())){
	    eta += 0.1;
	    continue;
	}

	double W0 = m_W01/eta +m_W02/(eta*eta);
	double I0d = W0*I0;

	// Rhdの計算
	double W1 = Rh*( m_Voc/(eta*m_NsVt)*( m_etad/eta -m_betaOc +1.0/m_Tkel ) -W0 );  // Tkel modified 2024.1.20 !!!!!
	double W2 = (Rh*m_alphaSc)*m_Isc +Rh*I0d -m_betaOc*m_Voc;  // +Rh*I0d modified 2024.1.20 !!!!

	double Rhd = -Rh/m_Voc*( I01/(1.0-emx)*W1 +W2 );

	this->m_B0 = I01*ZpXm1*( m_Vmp/(eta*m_NsVt)*(m_etad/eta-m_gammaMp+1.0/m_Tkel) -W0)
	    +m_Vmp/(Rh*Rh)*Rhd +m_alphaSc*m_Isc+I0d-m_gammaMp*m_Vmp/Rh;   // Tkel modified 2024.1.21 !!!!!
	this->m_absB0 = std::abs(m_B0);

	if(m_absB0<std::sqrt(std::numeric_limits<double>::min()) || m_absB0>std::sqrt(std::numeric_limits<double>::max())){
	    eta += 0.1;
	    continue;
	}

	this->m_eta00 = eta;

	return;  // ここまで来たら問題ないので、return
    }

    std::cout << "Error in adaptiveModel.calcA0B0" << std::endl;
    exit(0);
}

/** 現在のVoc, Isc, Vmp, Impに対するRs, eta, Rh, I0, Iphをメッシュサーチして設定する 
 *  入力：numThreads, numThreadsRs, numThreadsEta, Rs0, Rs1, NRs, eta0, eta1, NEta
 *  出力：Rs0, Rs1, eta0, eta1, Rs, eta, Rh, I0, Iph
 *  関数呼び出し：Search, calcI0, calcIph
 *  戻り値：Emin (==INFINITYのときは、Rs, etaなどが設定されていないので注意）   *****/
double AdaptiveModel::calc(){
    SearchResult sr0 {m_Isc, m_Voc, m_Imp, m_Vmp, m_Tkel, m_absA0, m_absB0,
		     m_Rs, m_eta, m_Rh, m_I0, m_Iph, 0.0};
    int iSr {m_gc.sh.haveSearchResult(sr0)};
    if(iSr >= 0){
	m_Rs = m_gc.sh.m_vsr[iSr].Rs;
	m_eta = m_gc.sh.m_vsr[iSr].eta;
	m_Rh = m_gc.sh.m_vsr[iSr].Rh;
	m_I0 = m_gc.sh.m_vsr[iSr].I0;
	m_Iph = m_gc.sh.m_vsr[iSr].Iph;
	// printf("calc: haveSearchResult\n");//////////////////////////
	return m_gc.sh.m_vsr[iSr].Emin;
    }
    
    double Emin {std::numeric_limits<double>::infinity()};
    
    if(!m_gc.haveGPU) Emin = calcNoGPU();  // GPUを使わないとき

#if !defined(UBUNTU)
    cudaError_t err = cudaSuccess;
    size_t sizeAll {m_numThreads*sizeof(double)};  // スレッドの総数に対応するdoubleメモリサイズ
    size_t sizeRs {m_numThreadsRs*sizeof(double)};
    size_t sizeEta {m_numThreadsEta*sizeof(double)};

    // 各スレッドのRs探索範囲の下限である配列Rs0nthの準備
    double *Rs0nth = (double *)malloc(sizeRs);
    if(Rs0nth==NULL){
	fprintf(stderr, "Failed to allocate host vector Rs0nth\n");
	exit(EXIT_FAILURE);
    }
  
    // 各スレッドのEta探索範囲の下限である配列eta0nthの準備
    double *eta0nth = (double *)malloc(sizeEta);
    if(eta0nth==NULL){
	fprintf(stderr, "Failed to allocate host vector eta0nth\n");
	exit(EXIT_FAILURE);
    }
  
    // 各スレッドが探索するRsの範囲の下限
    double *d_Rs0nth = NULL;
    err = cudaMalloc((void **)&d_Rs0nth, sizeRs);
    if (err != cudaSuccess){
	fprintf(stderr, "Failed to allocate device vector Rs0nth (error code %s)!\n", cudaGetErrorString(err));
	exit(EXIT_FAILURE);
    }

    // 各スレッドが探索するEtaの範囲の下限
    double *d_eta0nth = NULL;
    err = cudaMalloc((void **)&d_eta0nth, sizeEta);
    if (err != cudaSuccess){
	fprintf(stderr, "Failed to allocate device vector Eta0nth (error code %s)!\n", cudaGetErrorString(err));
	exit(EXIT_FAILURE);
    }

    // 各スレッドが探索した最小のEに対するRsの値
    double *d_RsMin = NULL;
    err = cudaMalloc((void **)&d_RsMin, sizeAll);
    if (err != cudaSuccess){
	fprintf(stderr, "Failed to allocate device vector RsMin (error code %s)!\n", cudaGetErrorString(err));
	exit(EXIT_FAILURE);
    }
  
    // 各スレッドが探索した最小のEに対するRhの値、 Rh<0のことがあるので
    double *d_RhMin = NULL;
    err = cudaMalloc((void **)&d_RhMin, sizeAll);
    if (err != cudaSuccess){
	fprintf(stderr, "Failed to allocate device vector RhMin (error code %s)!\n", cudaGetErrorString(err));
	exit(EXIT_FAILURE);
    }
  
    // 各スレッドが探索した最小のEに対するetaの値
    double *d_etaMin = NULL;
    err = cudaMalloc((void **)&d_etaMin, sizeAll);
    if (err != cudaSuccess){
	fprintf(stderr, "Failed to allocate device vector etaMin (error code %s)!\n", cudaGetErrorString(err));
	exit(EXIT_FAILURE);
    }
  
    // 各スレッドが探索した最小のEの値
    double *d_Emin = NULL;
    err = cudaMalloc((void **)&d_Emin, sizeAll);
    if (err != cudaSuccess){
	fprintf(stderr, "Failed to allocate device vector Emin (error code %s)!\n", cudaGetErrorString(err));
	exit(EXIT_FAILURE);
    }

    double *RsMinTh = (double *)malloc(sizeAll);  // コピー先の配列
    if(RsMinTh==NULL){
	fprintf(stderr, "Failed to allocate host vector RsMinTh\n");
	exit(EXIT_FAILURE);
    }
	
    double *RhMinTh = (double *)malloc(sizeAll);  // コピー先の配列////////////////////t
    if(RhMinTh==NULL){
	fprintf(stderr, "Failed to allocate host vector RhMinTh\n");
	exit(EXIT_FAILURE);
    }
  
    double *etaMinTh = (double *)malloc(sizeAll);  // コピー先の配列
    if(etaMinTh==NULL){
	fprintf(stderr, "Failed to allocate host vector etaMinTh\n");
	exit(EXIT_FAILURE);
    }
	
    double *EminTh = (double *)malloc(sizeAll);  // コピー先の配列
    if(EminTh==NULL){
	fprintf(stderr, "Failed to allocate host vector EminTh\n");
	exit(EXIT_FAILURE);
    }
  
    for(int count=0; count<3; count++){  // 探索範囲変更のループ
    
	double kRs = std::exp( std::log(m_Rs1/m_Rs0) /(m_NRs-1) );
	double dEta = (m_eta1-m_eta0)/(m_NEta-1);

	Rs0nth[0] = m_Rs0;
	double knRs = std::exp( std::log(m_Rs1/m_Rs0)*m_nRs/(m_NRs-1) );

	for(int n=1; n<m_numThreadsRs; n++){
	    Rs0nth[n] = Rs0nth[n-1]*knRs;
	}
	
	eta0nth[0] = m_eta0;
	double dnEta = dEta*m_nEta;

	for(int n=1; n<m_numThreadsEta; n++){
	    eta0nth[n] = eta0nth[n-1]+dnEta;
	}

	// Rs探索範囲の下限である配列Rs0nthのGPUへのコピー
	err = cudaMemcpy(d_Rs0nth, Rs0nth, sizeRs, cudaMemcpyHostToDevice);
	if (err != cudaSuccess){
	    fprintf(stderr, "Failed to copy vector Rs0nth from host to device (error code %s)!\n", cudaGetErrorString(err));
	    exit(EXIT_FAILURE);
	}

	// Eta探索範囲の下限である配列Eta0nthのGPUへのコピー
	err = cudaMemcpy(d_eta0nth, eta0nth, sizeEta, cudaMemcpyHostToDevice);
	if (err != cudaSuccess){
	    fprintf(stderr, "Failed to copy vector eta0nth from host to device (error code %s)!\n", cudaGetErrorString(err));
	    exit(EXIT_FAILURE);
	}

	// GPUカーネルの実行
	Search<<<m_numBlocks, m_threadsPerBlock>>>(d_Rs0nth, d_eta0nth,
						   d_RsMin, d_RhMin, d_etaMin, d_Emin,
						   m_Isc, m_Voc, m_Imp, m_Vmp,
						   m_alphaSc, m_betaOc, m_gammaMp,
						   m_NsVt, m_W01, m_W02, m_etad, 1.0/m_Tkel,
						   kRs, dEta, m_nRs, m_nEta, m_absA0, m_absB0);
	err = cudaGetLastError();

	if (err != cudaSuccess){
	    fprintf(stderr, "Failed to launch CUDA kernel (error code %s)!\n", cudaGetErrorString(err));
	    exit(EXIT_FAILURE);
	}else{
	    // fprintf(stderr, "# Succeeded to launch CUDA kernel (success code: %s)!\n", cudaGetErrorString(err));/////////////////t
	}    

	// 探索されたRsの値のホストへのコピー
	err = cudaMemcpy(RsMinTh, d_RsMin, sizeAll, cudaMemcpyDeviceToHost);
	if (err != cudaSuccess){
	    fprintf(stderr, "Failed to copy vector RsMin from device to host (error code %s)!\n", cudaGetErrorString(err));
	    exit(EXIT_FAILURE);
	}
  
	// 探索されたRhの値のホストへのコピー  ///////////////////t
	err = cudaMemcpy(RhMinTh, d_RhMin, sizeAll, cudaMemcpyDeviceToHost);
	if (err != cudaSuccess){
	    fprintf(stderr, "Failed to copy vector RhMin from device to host (error code %s)!\n", cudaGetErrorString(err));
	    exit(EXIT_FAILURE);
	}
  
	// 探索されたetaの値のホストへのコピー
	err = cudaMemcpy(etaMinTh, d_etaMin, sizeAll, cudaMemcpyDeviceToHost);
	if (err != cudaSuccess){
	    fprintf(stderr, "Failed to copy vector etaMin from device to host (error code %s)!\n", cudaGetErrorString(err));
	    exit(EXIT_FAILURE);
	}
  
	// 探索されたEminの値のホストへのコピー
	err = cudaMemcpy(EminTh, d_Emin, sizeAll, cudaMemcpyDeviceToHost);
	if (err != cudaSuccess){
	    fprintf(stderr, "Failed to copy vector Emin from device to host (error code %s)!\n", cudaGetErrorString(err));
	    exit(EXIT_FAILURE);
	}

	// 最小値（最適値）の決定
	Emin = std::numeric_limits<double>::infinity();
	for(int nth=0; nth<m_numThreads; nth++){
	    if(EminTh[nth] < Emin){
		// this->m_Rs = RsMinTh[nth];
		// this->m_eta = etaMinTh[nth];
		// this->m_Rh = RhMinTh[nth];//////////////////////t
		m_Rs = RsMinTh[nth];
		m_eta = etaMinTh[nth];
		m_Rh = RhMinTh[nth];//////////////////////t
		Emin = EminTh[nth];
	    }
	}

	// // 2024.3.20
	// // ひとつもRh>=0のような条件を満たすRs, etaがなかったら、探索を終わる
	// if(Emin == INFINITY){
	//     // std::cout << "AdaptiveModel::calc()$$ no suitable Rs and eta, Emin==INFINITY\n";///////////////////t
	//     break;
	// }
	
	// 最適値が範囲の端から2番目以内だったら、探索範囲を取り直して再探索。
	if(this->m_Rs <= m_Rs0*pow(kRs, 1.5)){
	    double a = 1.05*m_Rs0/m_Rs1;
	    m_Rs0 *= a;
	    m_Rs1 *= a;
	    // printf("New smaller  Rs0=%g, Rs1=%g\n", Rs0, Rs1);////////////////t
	}else if(m_Rs1/pow(kRs, 1.5) <= this->m_Rs){
	    double b = 0.95*m_Rs1/m_Rs0;
	    m_Rs0 *= b;
	    m_Rs1 *= b;
	    // printf("New larger Rs0=%g, Rs1=%g\n", Rs0, Rs1);////////////////t
	    // }else if(eta <= eta0+1.5*d){
	    //   double eee =eta0;
	    //   eta0 *= 0.1;
	    //   eta1 = eee +d;
	    //   printf("New smaller  eta0=%g, eta1=%g\n", eta0, eta1);
	}else if(m_eta1-1.5*dEta <= this->m_eta){
	    double eee = m_eta0;
	    m_eta0 = m_eta1-dEta;
	    m_eta1 = 2.0*m_eta1 -eee -dEta;
	    // printf("New larger  eta0=%g, eta1=%g\n", eta0, eta1);///////////////t
	}else{  // 端っこでなかったら、終わる
	    break;
	}
    }
  
    // パラメータRh, I0, Iphの値設定
    // Rh = calcRh(this->Rs, this->eta);///////////////////t GPUで計算したRhを使うことにする。非常に微妙で、わずかの数値誤差で大きく変わることがある ex. G=770.4000244, T=57.41066742
    // m_I0 = calcI0(this->m_Rs, this->m_eta, m_Rh);
    // m_Iph = calcIph(this->m_Rs, this->m_eta, m_Rh);
    m_I0 = calcI0(m_Rs, m_eta, m_Rh);
    m_Iph = calcIph(m_Rs, m_eta, m_Rh);

    // mallocで確保したメモリの解放
    free(Rs0nth);
    free(eta0nth);
    free(RsMinTh);
    free(RhMinTh);/////////////////t
    free(etaMinTh);
    free(EminTh);

    err = cudaFree(d_Rs0nth);
    if (err != cudaSuccess){
	fprintf(stderr, "Failed to free device vector Rs0nth (error code %s)!\n", cudaGetErrorString(err));
	exit(EXIT_FAILURE);
    }
  
    err = cudaFree(d_eta0nth);
    if (err != cudaSuccess){
	fprintf(stderr, "Failed to free device vector eta0nth (error code %s)!\n", cudaGetErrorString(err));
	exit(EXIT_FAILURE);
    }
  
    err = cudaFree(d_RsMin);
    if (err != cudaSuccess){
	fprintf(stderr, "Failed to free device vector RsMin (error code %s)!\n", cudaGetErrorString(err));
	exit(EXIT_FAILURE);
    }
  
    err = cudaFree(d_RhMin);/////////////////////t
    if (err != cudaSuccess){
	fprintf(stderr, "Failed to free device vector RhMin (error code %s)!\n", cudaGetErrorString(err));
	exit(EXIT_FAILURE);
    }
  
    err = cudaFree(d_etaMin);
    if (err != cudaSuccess){
	fprintf(stderr, "Failed to free device vector etaMin (error code %s)!\n", cudaGetErrorString(err));
	exit(EXIT_FAILURE);
    }
  
    err = cudaFree(d_Emin); 
    if (err != cudaSuccess){
	fprintf(stderr, "Failed to free device vector Emin (error code %s)!\n", cudaGetErrorString(err));
	exit(EXIT_FAILURE);
    }
#endif
    
    SearchResult sr {m_Isc, m_Voc, m_Imp, m_Vmp, m_Tkel, m_absA0, m_absB0,
		     m_Rs, m_eta, m_Rh, m_I0, m_Iph, Emin};
    if(m_gc.sh.haveSearchResult(sr) < 0){
	// printf("new search Result, m_vsr.size()=%zd\n", m_gc.sh.m_vsr.size());////////////////
	m_gc.sh.m_vsr.push_back(sr);  // save a search result in the history
	printf("Emin=%g, Rs=%g, eta=%g, Rh=%g, I0=%g, Iph=%g, m_vsr.size()=%zd\n", Emin, m_Rs, m_eta, m_Rh, m_I0, m_Iph, m_gc.sh.m_vsr.size());/////////////////////t
    }

    return Emin;
}

#if !defined(UBUNTU)
/**************************   Kernel definition   ***********************************/
__global__ void Search(double *d_Rs0nth, double *d_eta0nth,
		       double *d_RsMin, double *d_RhMin, double *d_etaMin, double *d_Emin,
		       double Isc, double Voc, double Imp, double Vmp,
		       double alphaSc, double betaOc, double gammaMp,
		       double NsVt, double W01, double W02,
		       double etad, double TInv,
		       double kRs, double dEta, int nRs, int nEta,
		       double absA0, double absB0){
    int idx_x = blockDim.x*blockIdx.x + threadIdx.x;
    int idx_y = blockDim.y*blockIdx.y + threadIdx.y;
    int idx = idx_x*(gridDim.y*blockDim.y) + idx_y;

    double VocMiVmp = Voc-Vmp;
    double IscMiImp = Isc-Imp;
  
    double RsMinVal = 0.0;
    double RhMinVal = 10.0;    /////////////////////t
    double etaMinVal = 1.0;  // 2024.2.25
    double EminVal = INFINITY;

    double Rs = d_Rs0nth[idx_x];  // Global Memroy内の配列からRs探索範囲の下限を読み取る
    double eta = d_eta0nth[idx_y];  // Global Memroy内の配列からRs探索範囲の下限を読み取る
    for(int i=0; i<nRs; i++){ // Rsのループ
	double RsIsc = Rs*Isc;
	double RsImp = Rs*Imp;
	double VmpPluRsImp = Vmp +RsImp;
	double RsImpMiVmp = RsImp-Vmp;
    
	for(int j=0; j<nEta; j++){ // etaのループ
	    
	    double etaNsVt = eta*NsVt;
	    double x =Voc/etaNsVt, y=RsIsc/etaNsVt, z=VmpPluRsImp/etaNsVt;

	    if(x==y) continue;
	    
	    // Rhの計算
	    double Rh;
	    double ZpXY = NAN;  // Z/(X-Y)もついでにここで計算しておく
	    double eyx, ezx, exy, ezy, exz, eyz;
	    
	    if(x>=y && x>=z){
		eyx = exp(y-x);
		ezx = exp(z-x);
		Rh = ( ( Vmp +VocMiVmp*eyx -Voc*ezx )
		       / ( IscMiImp +Imp*eyx -Isc*ezx ) ) - Rs;
		ZpXY = ezx/(1.0-eyx);
	    }else if(y >= z){
		exy = exp(x-y);
		ezy = exp(z-y);
		Rh = ( ( Vmp*exy +VocMiVmp -Voc*ezy )
		       / ( IscMiImp*exy +Imp -Isc*ezy ) ) - Rs;
		ZpXY = ezy/(exy-1.0);
	    }else{
		exz = exp(x-z);
		eyz = exp(y-z);
		Rh = ( ( Vmp*exz +VocMiVmp*eyz -Voc )
		       / ( IscMiImp*exz +Imp*eyz -Isc ) ) - Rs;
	    }

	    if( Rh<=0.0 || isnan(Rh) ) continue;

	    // I0の計算
	    // double I01 = (1.0 +Rs/Rh )*Isc -Voc/Rh;
	    double I01 = ( (Rh+Rs)*Isc -Voc )/Rh;
      
	    double oneXY, XoneXY, XpXY;
	    if(x>=y){
		double eyx1 = exp(y-x);
		double emx = exp(-x);
		oneXY = emx/(1.0-eyx1);
		XoneXY = (1.0-emx)/(1.0-eyx1);
		XpXY = 1.0/(1.0-eyx1);
	    }else{
		double exy1 = exp(x-y);
		double emy = exp(-y);
		oneXY = emy/(exy1-1.0);
		XoneXY = (exy1-emy)/(exy1-1.0);
		XpXY = exy1/(exy1-1.0);
	    }
      
	    double I0 = I01*oneXY;
	    if( !isfinite(I0) || I0<0.0) continue;

	    // Iphの計算
	    double Iph = I01*XoneXY +Voc/Rh;

	    if( !isfinite(Iph) || Iph<0.0) continue;

	    // absAの計算
	    double P0 = RsImpMiVmp/etaNsVt*I01;
	    double P1 = (RsImpMiVmp + Rh*Imp)/Rh;

	    double absA, lnZpXY, XYpZ;
	    if( !isnan(ZpXY) ){  // zが最大、ではなかったら
		absA = fabs(P0*ZpXY+P1);  // ZpXYは上のRhの計算のところで計算している、ｚが最大でない場合
	    }else{  // ｚが最大の場合
		if(x > y) lnZpXY = -log(exz-eyz);
		else lnZpXY = -log(eyz-exz);

		XYpZ = exz-eyz;

		double Atilder = lnZpXY +log( fabs(P0 +P1*XYpZ) +fabs(XYpZ) );
		absA = exp(Atilder)-1.0;
	    }
	    
	    if( !isfinite(absA) ) continue;
      
	    // I0dの計算
	    double W0 = (W01 +W02/eta)/eta;
	    double I0d = W0*I0;
	    double etadPerEta = etad/eta;

	    // Rsdの計算
	    double Rsd = (etadPerEta -alphaSc +TInv)*Rs +etaNsVt/Isc*(exp(-y)-1.0)*W0;

	    // Rhdの計算
	    double W1 = Rh*( Voc/etaNsVt*( etadPerEta -betaOc +TInv ) -W0 );
	    double W2 = (Rsd +(Rs+Rh)*alphaSc)*Isc +Rh*I0d -betaOc*Voc;  // +Rh*I0d modified 2024.1.21

	    double RhdPerRh = ( I01*XpXY*W1 +W2 )/(RsIsc-Voc);

	    // Iphdの計算
	    double Iphd = ( Rsd -Rs*RhdPerRh +(Rh+Rs)*alphaSc )*Isc/Rh;

	    // absBの計算
	    double xxx = -Imp*Rsd -gammaMp*Vmp;
	    double W3 = I01*( ( xxx +VmpPluRsImp*(etadPerEta +TInv) )/etaNsVt -W0 );
	    double W4 = ( xxx +VmpPluRsImp*RhdPerRh )/Rh +Iphd +I0d;

	    double absB;

	    if( !isnan(ZpXY) ){  // zが最大、ではなかったら
		absB = fabs(ZpXY*W3 +W4);   // ZpXYはｚが最大でないときは、Rhの計算のところで求めている
	    }else{  // ｚが最大のとき、このときlnZpXY, XYpZはAtilderのところで計算している
		double Btilder = lnZpXY +log( fabs(W3+W4*XYpZ)+fabs(XYpZ) );
		absB = exp(Btilder)-1.0;
	    }

	    if( !isfinite(absB) ) continue;
	    
	    // Eの計算
	    double E = absA/absA0 +absB/absB0;
	    if(E < EminVal){
		RsMinVal = Rs;
		etaMinVal = eta;
		EminVal = E;
		RhMinVal = Rh;
	    }
	    eta += dEta;
	}
	Rs *= kRs;
    }

    // 探索結果をGlobal Memory内の配列に出力する
    d_RsMin[idx] = RsMinVal;
    d_etaMin[idx] = etaMinVal;
    d_Emin[idx] = EminVal;
    d_RhMin[idx] = RhMinVal;
}
#endif

/** Rhを計算して返す, Rh==+infもあり得ることに注意（上のSearch()でNaNと<0だけ除外している)
 *  入力：Voc, Isc, Vmp, Imp, NsVt     **********************************************/
double AdaptiveModel::calcRh(double Rs, double eta) const{
    double Rh;
    double etaNsVt = eta*m_NsVt;
    double x =m_Voc/etaNsVt, y=Rs*m_Isc/etaNsVt, z=(m_Vmp+Rs*m_Imp)/etaNsVt;
    
    if(x>=y && x>=z){
	double eyx = std::exp(y-x);
	double ezx = std::exp(z-x);
	Rh = ( ( m_Vmp +(m_Voc-m_Vmp)*eyx -m_Voc*ezx )
	       / ( m_Isc-m_Imp +m_Imp*eyx -m_Isc*ezx ) ) - Rs;
    }else if(y >= z){
	double exy = std::exp(x-y);
	double ezy = std::exp(z-y);
	Rh = ( ( m_Vmp*exy +m_Voc-m_Vmp -m_Voc*ezy )
	       / ( (m_Isc-m_Imp)*exy +m_Imp -m_Isc*ezy ) ) - Rs;
    }else{
	double exz = std::exp(x-z);
	double eyz = std::exp(y-z);
	Rh = ( ( m_Vmp*exz +(m_Voc-m_Vmp)*eyz -m_Voc )
	       / ( (m_Isc-m_Imp)*exz +m_Imp*eyz -m_Isc ) ) - Rs;
    }

    return Rh;
}

/** I0を計算して返す   **************************************************************
 *  入力： Voc, Isc, NsVt  **********************************************************/
double AdaptiveModel::calcI0(double Rs, double eta, double Rh) const{
    double etaNsVt = eta*m_NsVt;
    double x =m_Voc/etaNsVt, y=Rs*m_Isc/etaNsVt;
  
    double I01 = ((Rh +Rs)*m_Isc -m_Voc)/Rh;
      
    double oneXY;
    if(x>=y){
	oneXY = std::exp(-x)/(1.0-std::exp(y-x));
    }else{
	oneXY = std::exp(-y)/(std::exp(x-y)-1.0);
    }
      
    double I0 = I01*oneXY;
    return I0;
}

/** Iphを計算して返す  **************************************************************
 *  入力：Voc, Isc, NsVt   **********************************************************/
double AdaptiveModel::calcIph(double Rs, double eta, double Rh) const{
    double etaNsVt = eta*m_NsVt;
    double x =m_Voc/etaNsVt, y=Rs*m_Isc/etaNsVt;
  
    double XoneXY;
    if(x>=y){
	XoneXY = (1.0-std::exp(-x))/(1.0-std::exp(y-x));
    }else{
	double exy = std::exp(x-y);
	XoneXY = (exy-std::exp(-y))/(exy-1.0);
    }
      
    double Iph = ( ((Rh +Rs )*m_Isc -m_Voc)*XoneXY +m_Voc )/Rh;
    return Iph;
}

/** I (A)に対する電圧V (V)を計算して返す   ******************************************
 *  入力：Rs, eta, Rh, I0, Iph, NsVt
 *  関数呼び出し： LambertWhara::LambertW0haraAexp     ******************************/
double AdaptiveModel::getV(double I) const{
    double bunbo = m_eta*m_NsVt;  // 桁落ち対策  ノート65-10頁
    
    if(std::isinf(m_Rh)){  // Rhが+infのとき 2024.3.28追加
	if(I>=(m_Iph+0.99*m_I0)) return -std::numeric_limits<double>::infinity();
	return bunbo*std::log1p((m_Iph - I)/m_I0) -m_Rs*I;
    }
	
    double a = m_Rh*m_I0/bunbo;
    double x = a*( (m_Iph - I)/m_I0 +1.0 );

    // printf("AM::getV a=%g, x=%g, Rh=%g, I0=%g, bunbo=%g\n", a, x, m_Rh, m_I0, bunbo);//////////t
    
    if(a<0.0){/////////////////////t
	// printf("AM::getV a=%g, Rh=%g, I0=%g, bunbo=%g\n", a, Rh, I0, bunbo);
	// exit(0);
	throw std::runtime_error("AM::getV a="+std::to_string(a)+ " Rh="+std::to_string(m_Rh)+" I0="+std::to_string(m_I0)+" bunbo="+std::to_string(bunbo));
    }

    return  bunbo*(x - utl::LambertW0haraAexp(a, x) ) - m_Rs*I;
}

/** V (V)に対する電流I (A)を計算して返す   ******************************************
 *  入力：Rs, eta, Rh, I0, Iph, NsVt   **********************************************/
double AdaptiveModel::getI(double V) const{
    if(std::isinf(m_Rh)){  // Rhが+infのとき 2024.3.28追加
	double bunbo = m_eta*m_NsVt;
	double aa = m_Rs*m_I0/bunbo;
	double xx = ((m_Iph+m_I0)*m_Rs+V)/bunbo;
	
	if(aa<0.0){/////////////////////t
	    printf("AM::getI a=%g, Rs=%g, I0=%g, bunbo=%g\n", aa, m_Rs, m_I0, bunbo);
	    exit(0);
	}
	
	return -bunbo/m_Rs*utl::LambertW0haraAexp(aa, xx) +m_Iph+m_I0;
    }

    double etaNsVt = m_eta*m_NsVt;
    double RsPlusRh = m_Rs+m_Rh;
    double a = m_Rs*m_Rh*m_I0/(etaNsVt*(RsPlusRh));
    double x = m_Rh*((m_Iph+m_I0)*m_Rs+V)/(etaNsVt*(RsPlusRh));
    
    if(a<0.0){/////////////////////t
	// if(true){//////////////////t
	printf("AM::getI a=%g, x=%g, Rs=%g, Rh=%g, I0=%g, etaNsVt=%g\n", a, x, m_Rs, m_Rh, m_I0, etaNsVt);
	exit(0);
    }

    return -etaNsVt/m_Rs*utl::LambertW0haraAexp(a, x) + ((m_Iph+m_I0)*m_Rh-V)/(RsPlusRh);
}

/** PCSデータによるパラメータ推定用  ************************************************
 *  関数呼び出し：setVars, calcA0B0
 *  関数setVarsを通して、現在の
 *  入力：Rs0init, Rs1init, eta0init, eta1init, VocSTC, IscSTC, VmpSTC, ImpSTC, Ns,
 *        etad, alphaSC, betaOc, gammaMp
 *  に対する、各G=GG[0...K-1], T=TT[0,...K-1]に対する5パラメータを計算する
 *  出力：G, T, aNsVt, aIsc, aVoc, aImp, aVmp, aRs, aEta, aRh, aI0, aIph  ***********/
double AdaptiveModel::calcPCS(double GG[], double TT[]){
    double Eret = 0.0;
    
    for(int k=0; k<GlobalConstants::K; k++){
	m_G = GG[k];
	m_T = TT[k];
	
	setVars();

	m_aG[k] = m_G;
	m_aT[k] = m_T;
	m_aNsVt[k] = m_NsVt;
	m_aIsc[k] = m_Isc;
	m_aVoc[k] = m_Voc;
	m_aImp[k] = m_Imp;
	m_aVmp[k] = m_Vmp;
	m_aLnGpGSTC[k] = m_lnGpGSTC;
	m_aW01[k] = m_W01;
	m_aW02[k] = m_W02;
	
	calcA0B0();

	// 現在のG,Tに対する５パラメータがセットされる
	double Emin = calc();
	if(Emin == std::numeric_limits<double>::infinity()) return Emin;
	else if(Eret < Emin) Eret = Emin;

	m_aRs[k] = m_Rs;
	m_aEta[k] = m_eta;
	m_aRh[k] = m_Rh;
	m_aI0[k] = m_I0;
	m_aIph[k] = m_Iph;
    }

    return Eret;
}

/** PCSデータによるパラメータ推定用,  StringGroup::getE用     ************************
 *  関数呼び出し：setVars, calcA0B0
 *  関数setVarsを通して、現在の
 *  入力：Rs0init, Rs1init, eta0init, eta1init, VocSTC, IscSTC, VmpSTC, ImpSTC, Ns,
 *        etad, alphaSC, betaOc, gammaMp
 *  に対する、各G=GG[0...K-1], T=TT[0,...K-1]に対する5パラメータを計算する
 *  出力：G, T, aNsVt, aIsc, aVoc, aImp, aVmp, aRs, aEta, aRh, aI0, aIph  ***********/
double AdaptiveModel::calcPCSGetE(){
    double Eret = 0.0;
    
    for(int k=0; k<GlobalConstants::K; k++){
	m_G = m_aG[k];
	m_T = m_aT[k];
	
	// setVarsGetE();
	m_Rs0 = m_Rs0init;
	m_Rs1 = m_Rs1init;
	m_eta0 = m_eta0init;
	m_eta1 = m_eta1init;
	
	m_GpGSTC = m_G/1000.0;
	m_lnGpGSTC = m_aLnGpGSTC[k];
	m_TmTSTC = m_T - 25.0;
	m_PmaxSTC = m_ImpSTC*m_VmpSTC;   // いらない？
	// Tkel = T+273.15;
	// double kT = 1.380649E-23*Tkel;
	m_NsVt = m_aNsVt[k];

	// double Eg = 1.8516E-19 -1.12488E-22*Tkel*Tkel/(Tkel+1108.0);
	// double Egd = -1.12488E-22*(Tkel*Tkel+2.0*Tkel*1108.0)/( (Tkel+1108.0)*(Tkel+1108.0) );
	// W01 = 3.0/Tkel -Egd/kT +Eg/(kT*Tkel);
	// W02 = (Eg/kT -3.0*logf(Tkel))*etad;
	m_W01 = m_aW01[k];
	m_W02 = m_aW02[k];

	m_Isc = m_IscSTC*m_GpGSTC*(1.0+m_alphaSc*m_TmTSTC);
	m_Voc = m_VocSTC/(1.0-m_deltaOc*m_lnGpGSTC)*(1.0+m_betaOc*m_TmTSTC);
	m_Vmp = m_VmpSTC*( 1.0/(1.0-m_deltaOc*m_lnGpGSTC) +m_betaOc*m_TmTSTC );
	double Pmax = m_PmaxSTC*m_GpGSTC/(1.0-m_deltaMp*m_lnGpGSTC)*(1.0+m_gammaMp*m_TmTSTC);
	m_Imp = Pmax/m_Vmp;

	calcA0B0();

	// 現在のG,Tに対する５パラメータがセットされる
	double Emin = calc();
	if(Emin == std::numeric_limits<double>::infinity()) return Emin;
	else if(Eret < Emin) Eret = Emin;

	m_aRs[k] = m_Rs;
	m_aEta[k] = m_eta;
	m_aRh[k] = m_Rh;
	m_aI0[k] = m_I0;
	m_aIph[k] = m_Iph;
    }

    return Eret;
}

/** 時点kに対する、getV, getIの準備をする    *****************************************
 *  入力： aRs[k], aEta[k], aRh[k], aI0[k], aIph[k], aNsVt[k]
 *  出力： Rs, eta, Rh, I0, Iph, NsVt       *****************************************/
void AdaptiveModel::getVIjunbi(int k){
    m_Rs = m_aRs[k];
    m_eta = m_aEta[k];
    m_Rh = m_aRh[k];
    m_I0 = m_aI0[k];
    m_Iph = m_aIph[k];
    m_NsVt = m_aNsVt[k];
}

/** 現在のVoc, Isc, Vmp, Impに対するRs, eta, Rh, I0, Iphをメッシュサーチして設定する 
 *  入力：Rs0, Rs1, NRs, eta0, eta1, NEta
 *  出力：Rs0, Rs1, eta0, eta1, I0, Iph
 *  関数呼び出し：SearchNoGPU, calcI0, calcIph  *************************************/
double AdaptiveModel::calcNoGPU(){
    // printf("calcNoGPU\n");///////////////t
    double Emin;
  
    for(int count=0; count<3; count++){  // 探索範囲変更のループ
    
	double kRs = std::exp( std::log(m_Rs1/m_Rs0) /(m_NRs-1) );
	double dEta = (m_eta1-m_eta0)/(m_NEta-1);

	// GPUを使わないサーチの実行
	Emin = SearchNoGPU(kRs, dEta);
	
	// 2024.3.20
	// ひとつもRh>=0のような条件を満たすRs, etaがなかったら、探索を終わる
	if(Emin == std::numeric_limits<double>::infinity()){
	    std::cout << "AdaptiveModel::calcNoGPU()$$ no suitable Rs and eta, Emin==INFINITY\n";///////////////////t
	    break;
	}
	
	// printf("count=%d, Emin=%g, Rs=%g, eta=%g\n", count, Emin, m_Rs, m_eta);/////////////////////t

	// 最適値が範囲の端から2番目以内だったら、探索範囲を取り直して再探索。
	if(this->m_Rs <= m_Rs0*pow(kRs, 1.5)){
	    double a = 1.05*m_Rs0/m_Rs1;
	    m_Rs0 *= a;
	    m_Rs1 *= a;
	    // printf("New smaller  Rs0=%g, Rs1=%g\n", m_Rs0, m_Rs1);////////////////t
	}else if(m_Rs1/pow(kRs, 1.5) <= this->m_Rs){
	    double b = 0.95*m_Rs1/m_Rs0;
	    m_Rs0 *= b;
	    m_Rs1 *= b;
	    // printf("New larger Rs0=%g, Rs1=%g\n", m_Rs0, m_Rs1);////////////////t
	}else if(m_eta1-1.5*dEta <= this->m_eta){
	    double eee = m_eta0;
	    m_eta0 = m_eta1-dEta;
	    m_eta1 = 2.0*m_eta1 -eee -dEta;
	    // printf("New larger  eta0=%g, eta1=%g\n", m_eta0, m_eta1);///////////////t
	}else{  // 端っこでなかったら、終わる
	    break;
	}
    }
  
    // パラメータRh, I0, Iphの値設定
    // Rh = calcRh(this->Rs, this->eta);
    this->m_I0 = calcI0(this->m_Rs, this->m_eta, this->m_Rh);
    this->m_Iph = calcIph(this->m_Rs, this->m_eta, this->m_Rh);

    printf("Emin=%g, Rs=%g, eta=%g, Rh=%g, I0=%g, Iph=%g\n", Emin, m_Rs, m_eta, m_Rh, m_I0, m_Iph);/////////////////////t
    
    return Emin;
}

/** [Rs0, Rs1]x[eta0, eta1]の範囲内でEを最小にするRs, etaを探索、設定してEminを返す **
 *  入力： Tkel, Voc, Vmp, Isc, Imp, Rs0, eta0, NRs, NEta, NsVt
 *  出力： Rs, eta, Rh     **********************************************************/
double AdaptiveModel::SearchNoGPU(double kRs, double dEta){
    double TInv = 1.0/m_Tkel;
    
    double VocMiVmp = m_Voc-m_Vmp;
    double IscMiImp = m_Isc-m_Imp;
  
    double Emin = std::numeric_limits<double>::infinity();

    double xRs = m_Rs0;  // Global Memroy内の配列からRs探索範囲の下限を読み取る
    for(int i=0; i<m_NRs; i++){ // Rsのループ
	
	double RsIsc = xRs*m_Isc;
	double RsImp = xRs*m_Imp;
	double VmpPluRsImp = m_Vmp +RsImp;
	double RsImpMiVmp = RsImp-m_Vmp;
    
	double xeta = m_eta0;  // Global Memroy内の配列からRs探索範囲の下限を読み取る
	for(int j=0; j<m_NEta; j++){ // etaのループ
	    
	    double etaNsVt = xeta*m_NsVt;
	    double x =m_Voc/etaNsVt, y=RsIsc/etaNsVt, z=VmpPluRsImp/etaNsVt;

	    if(x==y) continue;
	    
	    // Rhの計算
	    double Rh;
	    double ZpXY = NAN;  // Z/(X-Y)もついでにここで計算しておく
	    double eyx, ezx, exy, ezy, exz, eyz;
	    
	    if(x>=y && x>=z){
		eyx = std::exp(y-x);
		ezx = std::exp(z-x);
		Rh = ( ( m_Vmp +VocMiVmp*eyx -m_Voc*ezx )
		       / ( IscMiImp +m_Imp*eyx -m_Isc*ezx ) ) - xRs;
		ZpXY = ezx/(1.0-eyx);
	    }else if(y >= z){
		exy = std::exp(x-y);
		ezy = std::exp(z-y);
		Rh = ( ( m_Vmp*exy +VocMiVmp -m_Voc*ezy )
		       / ( IscMiImp*exy +m_Imp -m_Isc*ezy ) ) - xRs;
		ZpXY = ezy/(exy-1.0);
	    }else{
		exz = std::exp(x-z);
		eyz = std::exp(y-z);
		Rh = ( ( m_Vmp*exz +VocMiVmp*eyz -m_Voc )
		       / ( IscMiImp*exz +m_Imp*eyz -m_Isc ) ) - xRs;
	    }

	    if(std::isnan(Rh) || Rh<=0.0) continue;

	    // I0の計算
	    // double I01 = (1.0 +xRs/Rh )*m_Isc -m_Voc/Rh;
	    double I01 = ( (Rh+xRs)*m_Isc -m_Voc )/Rh;
      
	    double oneXY, XoneXY, XpXY;
	    if(x>=y){
		double eyx1 = std::exp(y-x);
		double emx = std::exp(-x);
		oneXY = emx/(1.0-eyx1);
		XoneXY = (1.0-emx)/(1.0-eyx1);
		XpXY = 1.0/(1.0-eyx1);
	    }else{
		double exy1 = std::exp(x-y);
		double emy = std::exp(-y);
		oneXY = emy/(exy1-1.0);
		XoneXY = (exy1-emy)/(exy1-1.0);
		XpXY = exy1/(exy1-1.0);
	    }
      
	    double I0 = I01*oneXY;
	    if(!std::isfinite(I0) || I0<0.0) continue;

	    // Iphの計算
	    double Iph = I01*XoneXY +m_Voc/Rh;
      
	    if(!std::isfinite(Iph) || Iph<0.0) continue;

	    // absAの計算
	    double P0 = RsImpMiVmp/etaNsVt*I01;
	    double P1 = (RsImpMiVmp +Rh*m_Imp )/Rh;

	    double absA, lnZpXY, XYpZ;
	    if( !std::isnan(ZpXY) ){  // zが最大、ではなかったら
		absA = std::abs(P0*ZpXY+P1);  // ZpXYは上のRhの計算のところで計算している、ｚが最大でない場合
	    }else{  // ｚが最大の場合
		if(x > y) lnZpXY = -std::log(exz-eyz);
		else lnZpXY = -std::log(eyz-exz);

		XYpZ = exz-eyz;

		double Atilder = lnZpXY +std::log( std::abs(P0 +P1*XYpZ) +std::abs(XYpZ) );
		absA = std::exp(Atilder)-1.0;
	    }
	    
	    if(!std::isfinite(absA)) continue;
      
	    // I0dの計算
	    double W0 = (m_W01 +m_W02/xeta)/xeta;
	    double I0d = W0*I0;
	    double etadPerEta = m_etad/xeta;

	    // Rsdの計算
	    double Rsd = (etadPerEta -m_alphaSc +TInv)*xRs +etaNsVt/m_Isc*(std::exp(-y)-1.0)*W0;

	    // Rhdの計算
	    double W1 = Rh*( m_Voc/etaNsVt*( etadPerEta -m_betaOc +TInv ) -W0 );
	    double W2 = (Rsd +(xRs+Rh)*m_alphaSc)*m_Isc +Rh*I0d -m_betaOc*m_Voc;  // +Rh*I0d modified 2024.1.21

	    double RhdPerRh = ( I01*XpXY*W1 +W2 )/(RsIsc-m_Voc);

	    // Iphdの計算
	    // double Iphd = Isc/Rh*Rsd -RsIsc/(Rh*Rh)*Rhd +(1.0+xRs/Rh)*alphaSc*Isc;
	    double Iphd = (Rsd -xRs*RhdPerRh +(xRs+Rh)*m_alphaSc)*m_Isc/Rh;

	    // absBの計算
	    double xxx = -m_Imp*Rsd -m_gammaMp*m_Vmp;
	    double W3 = I01*( ( xxx +VmpPluRsImp*(etadPerEta +TInv) )/etaNsVt -W0 );
	    double W4 = ( xxx +VmpPluRsImp*RhdPerRh )/Rh +Iphd +I0d;

	    double absB;
	    if( !std::isnan(ZpXY) ){  // zが最大、ではなかったら
		absB = std::abs(ZpXY*W3 +W4);   // ZpXYはｚが最大でないときは、Rhの計算のところで求めている
	    }else{  // ｚが最大のとき、このときlnZpXY, XYpZはAtilderのところで計算している
		double Btilder = lnZpXY +std::log( std::abs(W3+W4*XYpZ)+std::abs(XYpZ) );
		absB = std::exp(Btilder)-1.0;
	    }

	    if(!std::isfinite(absB)) continue;
	    
	    // Eの計算
	    double E = absA/m_absA0 +absB/m_absB0;
	    if(E < Emin){
		this->m_Rs = xRs;
		this->m_eta = xeta;
		this->m_Rh = Rh;
		Emin = E;
	    }
	    xeta += dEta;
	}
	xRs *= kRs;
    }
    
    return Emin;
}

/******************    *************************************************************/
void AdaptiveModel::set5params(double Rss, double etaa, double Rhh, double I00, double Iphh){
    m_Rs = Rss;
    m_eta = etaa;
    m_Rh = Rhh;
    m_I0 = I00;
    m_Iph = Iphh;
}

/********************     **********************************************************/
void AdaptiveModel::set5params(double Rss, double etaa, double Rhh, double I00, double Iphh, double NsVtt, int k){
    m_aRs[k] = Rss;
    m_aEta[k] = etaa;
    m_aRh[k] = Rhh;
    m_aI0[k] = I00;
    m_aIph[k] = Iphh;
    m_aNsVt[k] = NsVtt;
}

/******************************************************************************/
void AdaptiveModel::setGT(double G, double T){
    m_G = G;
    m_T = T;
    
    setVars();
    calcA0B0();
    calc();    // added Mar.4, 2026
    m_Imax = getI(m_Vmin);
    m_Vmax = getV(m_Imin);
}

/****** Suppose refreshment after set5params *******/
void AdaptiveModel::refreshImaxVmax(){
    m_Imax = getI(m_Vmin);
    m_Vmax = getV(m_Imin);
}    

/** STCパラメータをセットする   ****************************************************
 *  出力： IscSTC, VocSTC, ImpSTC, VmpSTC      *************************************/
void AdaptiveModel::setSTCparams(double X[4]){
    m_IscSTC = static_cast<double>(X[0]);
    m_VocSTC = static_cast<double>(X[1]);
    m_ImpSTC = static_cast<double>(X[2]*m_IscSTC);
    m_VmpSTC = static_cast<double>(X[3]*m_VocSTC);

    setVars();
    calcA0B0();
    calc();
    m_Imax = getI(m_Vmin);
    m_Vmax = getV(m_Imin);

    // printf("AdaptiveModel::setSTCparams, IscSTC=%g, VocSTC=%g, ImpSTC=%g, VmpSTC=%g\n", IscSTC, VocSTC, ImpSTC, VmpSTC);////////t
}

/**    ****************************************************************************/
void AdaptiveModel::print() const{
    printf("AdaptiveModel: IscSTC=%g, VocSTC=%g, ImpSTC=%g, VmpSTC=%g, alphaSc=%g, betaOc=%g, gammaMp=%g, G=%g, T=%g, Ns=%d\n", m_IscSTC, m_VocSTC, m_ImpSTC, m_VmpSTC, m_alphaSc, m_betaOc, m_gammaMp, m_G, m_T, m_Ns);
    for(int k=0; k<GlobalConstants::K; k++){
	printf("AdaptiveModel: k=%d, aG=%g, aT=%g, aNsVt=%g, aIsc=%g, aVoc=%g, aImp=%g, aVmp=%g, aRs=%g, aEta=%g, aRh=%g, aI0=%g, aIph=%g, aLnGpGSTC=%g, aW01=%g, aW02=%g\n",
	       k, m_aG[k], m_aT[k], m_aNsVt[k], m_aIsc[k], m_aVoc[k], m_aImp[k], m_aVmp[k], m_aRs[k], m_aEta[k], m_aRh[k], m_aI0[k], m_aIph[k], m_aLnGpGSTC[k], m_aW01[k], m_aW02[k]);
    }
}

/** Returen Five Parameter Rs, eta, Rh, I0, Iph    **********************************/
FiveParams AdaptiveModel::get5params() const{
    FiveParams fiveParams;

    for(int k=0; k<GlobalConstants::K; k++){
	fiveParams.aRs[k] = m_aRs[k];
	fiveParams.aEta[k] = m_aEta[k];
	fiveParams.aRh[k] = m_aRh[k];
	fiveParams.aI0[k] = m_aI0[k];
	fiveParams.aIph[k] = m_aIph[k];
    }
    
    return fiveParams;
}
    
/**    *******************************************************************************/
void FiveParams::print() const{
    for(int k=0; k<GlobalConstants::K; k++){
	printf("k=%d: Rs=%g, eta=%g, Rh=%g, I0=%g, Iph=%g\n", k, aRs[k], aEta[k], aRh[k],
	       aI0[k], aIph[k]);
    }
}
