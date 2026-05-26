// -*- C++ -*-
/************************************************************************************
 *    main.cu
 ************************************************************************************/
#include <iostream>
#include <iomanip>  // setprecision
#include <chrono>
#include <ctime>
#include <cuda_runtime.h>
#include <fstream>
#include <ios>

#include "PCSData.hh"
#include "adaptiveModel.hh"  /////////////////t
// #include "LambertWhara.hh"
#include "globalConstants.hh"
// #include "substring.hh"  ///////////////////t
// #include "stringModel.hh"
#include "stringGroup.hh"

/**  PCS測定データを用いたストリング群のパラメータ推定 */
void stringGroupParamEstim(GlobalConstants gc);

/*****************************   メイン関数   ***************************************/
int main(void){
    GlobalConstants gc;
    
    // GPU関係の定数
    cudaDeviceProp prop;
  
    cudaError_t cudaError = cudaGetDeviceProperties(&prop, 0);
    if(cudaError == cudaSuccess){
	std::cout << "This PC has a GPU !" << std::endl;
	gc.haveGPU = true;
    }else{
	std::cout << "This PC has no GPU !" << std::endl;
	gc.haveGPU = false;
    }

    if(gc.haveGPU){
	std::cout << "# Device: " << prop.name << std::endl;
	std::cout << "#   totalGlobalMem: " << prop.totalGlobalMem << std::endl;
	std::cout << "#   multiProcessorCount: " << prop.multiProcessorCount << std::endl;
	std::cout << "#   maxBlocksPerMultiProcessor: " << prop.maxBlocksPerMultiProcessor << std::endl;
	std::cout << "#   maxThreadsPerMultiProcessor: " << prop.maxThreadsPerMultiProcessor << std::endl;
	std::cout << "#   maxThreadsPerBlock: " << prop.maxThreadsPerBlock << std::endl;
    }
    
    // 実行時間計測の準備
    auto start = std::chrono::system_clock::now();
    std::cout << std::setprecision(20); //////////////////t

    stringGroupParamEstim(gc);
    
    // 実行時間計測
    auto end = std::chrono::system_clock::now();
    auto nanosec = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();
    std::cout << "# Execution time of AdaptiveModel.calc()= " << nanosec/1.0E9 << " sec" << std::endl;
}

/**************  PCS測定データを用いたストリング群のパラメータ推定 ******************/
void stringGroupParamEstim(GlobalConstants gc){
    PCSData pd;
    pd.readFile();
    pd.averageGTIV();

    double G[GlobalConstants::K], T[GlobalConstants::K];
    double I[GlobalConstants::K], V[GlobalConstants::K];
    std::string Time[GlobalConstants::K];

    StringGroup sg {gc, 1000.0, 25.0};

    for(int i=0; i<=pd.G.size()-gc.K; i+=gc.K){
	for(int k=0; k<gc.K; k++){
	    G[k] = pd.G[i+k];
	    T[k] = pd.T[i+k];
	    I[k] = pd.I[i+k];
	    V[k] = pd.V[i+k];
	    Time[k] = pd.Time[i+k];
	    
	    printf("k=%d, %s G=%g, T=%g, I=%g, V=%g\n", k, Time[k].c_str(), G[k], T[k], I[k], V[k]);//////////t
	}
	
	sg.paramEstim(G, T, I, V, Time);
	
	break;////////////////////t
    }	
}

    
// /** メイン関数 */
// int main(void){
//     GlobalConstants gc;
    
//      // GPU関係の定数
//     cudaDeviceProp prop;
  
//     cudaError_t cudaError = cudaGetDeviceProperties(&prop, 0);
//     if(cudaError == cudaSuccess){
// 	std::cout << "This PC has a GPU !" << std::endl;
// 	gc.haveGPU = true;
//     }else{
// 	std::cout << "This PC has no GPU !" << std::endl;
// 	gc.haveGPU = false;
//     }

//     if(gc.haveGPU){
// 	std::cout << "# Device: " << prop.name << std::endl;
// 	std::cout << "#   totalGlobalMem: " << prop.totalGlobalMem << std::endl;
// 	std::cout << "#   multiProcessorCount: " << prop.multiProcessorCount << std::endl;
// 	std::cout << "#   maxBlocksPerMultiProcessor: " << prop.maxBlocksPerMultiProcessor << std::endl;
// 	std::cout << "#   maxThreadsPerMultiProcessor: " << prop.maxThreadsPerMultiProcessor << std::endl;
// 	std::cout << "#   maxThreadsPerBlock: " << prop.maxThreadsPerBlock << std::endl;
//     }
    
//     // 実行時間計測の準備
//     auto start = std::chrono::system_clock::now();

//    std::cout << std::setprecision(20); //////////////////t

//     // AdaptiveModel am {gc, 500.0, 41.2};
//     // // am.set5params(1.21594e-07, 1.50072, 6.953e+07, 4.12134e-06, 4.26915);
//     // // printf("am.getV(0)=%g\n", am.getV(0.0));
//     // am.set5params(1.21594e-07, 1.50072, 6.953e+07, 4.12134e-06, 4.26915, 20.0*1.380649E-23*(41.2113+273.15)/1.602176634E-19, 2);
//     // printf("am.getV(0, 2)=%g\n", am.getV(0.0, 2));
//     // // exit(0);

//     stringGroupParamEstim(gc);
    
//     // const int max = 80;
//     // char str[max];
//     // time_t t = time(nullptr);
//     // tm* pt = localtime(&t);
//     // strftime(str, max, "%D, %H:%M\n", pt);
//     // printf(str);

//     // tm startTime {0, 0, 8, 30, 3, 118};
//     // strftime(str, max, "%D, %H:%M\n", &startTime);
//     // printf(str);

//     // startTime.tm_mday++;
//     // strftime(str, max, "%D, %H:%M\n", &startTime);
//     // printf(str);

//     // /** StringGroup */
//     // StringGroup sg {gc, 1000.0, 25.0};
    
//     // std::fstream ofs {"IV_getV.dat", std::ios_base::out};
//     // for(double I=-30.0; I<= 1000.0; I += 0.1){
//     // // for(double I=18.40; I<= 20.0; I += 0.01){
//     // 	// std::cout << std::setprecision(30) << ss.getV(I) << " " << I << std::endl;
//     // 	ofs << std::scientific << sg.getV(I) << " " << I << std::endl;
//     // }
//     // ofs.close();

//     // std::fstream ofs1 {"IV_getI.dat", std::ios_base::out};
//     // // for(double V=-0.608; V<= 13.0; V += 0.01){
//     // for(double V= -25.0; V<= 600.0; V += 0.1){
//     // // for(double V=-50; V<= 600.0; V += 0.1){
//     // 	// std::cout << std::scientific << V << " " << sm.getI(V) << std::endl;
//     // 	ofs1 << std::scientific << V << " " << sg.getI(V) << std::endl;
//     // }
//     // ofs1.close();

//     // double W = utl::LambertW0hara(2.0);
//     // std::cout << std::setprecision(20) << W << std::endl;
//     // W = utl::LambertW0haraAexp(2.3, 26.0);
//     // std::cout << W << std::endl;
//     // exit(0);

//     /** StringModel */
//     // StringModel sm {gc, 1000.0, 25.0};

//     // sm.addSS(700.0, 25.0, gc.IscSTC, gc.VocSTC/gc.Ns*gc.Ns_sub, gc.ImpSTC, gc.VmpSTC/gc.Ns*gc.Ns_sub);
//     // sm.addSS(700.0, 25.0, gc.IscSTC, gc.VocSTC/gc.Ns*gc.Ns_sub, gc.ImpSTC, gc.VmpSTC/gc.Ns*gc.Ns_sub);
//     // sm.addSS(700.0, 25.0, gc.IscSTC, gc.VocSTC/gc.Ns*gc.Ns_sub, gc.ImpSTC, gc.VmpSTC/gc.Ns*gc.Ns_sub);
//     // sm.addSS(700.0, 25.0, gc.IscSTC, gc.VocSTC/gc.Ns*gc.Ns_sub, gc.ImpSTC, gc.VmpSTC/gc.Ns*gc.Ns_sub);
//     // sm.addSS(700.0, 25.0, gc.IscSTC, gc.VocSTC/gc.Ns*gc.Ns_sub, gc.ImpSTC, gc.VmpSTC/gc.Ns*gc.Ns_sub);

//     // std::fstream ofs {"IV_getV.dat", std::ios_base::out};
//     // for(double I=-10.0; I<= 20.0; I += 0.01){
//     // // for(double I=18.40; I<= 20.0; I += 0.01){
//     // 	// std::cout << std::setprecision(30) << ss.getV(I) << " " << I << std::endl;
//     // 	ofs << std::scientific << sm.getV(I) << " " << I << std::endl;
//     // }
//     // ofs.close();
    
//     // // double I = sm.getI(V);
//     // // std::cout << std::setprecision(20) << "I= " << I << std::endl;

//     // // std::cout << sm.getI(375.2846) << std::endl;/////////////////////
    
//     // std::fstream ofs1 {"IV_getI.dat", std::ios_base::out};
//     // // for(double V=-0.608; V<= 13.0; V += 0.01){
//     // for(double V= sm.getV(30.0); V<= 600.0; V += 0.1){
//     // // for(double V=-50; V<= 600.0; V += 0.1){
//     // 	// std::cout << std::scientific << V << " " << sm.getI(V) << std::endl;
//     // 	ofs1 << std::scientific << V << " " << sm.getI(V) << std::endl;
//     // }
//     // ofs1.close();

//     /** Substring */
//     // // Substring ss {gc, 1000.0, 25.0};
//     // Substring ss {gc, 1000.0f, 25.0f, gc.IscSTC*0.8f, gc.VocSTC, gc.ImpSTC*0.8f, gc.VmpSTC};
    
//     // //std::cout << std::scientific << "V=" << ss.getV(18.4199999999999128874605958117) << std::endl;
//     // // std::cout << std::scientific << "V=" << ss.getV(18.43) << std::endl;

//     // // double x2= 8.6390624999999587885;
//     // // double xmid= 8.6030664062499582201;
//     // // // std::cout << std::setprecision(30) << (x2+xmid)*0.5 << " " << ss.func((x2+xmid)*0.5) << std::endl;
//     // // std::cout << std::setprecision(30) << "x= " << (x2+xmid)*0.5 << std::endl;
//     // // std::cout << ss.am.getV((x2+xmid)*0.5) << std::endl;
//     // // exit(0);
    
//     // std::fstream ofs {"IV_getV.dat", std::ios_base::out};
//     // for(double I=-10.0; I<= 20.0; I += 0.01){
//     // // for(double I=18.40; I<= 20.0; I += 0.01){
//     // 	// std::cout << std::setprecision(30) << ss.getV(I) << " " << I << std::endl;
//     // 	ofs << std::scientific << ss.getV(I) << " " << I << std::endl;
//     // 	// ofs << std::scientific << sm.getV(I) << " " << I << std::endl;
//     // }
//     // ofs.close();
    
//     // std::fstream ofs1 {"IV_getI.dat", std::ios_base::out};
//     // // for(double V=-0.608; V<= 13.0; V += 0.01){
//     // for(double V=-0.608; V<= 40.0; V += 0.01){
//     // // for(double V=-50; V<= 600.0; V += 0.1){
//     // 	// std::cout << std::scientific << ss.getV(I) << " " << I << std::endl;
//     // 	ofs1 << std::scientific << V << " " << ss.getI(V) << std::endl;
//     // 	// ofs1 << std::scientific << V << " " << sm.getI(V) << std::endl;
//     // }
//     // ofs1.close();
	
//     //    std::cout << "ss.getV(-1.0)= " << ss.getV(-1.0) << std::endl;
  
//     /** Adaptive Model  */
//     // // AdaptiveModel am {8.61, 37.05, 8.16, 29.42, 0.00038,
//     // // 		    -0.00329, -0.00440, 1000.0, 25.0, 60};
//     // AdaptiveModel am {gc.IscSTC, gc.VocSTC/gc.Ns*gc.Ns_sub, gc.ImpSTC, gc.VmpSTC/gc.Ns*gc.Ns_sub,
//     // 		      gc.alphaSc, gc.betaOc, gc.gammaMp, 1000.0, 25.0, gc.Ns_sub};

//     // // AdaptiveModel am {gc, 700.0, 25.0};
    
//     // float Emin = am.calc();
    
//     // printf("Rs=%g, eta=%g, Rh=%g, I0=%g, Iph=%g, Emin=%g\n",
//     // 	   am.getRs(),am.getEta(), am.getRh(), am.getI0(), am.getIph(), Emin);
    
//     // // // std::fstream ofs {"IV_getV.dat", std::ios_base::out};
//     // // // for(double I=-10.0; I<= am.getIph(); I += 0.01){
//     // // // // for(double I=18.40; I<= 20.0; I += 0.01){
//     // // // 	// std::cout << std::setprecision(30) << ss.getV(I) << " " << I << std::endl;
//     // // // 	ofs << std::scientific << am.getV(I) << " " << I << std::endl;
//     // // // 	// ofs << std::scientific << sm.getV(I) << " " << I << std::endl;
//     // // // }
//     // // // ofs.close();

//     // std::fstream ofs1 {"IV_getI.dat", std::ios_base::out};
//     // // for(double V=-0.608; V<= 13.0; V += 0.01){
//     // for(double V=-0.608; V<= am.getV(0.0)+0.1; V += 0.01){
//     // // for(double V=-50; V<= 600.0; V += 0.1){
//     // 	// std::cout << std::scientific << ss.getV(I) << " " << I << std::endl;
//     // 	ofs1 << std::scientific << V << " " << am.getI(V) << std::endl;
//     // 	// ofs1 << std::scientific << V << " " << sm.getI(V) << std::endl;
//     // }
//     // ofs1.close();

//     // 実行時間計測
//     auto end = std::chrono::system_clock::now();
//     auto nanosec = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();
//     std::cout << "# Execution time of AdaptiveModel.calc()= " << nanosec/1.0E9 << " sec" << std::endl;
  
//     return 0;
// }

// /**  PCS測定データを用いたストリング群のパラメータ推定 */
// void stringGroupParamEstim(GlobalConstants gc){
// //     k=2 G=770.4000244 T=57.41066742 Rs=7.880382327e-06 eta=1.671120167 Rh=-127617120 I0=7.309637294e-05 Iph=6.714838505
// // Voc=10.88013649 Isc=6.714838505 Vmp=8.625005722 Imp=6.086348534
    
//     // // AdaptiveModel am(gc, 770.4000244, 57.41066742);
//     // AdaptiveModel am(gc, 800.0, 35.0);
//     // am.calc();
//     // for(double V=-5.0; V<=10.0; V+=0.1){
//     // 	double I = am.getI(V);
//     // 	double V1 = am.getV(I);
//     // 	printf("V=%g, I=%g, V1=%g\n", V, I, V1);
//     // }
//     // exit(0);//////////////////////t
    
//     // // AdaptiveModel am1(gc, 770.4000244, 57.41066742);
//     // AdaptiveModel am1(gc, 800.0, 35.0);
//     // am1.calcNoGPU();
//     // // am.calcRh(7.880382327e-06, 1.671120167);
//     // exit(0);

//     PCSData pd;
//     pd.readFile();

//     pd.averageGTIV();
    
//     // for(int i=0; i<pd.G.size(); i++){
//     // 	// for(int i=0; i<200; i++){
//     // 	std::cout << i << " " << pd.G[i] << " " << pd.T[i] << " " << pd.I[i] << " " << pd.V[i] << " " << pd.Time[i] << std::endl;
//     //    }
//     // exit(0);
    
//     // StringGroup sg {gc, 1000.0, 25.0};
//     // // StringGroup sg {gc, 324.8, 30.2547};//408.563
//     // // printf("sg I=%g, %g\n", sg.getI(408.563), sg.getI(408.563)/gc.Sg);//////////////
//     // exit(0);/////////////////

//     // gc.haveGPU = false;///////////////////////////////t
    
//     // StringModel sm {gc, 1000.0, 25.0};
//     StringGroup sg {gc, 1000.0, 25.0};
   
//     for(int i=0; i<=pd.G.size()-gc.K; i+=gc.K){
//     // for(int i=0; i<=30; i+=gc.K){////////////////t
//     // for(int i=980; i<=pd.G.size()-gc.K; i+=gc.K){//////////////////////t
//     // for(int i=1240; i<=pd.G.size()-gc.K; i+=gc.K){//////////////////////t
//     // for(int i=64; i<=pd.G.size()-gc.K; i+=pd.G.size()){
// 	double G[GlobalConstants::K], T[GlobalConstants::K];
// 	double I[GlobalConstants::K], V[GlobalConstants::K];
// 	std::string Time[GlobalConstants::K];

// 	for(int k=0; k<gc.K; k++){
// 	    G[k] = pd.G[i+k];
// 	    T[k] = pd.T[i+k];
// 	    I[k] = pd.I[i+k];
// 	    V[k] = pd.V[i+k];
// 	    Time[k] = pd.Time[i+k];
// 	    // printf("k=%d, G=%g, T=%g, I=%g, V=%g, ", k, G[k], T[k], I[k], V[k]);//////////t
// 	    // std::cout << Time[k] << std::endl;/////////////t
// 	}

// 	// AdaptiveModel am {gc, 1000.0, 25.0};/////////////////t
// 	// am.calcPCS(G, T);/////////////////t
// 	// for(int k=0; k<gc.K; k++){///////////////////t
// 	//     double Vk = V[k]/gc.L;
// 	//     double Ik = I[k]/gc.S;
// 	//     double Isim = am.getI(Vk, k);
// 	//     double Vsim = am.getV(Ik, k);
// 	//     printf("k=%d, G=%g, T=%g, Vk=%g, Vsim=%g, Ik=%g, Isim=%g\n", k, G[k], T[k], Vk, Vsim, Ik, Isim);
// 	// }
// 	// exit(0);/////////////////////t

// 	// Substring ss {gc, 1000.0, 25.0};
// 	// ss.setGTetc(G, T);

// 	// for(int k=0; k<gc.K; k++){
// 	//     // double I = ss.getI(V);
// 	//     ss.getVIjunbi(k);
// 	//     double Vsim = ss.getV(I[k]/gc.S);
// 	//     ss.getVIafter(k);
// 	//     double Isim = ss.getI(Vsim);
// 	//     // printf("V=%g, I=%g, bypassDiodeOn=%d, Vsim=%g\n", V[k], I[k], ss.getBypassDiodeOn(), Vsim);
// 	//     printf("k=%d, %s G=%g, T=%g, V=%g, I/S=%g, Vsim=%g, bypassDiodeOn=%d, Isim=%g\n", k, Time[k].c_str(), G[k], T[k], V[k], I[k]/gc.S, Vsim, ss.getBypassDiodeOn(k), Isim);
// 	// }
// 	// exit(0);

// 	// sm.setVars_k(G, T);
// 	// // for(double I = -1.0; I < 10.0; I+=0.1){
// 	// //     printf("I=%g, Vsim=%g\n", I, sm.getV(I));
// 	// // }
// 	// for(int k=0; k<GlobalConstants::K; k++){
// 	//     sm.getVIjunbi(k);
// 	//     double Vsim = sm.getV(I[k]/gc.S);
// 	//     double Isim = sm.getI(Vsim);
// 	//     printf("k=%d, %s G=%g, T=%g, I=%g, V=%g, Vsim=%g, Isim=%g\n", k, Time[k].c_str(), G[k], T[k], I[k]/gc.S, V[k], Vsim, Isim);
// 	// }
// 	// if(i>=0) exit(0);
	
// 	for(double V = -50.0; V < 800.0; V+=2.0){
// 	    double Isim = sg.getI(V);
// 	    double Vsim = sg.getV(Isim);
// 	    printf("V=%g, Isim=%g, Vsim=%g\n", V, Isim, Vsim);
// 	}
// 	exit(0);
    
// 	if(i%100==0){
// 	    // if(i>=0){
// 	    std::cout << "i=" << i << " " << Time[0] << std::endl;///////////////t
// 	}
	
// 	// sg.paramEstim(G, T, I, V, Time);
//     }
// }
