// -*- C++ -*-
/************************************************************************************
 *    main.cu
 ************************************************************************************/
#include <iostream>
#include <iomanip>  // setprecision
#include <chrono>
#include <ctime>
#if !defined(UBUNTU)
#include <cuda_runtime.h>
#endif
#include <fstream>
#include <ios>

#include "globalConstants.hh"
#include "PCSData.hh"
// #include "adaptiveModel.hh"  /////////////////t
// #include "LambertWhara.hh"
// #include "substring.hh"  ///////////////////t
// #include "stringModel.hh"
// #include "stringGroup.hh"

/**************************************************/
void stringGroupParamEstim(GlobalConstants gc);

/*****************************   ***************************************/
int main(void){
    bool haveGPU {false};
    
#if !defined(UBUNTU)
    // GPU information
    cudaDeviceProp prop;
  
    cudaError_t cudaError = cudaGetDeviceProperties(&prop, 0);
    if(cudaError == cudaSuccess){
	std::cout << "This PC has a GPU !" << std::endl;
	haveGPU = true;
    }else{
	std::cout << "This PC has no GPU !" << std::endl;
	haveGPU = false;
    }

    if(haveGPU){
	std::cout << "# Device: " << prop.name << std::endl;
	std::cout << "#   totalGlobalMem: " << prop.totalGlobalMem << std::endl;
	std::cout << "#   multiProcessorCount: " << prop.multiProcessorCount << std::endl;
	std::cout << "#   maxBlocksPerMultiProcessor: " << prop.maxBlocksPerMultiProcessor << std::endl;
	std::cout << "#   maxThreadsPerMultiProcessor: " << prop.maxThreadsPerMultiProcessor << std::endl;
	std::cout << "#   maxThreadsPerBlock: " << prop.maxThreadsPerBlock << std::endl;
    }
#endif
    
    GlobalConstants gc(haveGPU);
    
    auto start = std::chrono::system_clock::now();
    std::cout << std::setprecision(20); //////////////////t

    // double G {1000.0};
    // double Tdeg {25.0};

    PCSData pd(gc);
    SimPCSData spd;
    pd.readFile("Dat/GTIVTime201905.dat");
    pd.makeSimPCSData05(spd);

    // gc.haveGPU = false;///////////////////////////
    
    // AdaptiveModel am(gc, G, Tdeg);
    // am.set5params(0.1, am.getEta(), am.getRh()*0.01, am.getI0(),
    // 		  am.getIph()*0.5);

    // std::fstream fs {"Graphs/am_IV.dat", std::ios_base::out};
    // for(double V=-0.5; V <= 0.8; V += 0.01){
    // 	fs << std::scientific << std::setprecision(20);
    // 	fs << "V= " << V << " ,I= " << am.getI(V) << std::endl;
    // }
    // // for(double I=-1.0; I < am.getImax(); I += 0.01){
    // // 	fs << std::scientific << std::setprecision(20);
    // // 	fs << "V= " << am.getV(I) << " ,I= " << I << std::endl;
    // // }
    // fs.close();

    // Substring ss(gc, G, Tdeg);
    // // double X[4]={gc.IscSTC*0.8, gc.VocSTC/gc.Ns, 0.7, 0.8};
    // // double X[4]={gc.IscSTC*0.5, gc.VocSTC/gc.Ns, gc.ImpSTC/gc.IscSTC,
    // // 		 gc.VmpSTC/gc.VocSTC};
    // // ss.addAM(X, 20);
    // ss.addAM(am, 10);

    // std::fstream fs_ss {"Graphs/ss_IV.dat", std::ios_base::out};
    // // for(double I=-1.0; I <= 10.0; I += 0.01){
    // // 	fs_ss << std::scientific << std::setprecision(20);
    // // 	fs_ss << "V= " << ss.getV(I) << " ,I= " << I << std::endl;
    // // }
    // // fs_ss.close();
    // for(double V=-1.0; V <= 15.0; V += 0.01){
    // 	fs_ss << std::scientific << std::setprecision(20);
    // 	fs_ss << "V= " << V << " ,I= " << ss.getI(V) << std::endl;
    // }
    // fs_ss.close();
    
    // StringModel sm(gc, G, Tdeg);
    // sm.addSS(ss, 3*7);

    // std::fstream fs {"Graphs/sm_IV.dat", std::ios_base::out};
    // // for(double I=-1.0; I <= 12.0; I += 0.01){
    // // 	fs << std::scientific << std::setprecision(20);
    // // 	fs << "V= " << sm.getV(I) << " ,I= " << I << std::endl;
    // // }
    // for(double V=sm.getV(20.0); V <= sm.getV(-sm.getIb0()*0.9); V += 0.1){
    // 	fs << std::scientific << std::setprecision(20);
    // 	fs << "V= " << V << " ,I= " << sm.getI(V) << std::endl;
    // }
    // fs.close();
    
    // StringGroup sg(gc, G, Tdeg);
    // // sg.addSM(sm, 40);

    // double Vmp, Imp;
    // double Pmax = sg.getPmax(Vmp, Imp);
    // printf("StringGroup Vmp=%g, Imp=%g, Pmax=%g\n", Vmp, Imp, Pmax);

    // // std::fstream fs {"Graphs/sg_IV.dat", std::ios_base::out};
    // std::fstream fs {"Graphs/sg_PV.dat", std::ios_base::out};
    // // for(double I=sg.getImin(); I <= sg.getImax(); I += 1.0){
    // // 	fs << std::scientific << std::setprecision(20);
    // // 	fs << "V= " << sg.getV(I) << " ,I= " << I << std::endl;
    // // }
    // // for(double V=sg.getVmin(); V <= sg.getVmax()*1.2; V += 1.0){
    // // 	fs << std::scientific << std::setprecision(20);
    // // 	fs << "V= " << V << " ,I= " << sg.getI(V) << std::endl;
    // // }
    // for(double V=sg.getVmin(); V <= sg.getVmax()*1.2; V += 1.0){
    // 	fs << std::scientific << std::setprecision(20);
    // 	fs << "V= " << V << " ,P= " << sg(V) << std::endl;
    // }
    // fs.close();
    
    auto end = std::chrono::system_clock::now();
    auto nanosec = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();
    std::cout << "# Execution time of AdaptiveModel.calc()= " << nanosec/1.0E9 << " sec" << std::endl;
}

/**************  ******************/
void stringGroupParamEstim(GlobalConstants gc){
    PCSData pd(gc);
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
