// -*- C++ -*-
/************************************************************************************
 *              PCSData.cu
 ************************************************************************************/
#include <iostream>
#include <fstream>
#include <string>
#include "PCSData.hh"

/*************** Constructor *************************/
PCSData::PCSData(GlobalConstants &gc): m_gc(gc){
}

void PCSData::readFile(){
    readFile("Dat/GTIVTime.dat");
}

void PCSData::readFile(std::string fileName){
    std::fstream ifs {fileName, std::ios_base::in};
    if(!ifs){
	std::cout << "couldn't open '" << fileName << "' for reading" << std::endl;
	exit(-1);
    }

    std::string str;
    while(std::getline(ifs, str)){
	// std::cout << str << std::endl;

	std::vector<std::string> vStr = splitString(str, ' ');

	for(int i=0; i<vStr.size(); i++){
	    // std::cout << vStr[i] << " ";
	    double Ii = stod(vStr[5]);
	    double Vi = stod(vStr[7]);
	    if(Ii < 10.0 || Vi < 10.0) continue;  // I==0, V==0のことがあるので、読み飛ばす

	    if(i==1) G.push_back(stod(vStr[i])*1000.0);   // W/m^2単位にする
	    else if(i==3) T.push_back(stod(vStr[i])); 
	    else if(i==5) I.push_back(Ii); 
	    else if(i==7) V.push_back(Vi); 
	    else if(i==9) Time.push_back(vStr[i]); 
	}
	// std::cout << std::endl;
    }

    // for(int i=0; i<G.size(); i++){
    // 	std::cout << i << " " << G[i] << " " << T[i] << " " << I[i] << " " << V[i] << " " << Time[i] << std::endl;
    // }
}

/** G,T,I,Vの３０点平均をとる   *****************************************************
 *  ただし、日はまたがないようにする ***********************************************/
void PCSData::averageGTIV(){
    std::cout << "PCSData::averageGTIV$$ Take 30 average with PCSData::averageGTIV()\n";
    
    std::vector<double> Gave;
    std::vector<double> Tave;
    std::vector<double> Iave;
    std::vector<double> Vave;
    std::vector<std::string> TimeAve;

    int count = 0;
    std::string day, dayNow, TimeNow;
    double Gsum, Tsum, Isum, Vsum;
    
    for(int i=0; i<G.size(); i++){
	if(count==0){
	    Gsum = 0.0;
	    Tsum = 0.0;
	    Isum = 0.0;
	    Vsum = 0.0;

	    dayNow = (splitString(Time[i], '/'))[2];
	    // std::cout << "dayNow=" << dayNow << std::endl;
	    TimeNow = Time[i];
	}

	day = (splitString(Time[i], '/'))[2];
	if(day != dayNow){
	    if(count >= 20){
		Gave.push_back(Gsum/count);
		Tave.push_back(Tsum/count);
		Iave.push_back(Isum/count);
		Vave.push_back(Vsum/count);
		TimeAve.push_back(TimeNow);
	    }

	    count = 0;
	    Gsum = 0.0;
	    Tsum = 0.0;
	    Isum = 0.0;
	    Vsum = 0.0;

	    dayNow = (splitString(Time[i], '/'))[2];
	    TimeNow = Time[i];
	}	    
	
	Gsum += G[i];
	Tsum += T[i];
	Isum += I[i];
	Vsum += V[i];

	count++;

	// printf("i=%d count=%d Time=%s TimeNow=%s day=%s G=%g T=%g I=%g V=%g\n", i, count, Time[i], TimeNow, day, G[i], T[i], I[i], V[i]);////////////////////////t 

	if(count == 30){
	    Gave.push_back(Gsum/30.0);
	    Tave.push_back(Tsum/30.0);
	    Iave.push_back(Isum/30.0);
	    Vave.push_back(Vsum/30.0);
	    TimeAve.push_back(TimeNow);

	    // printf("G=%g T=%g I=%g V=%g Time=%s\n", Gave[Gave.size()-1], Tave[Tave.size()-1], Iave[Iave.size()-1], Vave[Vave.size()-1], TimeNow);////////////////////t 

	    count = 0;
	}
    }

    G.clear();
    T.clear();
    I.clear();
    V.clear();
    Time.clear();

    G = Gave;
    T = Tave;
    I = Iave;
    V = Vave;
    Time = TimeAve;
}

/******************************* Apr. 9, 2026-- ******************************/
void PCSData::makeSimPCSDataPmax(SimPCSData &sd){
    std::fstream fs {"Dat/GTIVTimeSim.dat", std::ios_base::out};
    fs << std::scientific << std::setprecision(20);
    
    StringGroup sg(m_gc, 1000.0, 25.0);
    AdaptiveModel* p_am = nullptr;

    std::string timeEvent {"2019/5/7/12:0"};
    // std::string timeEvent {"2019/5/7/9:5"};
    bool beforeEvent {true};
    for(int n=0; n<G.size(); n++){
	printf("%s\n", Time[n].c_str());///////////////////////////
	sg.setGT(G[n], T[n]);
	
	if(beforeEvent && Time[n]==timeEvent){
	    beforeEvent = false;
	    addFault(sg, G[n], T[n]);

	    p_am = &(sg.m_vsm[0].m_vss[0].m_vam[0]);
	}
	
	if(beforeEvent == false){
	    p_am->set5params(10.0, p_am->getEta(), p_am->getRh(), p_am->getI0(),
			     p_am->getIph());
	    // p_am->set5params(0.1, p_am->getEta(), p_am->getRh(), p_am->getI0(),
	    // 		     p_am->getIph());
	    sg.refreshMaxMin();
	    // printf("am.Rs=%g\n", p_am->getRs());///////////////////
	}
	
	double VmpSim, ImpSim;
	sg.getPmax(VmpSim, ImpSim);

	fs << G[n] << ", " << T[n] << ", " << I[n] << ", " << V[n]
	   << ", " << Time[n] << ", " << ImpSim << ", " << VmpSim << std::endl;

	// if(n==10) break; //////////////////////////////////
    }
    fs.close();
    printf("file 'Dat/GTIVTimeSim.dat' saved.\n");////////////////////////
}

/****************************************************************/
void PCSData::addFault(StringGroup &sg, const double G, const double T){
    AdaptiveModel am(m_gc, G, T);
    Substring ss(m_gc, G, T);
    ss.addAM(am, 1);
    StringModel sm(m_gc, G, T);
    sm.addSS(ss, 1);
    sg.addSM(sm, 1);

    // // Broken or deteriorated cells, substrings, strings
    // AdaptiveModel am(m_gc, G, T);
    // Substring ss(m_gc, G, T);
    // ss.addAM(am, 10);
    // StringModel sm(m_gc, G, T);
    // sm.addSS(ss, 3*7);
    
    // sg.addSM(sm, 20);

    // AdaptiveModel &amx {sg.m_vsm[0].m_vss[0].m_vam[0]};
    // amx.set5params(0.1, amx.getEta(), amx.getRh(), amx.getI0(),
    // 		   amx.getIph());
}

/********    文字列strを文字cで分割する 連続したcはひとつとみなす    ***************/
std::vector<std::string> PCSData::splitString(std::string str, char c){
    std::vector<std::string> ret;

    unsigned int i0 = 0;
    bool final = false;
    
    while(true){
	int i1 = str.find_first_of(c, i0);

	if(i0==i1){   // 最初の文字がスペース
	    i0 = i1+1;
	    continue;
	}else if(i1 == std::string::npos){
	    i1 = str.size();
	    final = true;
	}
	
	std::string subStr {str, i0, i1-i0};
	ret.push_back(subStr);
	
	// std::cout << "i0= " << i0 << " i1= " << i1 << " " << subStr << std::endl;

	i0 = i1+1;

	if(final) break;
    }

    return ret;
}
