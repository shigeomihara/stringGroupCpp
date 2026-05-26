// -*- C++ -*-
/*********************************************************************
 *      searchHistoryAM.hh
 ********************************************************************/
#include "searchHistoryAM.hh"
#include <iostream>
#include <fstream>
#include <string>
#include "PCSData.hh"

/***************** Constructor *****************************/
SearchHistoryAM::SearchHistoryAM(bool haveGPU){
    if(haveGPU) m_fileName = "Dat/searchHistoryAM_GPU.dat";
    else m_fileName = "Dat/searchHistoryAM_NoGPU.dat";
    readFile();
}

/***************** Destructor *****************************/
SearchHistoryAM::~SearchHistoryAM(){
    saveFile();
}

/**********************************************************/
int SearchHistoryAM::haveSearchResult(SearchResult &sr){
    for(int i=0; i<m_vsr.size(); i++){
	if(m_vsr[i].equal(sr)) return i;
    }
    return -1;
}

/**********************************************************/
bool SearchResult::equal(SearchResult &sr){
    // printf("Isc=%g, sr.Isc=%g\n", Isc, sr.Isc);//////////////////////
    if(Isc==sr.Isc && Voc==sr.Voc && Imp==sr.Imp
       && Vmp==sr.Vmp && Tkel==sr.Tkel
       && absA0==sr.absA0 && absB0==sr.absB0) return true;
    else return false;
}

/***************************************************************/
void SearchHistoryAM::readFile(){
    std::fstream fs {m_fileName, std::ios_base::in};
    if(!fs){
	std::cout << "couldn't open '"<< m_fileName << "' for reading" << std::endl;
	return;
    }

    std::string str;
    while(std::getline(fs, str)){
	std::vector<std::string> vStr = PCSData::splitString(str, ',');

	SearchResult sr;
	sr.Isc = stod(vStr[0]);
	sr.Voc = stod(vStr[1]);
	sr.Imp = stod(vStr[2]);
	sr.Vmp = stod(vStr[3]);
	sr.Tkel = stod(vStr[4]);
	sr.absA0 = stod(vStr[5]);
	sr.absB0 = stod(vStr[6]);
	sr.Rs = stod(vStr[7]);
	sr.eta = stod(vStr[8]);
	sr.Rh = stod(vStr[9]);
	sr.I0 = stod(vStr[10]);
	sr.Iph = stod(vStr[11]);
	sr.Emin = stod(vStr[12]);

	m_vsr.push_back(sr);
	m_fileLineNum++;
    }
    printf("SearchHistoryAM::readFile(), file '%s' read.\n",
	   m_fileName.c_str());///////////////////
}

/***************************************************************/
void SearchHistoryAM::saveFile(){
    std::fstream fs {m_fileName, std::ios_base::app};
    fs << std::scientific << std::setprecision(20);
	
    for(int i=m_fileLineNum; i<m_vsr.size(); i++){
	fs << m_vsr[i].Isc << ", " << m_vsr[i].Voc << ", " << m_vsr[i].Imp << ", "
	    << m_vsr[i].Vmp << ", " << m_vsr[i].Tkel << ", " << m_vsr[i].absA0 << ", "
	   << m_vsr[i].absB0 << ", " << m_vsr[i].Rs << ", " << m_vsr[i].eta
	   << ", " << m_vsr[i].Rh << ", " << m_vsr[i].I0 << ", " << m_vsr[i].Iph << ", "
	   << m_vsr[i].Emin <<std::endl; 
    }
    fs.close();
    printf("SearchHistoryAM::saveFile(), %zd lines added to file '%s'.\n", m_vsr.size()-m_fileLineNum, m_fileName.c_str());///////////////////
}
