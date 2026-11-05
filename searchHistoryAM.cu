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

/***************** Constructor *****************************/
SearchHistoryLM::SearchHistoryLM(){
    readFile();
}

/***************** Destructor *****************************/
SearchHistoryAM::~SearchHistoryAM(){
    saveFile();
}

/***************** Destructor *****************************/
SearchHistoryLM::~SearchHistoryLM(){
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
void SearchHistoryLM::pushSearchResult(SearchResultLM &sr){
    bool be_equal {false};
    int i1, i2;
    for(int i ; i<m_vsr.size(); i++){
	if(m_vsr[i].g == sr.g){
	    i1 = i;
	    be_equal = true;
	    for(int j=i+1; j<m_vsr.size(); j++){
		if(sr.g > m_vsr[j].g){
		    i2 = j-1;
		    break;
		}
	    }
	}
	if(sr.g > m_vsr[i].g){
	    i2 = i;
	    if(i==0) i1 = 0;
	    else i1 = i-1;
	    break;
	}
    }  // From here !! 2026.11.3

}

/**********************************************************/
int SearchHistoryLM::nearestSearchResult(SearchResultLM &sr){
    for(int i=0; i<m_vsr.size(); i++){
	// if(m_vsr[i].equal(sr)) return i;
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
void SearchHistoryLM::readFile(){
    std::fstream fs {m_fileName, std::ios_base::in};
    if(!fs){
	std::cout << "couldn't open '"<< m_fileName << "' for reading" << std::endl;
	return;
    }

    std::string str;
    while(std::getline(fs, str)){
	std::vector<std::string> vStr = PCSData::splitString(str, ',');

	SearchResultLM sr;
	sr.g = stod(vStr[0]);
	sr.t = stod(vStr[1]);
	sr.rVocSTC = stod(vStr[2]);
	sr.rIscSTC = stod(vStr[3]);
	sr.rVmpSTC = stod(vStr[4]);
	sr.rImpSTC = stod(vStr[5]);
	sr.Rs = stod(vStr[6]);
	sr.eta = stod(vStr[7]);
	sr.Iph = stod(vStr[8]);
	sr.I0 = stod(vStr[9]);
	sr.Rh = stod(vStr[10]);
	sr.Emin = stod(vStr[11]);

	m_vsr.push_back(sr);
	m_fileLineNum++;
    }
    printf("SearchHistoryLM::readFile(), file '%s' read.\n",
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

/***************************************************************/
void SearchHistoryLM::saveFile(){
    std::fstream fs {m_fileName, std::ios_base::app};
    fs << std::scientific << std::setprecision(20);
	
    for(int i=m_fileLineNum; i<m_vsr.size(); i++){
	fs << m_vsr[i].g << ", " << m_vsr[i].t << ", " << m_vsr[i].rVocSTC << ", "
	    << m_vsr[i].rIscSTC << ", " << m_vsr[i].rVmpSTC << ", "
	   << m_vsr[i].rImpSTC << ", "
	   << m_vsr[i].Rs << ", " << m_vsr[i].eta << ", " << m_vsr[i].Iph
	   << ", " << m_vsr[i].I0 << ", " << m_vsr[i].Rh << ", "
	   << m_vsr[i].Emin <<std::endl; 
    }
    fs.close();
    printf("SearchHistoryAM::saveFile(), %zd lines added to file '%s'.\n", m_vsr.size()-m_fileLineNum, m_fileName.c_str());///////////////////
}
