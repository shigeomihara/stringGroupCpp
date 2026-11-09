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
    // rVocSTC
    int i1=-1, i2=0; // sr.rVocSTC==m_vsr[i].rVocSTC for i in [i1, i2)
    for(int i=0; i<m_vsr.size(); i++){
	if(m_vsr[i].rVocSTC == sr.rVocSTC){
	    i1 = i;
            int j {i+1};
	    for( ; j<m_vsr.size(); j++){
		if(sr.rVocSTC > m_vsr[j].rVocSTC) break;
	    }
            i2 = j;
            break;
	}
	if(sr.rVocSTC > m_vsr[i].rVocSTC){
            m_vsr.insert(m_vsr.begin()+i, sr);
	    return;
	}
    }
    if(i1 < 0){
        m_vsr.push_back(sr); 
        return;
    }

    // rIscSTC
    int i3=-1, i4=0; // sr.rIscSTC==m_vsr[i].rIscSTC for i in [i3, i4)
    for(int i=i1; i<i2; i++){
	if(m_vsr[i].rIscSTC == sr.rIscSTC){
	    i3 = i;
            int j {i+1};
	    for( ; j<i2; j++){
		if(sr.rIscSTC > m_vsr[j].rIscSTC) break;
	    }
            i4 = j;
            break;
	}
	if(sr.rIscSTC > m_vsr[i].rIscSTC){
            m_vsr.insert(m_vsr.begin()+i, sr);
	    return;
	}
    }
    if(i3 < 0){
        m_vsr.insert(m_vsr.begin()+i2, sr);
        return;
    }
    
    // rVmpSTC
    int i5=-1, i6=0; // sr.rVmpSTC==m_vsr[i].rVmpSTC for i in [i5, i6)
    for(int i=i3; i<i4; i++){
	if(m_vsr[i].rVmpSTC == sr.rVmpSTC){
	    i5 = i;
            int j {i+1};
	    for( ; j<i4; j++){
		if(sr.rVmpSTC > m_vsr[j].rVmpSTC) break;
	    }
            i6 = j;
            break;
	}
	if(sr.rVmpSTC > m_vsr[i].rVmpSTC){
            m_vsr.insert(m_vsr.begin()+i, sr);
	    return;
	}
    }
    if(i5 < 0){
        m_vsr.insert(m_vsr.begin()+i4, sr);
        return;
    }
    
    // rImpSTC
    int i7=-1, i8=0; // sr.rImpSTC==m_vsr[i].rImpSTC for i in [i7, i8)
    for(int i=i5; i<i6; i++){
	if(m_vsr[i].rImpSTC == sr.rImpSTC){
	    i7 = i;
            int j {i+1};
	    for( ; j<i6; j++){
		if(sr.rImpSTC > m_vsr[j].rImpSTC) break;
	    }
            i8 = j;
            break;
	}
	if(sr.rImpSTC > m_vsr[i].rImpSTC){
            m_vsr.insert(m_vsr.begin()+i, sr);
	    return;
	}
    }
    if(i7 < 0){
        m_vsr.insert(m_vsr.begin()+i6, sr);
        return;
    }
    
    // g
    int i9=-1, i10=0; // sr.g==m_vsr[i].g for i in [i9, i10)
    for(int i=i7; i<i8; i++){
	if(m_vsr[i].g == sr.g){
	    i9 = i;
            int j {i+1};
	    for( ; j<i8; j++){
		if(sr.g > m_vsr[j].g) break;
	    }
            i10 = j;
            break;
	}
	if(sr.g > m_vsr[i].g){
            m_vsr.insert(m_vsr.begin()+i, sr);
	    return;
	}
    }
    if(i9 < 0){
        m_vsr.insert(m_vsr.begin()+i8, sr);
        return;
    }
    
    // t
    int i11=-1, i12=0; // sr.t==m_vsr[i].t for i in [i11, i12)
    for(int i=i9; i<i10; i++){
	if(m_vsr[i].t == sr.t){
	    // All rVocSTC, rIscSTC, rVmpSTC, rImpSTC, g, t are the same as sr
	    // ==> don't add 
            return;
	}
	if(sr.t > m_vsr[i].t){
            m_vsr.insert(m_vsr.begin()+i, sr);
	    return;
	}
    }
    if(i11 < 0){
        m_vsr.insert(m_vsr.begin()+i10, sr);
        return;
    }
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
	sr.rVocSTC = stod(vStr[0]);
	sr.rIscSTC = stod(vStr[1]);
	sr.rVmpSTC = stod(vStr[2]);
	sr.rImpSTC = stod(vStr[3]);
	sr.g = stod(vStr[4]);
	sr.t = stod(vStr[5]);
	sr.Rs = stod(vStr[6]);
	sr.eta = stod(vStr[7]);
	sr.Iph = stod(vStr[8]);
	sr.I0 = stod(vStr[9]);
	sr.Rh = stod(vStr[10]);
	sr.Emin = stod(vStr[11]);

	m_vsr.push_back(sr);
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
    std::fstream fs {m_fileName, std::ios_base::out};  // out
    fs << std::scientific << std::setprecision(20);
	
    for(int i=0; i<m_vsr.size(); i++){  // from i=0
	fs << m_vsr[i].rVocSTC << ", " << m_vsr[i].rIscSTC << ", " << m_vsr[i].rVmpSTC << ", "
	   << m_vsr[i].rImpSTC << ", " << m_vsr[i].g << ", " << m_vsr[i].t << ", "
	   << m_vsr[i].Rs << ", " << m_vsr[i].eta << ", " << m_vsr[i].Iph
	   << ", " << m_vsr[i].I0 << ", " << m_vsr[i].Rh << ", "
	   << m_vsr[i].Emin <<std::endl; 
    }
    fs.close();
    printf("SearchHistoryLM::saveFile(), %zd lines saved to file '%s'.\n", m_vsr.size(), m_fileName.c_str());///////////////////
}
