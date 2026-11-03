/*********************************************************************
 *      searchHistoryAM.hh
 ********************************************************************/
#if !defined(_SEARCH_HISTORY_AM_HH)
#define _SEARCH_HISTORY_AM_HH

#include <vector>
#include <string>

/***********************************************/
struct SearchResult{
    // inputs
    double Isc, Voc, Imp, Vmp, Tkel, absA0, absB0;
    // outputs
    double Rs, eta, Rh, I0, Iph, Emin;

    bool equal(SearchResult &sr);
};

/***********************************************/
struct SearchResultLM{
    // inputs
    double g, t, rVocSTC, rIscSTC, rVmpSTC, rImpSTC;
    // outputs
    double Rs, eta, Iph, I0, Rh, Emin;
};

/************************************************/
class SearchHistoryAM{
private:
    std::string m_fileName;
    int m_fileLineNum {0};

public:
    SearchHistoryAM(bool haveGPU);
    ~SearchHistoryAM();
    
    // A history of search results
    std::vector<SearchResult> m_vsr;
    
    int haveSearchResult(SearchResult &sr);
    void readFile();
    void saveFile();
};

/******* 2026.11.3--
	 Search History for LM algorithm for adaptive model
*************************************/
class SearchHistoryLM{
private:
    std::string m_fileName {"Dat/searchHistoryLM.dat"};
    int m_fileLineNum {0};

public:
    SearchHistoryLM();
    ~SearchHistoryLM();
    
    // A history of search results
    std::vector<SearchResultLM> m_vsr;

    void pushSearchResult(SearchResultLM &sr);
    int nearestSearchResult(SearchResultLM &sr);
    void readFile();
    void saveFile();
};

#endif

