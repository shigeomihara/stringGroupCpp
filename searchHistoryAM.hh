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

// /***********************************************/
// struct SearchResultLM{
//     // inputs
//     double rVocSTC, rIscSTC, rVmpSTC, rImpSTC, g, t;
//     // outputs
//     double Rs, eta, Iph, I0, Rh, Emin;

//     double absDist(SearchResultLM &sr);
// };

/************************************************/
class SearchHistoryAM{
private:
    std::string m_fileName;
    int m_fileLineNum {0};

    float m_IscSTC_MD;
    float m_VocSTC_MD;
    float m_ImpSTC_MD;
    float m_VmpSTC_MD;

public:
    SearchHistoryAM(bool haveGPU);
    ~SearchHistoryAM();

    void setSTCValues(float IscSTC, float VocSTC, float ImpSTC, float VmpSTC); 
    
    // A history of search results
    std::vector<SearchResult> m_vsr;
    
    int haveSearchResult(SearchResult &sr);
    double absDist(SearchResult &sr1, SearchResult &sr2);
    int nearestSearchResult(SearchResult &sr, double &minAbsDist);
    
    void readFile();
    void saveFile();
};

// /******* 2026.11.3--
// 	 Search History for LM algorithm for adaptive model
// *************************************/
// class SearchHistoryLM{
// private:
//     std::string m_fileName {"Dat/searchHistoryLM.dat"};
//     int m_fileLineNum {0};

// public:
//     SearchHistoryLM();
//     ~SearchHistoryLM();
    
//     // A history of search results
//     std::vector<SearchResultLM> m_vsr;

//     int nearestSearchResult(SearchResultLM &sr, double &minAbsDist);
//     void addSearchResult(SearchResultLM &sr);
//     void readFile();
//     void saveFile();
// };

#endif

