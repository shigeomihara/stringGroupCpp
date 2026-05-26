/*********************************************************************
 *      searchHistoryAM.hh
 ********************************************************************/
#if !defined(_SEARCH_HISTORY_AM_HH)
#define _SEARCH_HISTORY_AM_HH

#include <vector>
#include <string>

struct SearchResult{
    // inputs
    double Isc, Voc, Imp, Vmp, Tkel, absA0, absB0;
    // outputs
    double Rs, eta, Rh, I0, Iph, Emin;

    bool equal(SearchResult &sr);
};

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

#endif

