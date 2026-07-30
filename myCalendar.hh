/**************************************************************************
 *        myCalendar.hh
 ************************************************************************/
#if !defined(__MYCALENDAR_HH)
#define __MYCALENDAR_HH
#include <ctime>
#include "PCSData.hh"

class MyCalendar: public tm{
public:
    MyCalendar(){}
    void set(std::string timeStr);

};

void MyCalendar::set(std::string timeStr){
    printf("%s\n", timeStr.c_str());
    std::vector<std::string> splitStr = PCSData::splitString(timeStr, '/');
    for(int i=0; i<splitStr.size(); i++){
	printf("%s\n", splitStr[i].c_str());
    }

    // int year {stoi(splitStr[0])};
    tm_year = stoi(splitStr[0])-1900;
    tm_mon = stoi(splitStr[1])-1;  // 1月が0
    tm_mday = stoi(splitStr[2]);
    std::vector<std::string> splitStrTime = PCSData::splitString(splitStr[3], ':');
    tm_hour = stoi(splitStrTime[0]);
    tm_min = stoi(splitStrTime[1]);
    printf("year=%d, month=%d, day=%d, hour=%d, min=%d\n",
	   tm_year, tm_mon, tm_mday, tm_hour, tm_min);
    char str[50];
    // std::tm time;
    strftime(str, 50, "%Y/%m/%d %H:%M", this);
    std::cout << str << std::endl;
    exit(0);
    
}

#endif
