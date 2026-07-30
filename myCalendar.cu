// -*- C++ -*-
/**************************************************************************
 *        myCalendar.cu
 ************************************************************************/
#include "myCalendar.hh"
#include "PCSData.hh"

MyCalendar::MyCalendar(std::string timeStr){
    set(timeStr);
}
    
void MyCalendar::set(std::string timeStr){
    // printf("%s\n", timeStr.c_str());
    std::vector<std::string> splitStr = PCSData::splitString(timeStr, '/');
    // for(int i=0; i<splitStr.size(); i++){
    //     printf("%s\n", splitStr[i].c_str());
    // }

    tm_year = stoi(splitStr[0])-1900;
    tm_mon = stoi(splitStr[1])-1;  // 1月が0
    tm_mday = stoi(splitStr[2]);
    std::vector<std::string> splitStrTime = PCSData::splitString(splitStr[3], ':');
    tm_hour = stoi(splitStrTime[0]);
    tm_min = stoi(splitStrTime[1]);
    // printf("year=%d, month=%d, day=%d, hour=%d, min=%d\n",
    //        tm_year, tm_mon, tm_mday, tm_hour, tm_min);
    // char str[50];
    // strftime(str, 50, "%Y/%m/%d %H:%M", this);
    // std::cout << str << std::endl;
    // exit(0);    
}

void MyCalendar::set(int y, int m, int d, int h, int min){
    tm_year = y-1900;
    tm_mon = m-1;  // 1月が0
    tm_mday = d;
    tm_hour = h;
    tm_min = min;
}

void MyCalendar::set(int h, int min){
    tm_hour = h;
    tm_min = min;
}

void MyCalendar::set(MyCalendar myCal){
    set(myCal.tm_year+1900, myCal.tm_mon+1, myCal.tm_mday, myCal.tm_hour,
        myCal.tm_min);
}

void MyCalendar::addHour(double hour){
     int h = static_cast<int>(hour);
     int min = static_cast<int>((hour-h)*60.0);
     // printf("hour=%f, h=%d, min=%d\n", hour, h, min);/////////
     tm_min += min;
     if(tm_min > 60){
         tm_hour += 1;
         tm_min -= 60;
     }
     tm_hour += h;
}

void MyCalendar::addMin(int min){
     tm_min += min;
     if(tm_min > 60){
         tm_hour += 1;
         tm_min -= 60;
     }
}

bool MyCalendar::before(MyCalendar myCal){
    if(tm_year < myCal.tm_year) return true;
    if(tm_year > myCal.tm_year) return false;
    if(tm_mon < myCal.tm_mon) return true;
    if(tm_mon > myCal.tm_mon) return false;
    if(tm_mday < myCal.tm_mday) return true;
    if(tm_mday > myCal.tm_mday) return false;
    if(tm_hour < myCal.tm_hour) return true;
    if(tm_hour > myCal.tm_hour) return false;
    if(tm_min <= myCal.tm_min) return true;
    return false;
}

bool MyCalendar::after(MyCalendar myCal){
    if(tm_year > myCal.tm_year) return true;
    if(tm_year < myCal.tm_year) return false;
    if(tm_mon > myCal.tm_mon) return true;
    if(tm_mon < myCal.tm_mon) return false;
    if(tm_mday > myCal.tm_mday) return true;
    if(tm_mday < myCal.tm_mday) return false;
    if(tm_hour > myCal.tm_hour) return true;
    if(tm_hour < myCal.tm_hour) return false;
    if(tm_min >= myCal.tm_min) return true;
    return false;
}

void MyCalendar::print(){
    char str[50];
    strftime(str, 50, "%Y/%m/%d %H:%M", this);
    std::cout << str << std::endl;
}
