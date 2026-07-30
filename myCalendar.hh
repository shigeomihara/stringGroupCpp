/**************************************************************************
 *        myCalendar.hh
 ************************************************************************/
#if !defined(__MYCALENDAR_HH)
#define __MYCALENDAR_HH
#include <ctime>
#include <string>

class MyCalendar: public tm{
public:
    MyCalendar(){}
    MyCalendar(std::string timeStr);
    void set(std::string timeStr);
    void set(int y, int m, int d, int h, int min);
    void set(int h, int min);
    void set(MyCalendar myCal);
    void addHour(double hour);
    void addMin(int min);
    bool before(MyCalendar myCal);
    bool after(MyCalendar myCal);
    void print();
};

#endif
