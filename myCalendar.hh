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

};

#endif
