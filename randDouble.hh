/**************************************************************************
 *        randDouble.hh
 ************************************************************************/
#if !defined(__RANDDOUBLE_HH)
#define __RANDDOUBLE_HH
#include <random>
#include <functional>

class RandDouble{
public:
    RandDouble(double low, double high)
	:r(std::bind(std::uniform_real_distribution<double>(low, high), std::default_random_engine())){
    }
    double operator()(){ return r(); } 

private:
    function<double()> r;
};

#endif
