/*******************************************************************
 *      LambertWhara.hh
 ******************************************************************/
#ifndef _utl_LambertWhara_h_
#define _utl_LambertWhara_h_

#include <cmath>
#include "HornerHara.hh"

namespace utl {
    double LambertW0haraAexp(const double a, const double x);
    
    template<typename Double, int n> struct Pade {
	static inline Double Approximation(const Double x);
    };

    template<typename Double> struct Pade<Double, 1>{
	static inline Double Approximation(const Double x){
	    return x * HORNER4(Double, x, 0.07066247420543414, 2.4326814530577687, 6.39672835731526, 4.663365025836821, 0.99999908757381) /
		HORNER4(Double, x, 1.2906660139511692, 7.164571775410987, 10.559985088953114, 5.66336307375819, 1);
	}
    };

    template<typename Double> struct Pade<Double, 2>{
	static inline Double Approximation(const Double x){
	    const Double y = std::log(Double(0.5)*x) - 2;
	    return 2 + y * HORNER3(Double, y, 0.00006979269679670452, 0.017110368846615806, 0.19338607770900237, 0.6666648896499793) /
		HORNER2(Double, y, 0.0188060684652668, 0.23451269827133317, 1);
	}
    };

    class BranchPoint { };

    template<typename Double> class Branch{
    public:
	template<int order> static Double AsymptoticExpansionAexp(const Double a, const Double x);
    };

    template<typename Double> inline Double HalleyStep(const Double x, const Double w);

    template<typename Double> inline Double HalleyStepAexp(const Double a, const Double y, const Double w);

    template<typename Double, int order> inline Double AsymptoticExpansionImpl(const Double a, const Double b);
    
    template<unsigned int order> class AsymptoticPolynomialB { };
    class AsymptoticPolynomialA { };
}

#endif
