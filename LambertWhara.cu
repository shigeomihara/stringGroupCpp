// -*- C++ -*-
/************************************************************************
 *             LambertWhara.cu
 ************************************************************************/
#include "HornerHara.hh"
#include "LambertWhara.hh"
#include <iostream> ///////////count
#include <cmath>

namespace utl {
    double log1p38 = std::log(1.38);
    double log236 = std::log(236.0);
    
    // S.Hara Feb. 18, 2024--
    // return W(a*exp(x)),  a>0
    double LambertW0haraAexp(const double a, const double x){
	if(a==0.0) return 0.0;     // S.Hara  May 6, 2024
	
	double w {0.0};
	double wNext {0.0};

	double xPlusLogA {x + std::log(a)};
	double aExpX;
	
	if(xPlusLogA < log1p38){
	    aExpX = a*std::exp(x);
	    if(aExpX == 0.0) return 0.0;    // S.Hara
	    w = Pade<double, 1>::Approximation(aExpX);
	}else if(xPlusLogA < log236){
	    aExpX = a*std::exp(x);
	    w = Pade<double, 2>::Approximation(aExpX);
	}else{
	    w = Branch<double>::AsymptoticExpansionAexp<5>(a, x);
	}

	// if(!isfinite(w)){//////////////t
	//     printf("a=%g, x=%g\n", a, x);
	//     std::cout << "initial w="<<w<<std::endl;
	//     exit(0);
	// }
	    
	// int count = 0;////////////////////t
	
	if(xPlusLogA < log236){
	    while(true){
		wNext = HalleyStep<double>(aExpX, w);
		// if(std::abs(wNext-w)<1.0E-12) break;
		if(std::abs(wNext-w)<=1.0E-10*std::abs(w)) break;

		// count++;////////////////////////////t
		// if(count%1000==0){
		//     std::cout << "LambertW0haraAexp-0, count="<< count << ", w=" << w << ", wNext=" << wNext << std::endl;////////////////t
		// }
		
		w = wNext;
	    }
	}else{
	    while(true){
		wNext = HalleyStepAexp<double>(a, x, w);
		// if(isnan(wNext)){////////////////////t
		//     printf("a=%g, x=%g, count=%d\n", a, x, count);
		//     std::cout << "w=" << w << ", wNext=" << wNext << std::endl;
		//     exit(0);
		// }
		
		// if(std::abs(wNext-w)<1.0E-10) break;
		if(std::abs(wNext-w)<=1.0E-10*std::abs(w)) break;
		
		// count++;////////////////////////////t
		// if(count%1000==0){
		//     std::cout << "LambertW0haraAexp-1, count="<< count << ", w=" << w << ", wNext=" << wNext << std::endl;////////////////t
		// }
		
		w = wNext;
	    }
	}
	
	return wNext;
    }

    // SH  Feb. 18, 2024
    template<typename Double> template<int order> Double Branch<Double>::AsymptoticExpansionAexp(const Double a, const Double x)
    {
	const Double logsx = std::log(a)+x;
	const Double logslogsx = std::log(logsx);
	return AsymptoticExpansionImpl<Double, order>(logsx, logslogsx);
    }
    
    template<typename Double, int order> inline Double AsymptoticExpansionImpl(const Double a, const Double b){
	return a + Horner<Double, AsymptoticPolynomialA, order>::Eval(1/a, b);
    }

    template<typename Double> inline Double HalleyStep(const Double x, const Double w){
	const Double ew = std::exp(w);
	const Double wew = w * ew;
	const Double wewx = wew - x;
	const Double w1 = w + 1;
	return w - wewx / (ew * w1 - (w + 2) * wewx/(2*w1));
    }

    // S. Hara Mar. 8, 2024
    template<typename Double> inline Double HalleyStepAexp(const Double a, const Double y, const Double w){
	const Double ewy = std::exp(w-y);
	const Double wewy = w * ewy;
	const Double wewya = wewy - a;
	const Double w1 = w + 1;
	return w - wewya / (ewy * w1 - (w + 2) * wewya/(2*w1));
    }
    
    template<typename Double> struct Polynomial<Double, AsymptoticPolynomialB<2>, 0>{
	static Double Coeff(){ return Double(0); }
    };
    template<typename Double> struct Polynomial<Double, AsymptoticPolynomialB<2>, 1>{
	static Double Coeff(){ return Double(-1); }
    };
    template<typename Double> struct Polynomial<Double, AsymptoticPolynomialB<2>, 2>{
	static Double Coeff(){ return Double(1./2); }
    };
    template<typename Double> struct Polynomial<Double, AsymptoticPolynomialB<3>, 0>{
	static Double Coeff(){ return Double(0); }
    };
    template<typename Double> struct Polynomial<Double, AsymptoticPolynomialB<3>, 1>{
	static Double Coeff(){ return Double(1); }
    };
    template<typename Double> struct Polynomial<Double, AsymptoticPolynomialB<3>, 2>{
	static Double Coeff(){ return Double(-3./2); }
    };
    template<typename Double> struct Polynomial<Double, AsymptoticPolynomialB<3>, 3>{
	static Double Coeff(){ return Double(1./3); }
    };
    template<typename Double> struct Polynomial<Double, AsymptoticPolynomialB<4>, 0>{
	static Double Coeff(){ return Double(0); }
    };
    template<typename Double> struct Polynomial<Double, AsymptoticPolynomialB<4>, 1>{
	static Double Coeff(){ return Double(-1); }
    };
    template<typename Double> struct Polynomial<Double, AsymptoticPolynomialB<4>, 2>{
	static Double Coeff(){ return Double(3); }
    };
    template<typename Double> struct Polynomial<Double, AsymptoticPolynomialB<4>, 3>{
	static Double Coeff(){ return Double(-11./6); }
    };
    template<typename Double> struct Polynomial<Double, AsymptoticPolynomialB<4>, 4>{
	static Double Coeff(){ return Double(1./4); }
    };
    template<typename Double> struct Polynomial<Double, AsymptoticPolynomialB<5>, 0>{
	static Double Coeff(){ return Double(0); }
    };
    template<typename Double> struct Polynomial<Double, AsymptoticPolynomialB<5>, 1>{
	static Double Coeff(){ return Double(1); }
    };
    template<typename Double> struct Polynomial<Double, AsymptoticPolynomialB<5>, 2>{
	static Double Coeff(){ return Double(-5); }
    };
    template<typename Double> struct Polynomial<Double, AsymptoticPolynomialB<5>, 3>{
	static Double Coeff(){ return Double(35./6); }
    };
    template<typename Double> struct Polynomial<Double, AsymptoticPolynomialB<5>, 4>{
	static Double Coeff(){ return Double(-25./12); }
    };
    template<typename Double> struct Polynomial<Double, AsymptoticPolynomialB<5>, 5>{
	static Double Coeff(){ return Double(-1./5); }
    };

    template<typename Double> struct Polynomial<Double, AsymptoticPolynomialA, 0> {
	static Double Coeff(const Double y) { return Double(-y); }
    };
    template<typename Double> struct Polynomial<Double, AsymptoticPolynomialA, 1> {
	static Double Coeff(const Double y) { return Double(y); }
    };
    template<typename Double> struct Polynomial<Double, AsymptoticPolynomialA, 2> {
	static Double Coeff(const Double y) { return Double(Horner<Double,  AsymptoticPolynomialB<2>, 2>::Eval(y)); }
    };
    template<typename Double> struct Polynomial<Double, AsymptoticPolynomialA, 3> {
	static Double Coeff(const Double y) { return Double(Horner<Double,  AsymptoticPolynomialB<3>, 3>::Eval(y)); }
    };
    template<typename Double> struct Polynomial<Double, AsymptoticPolynomialA, 4> {
	static Double Coeff(const Double y) { return Double(Horner<Double,  AsymptoticPolynomialB<4>, 4>::Eval(y)); }
    };
    template<typename Double> struct Polynomial<Double, AsymptoticPolynomialA, 5> {
	static Double Coeff(const Double y) { return Double(Horner<Double,  AsymptoticPolynomialB<5>, 5>::Eval(y)); }
    };

}
