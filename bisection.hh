/*************************************************************************************
 *           bisection.hh
 ************************************************************************************/
#if !defined(__BISECTION_HH)
#define __BISECTION_HH
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>   // 2026.4.13

template<typename T>
class Bisection{
    T* p;
public:
    Bisection(T* pp): p{pp} { };
    
    /** ÇQï™ñ@Ç…ÇÊÇÈóÎì_íTçı  x1 < x2,  f(x1)*f(x2)<0 Ç≈Ç†ÇÈÇ±Ç∆
     *  fÇ™É[ÉçÇ…Ç»ÇÈxÇï‘Ç∑ */
    double search(double x1, double x2);
    double search2(double x1, double x2);
};

/** ÇQï™ñ@Ç…ÇÊÇÈóÎì_íTçı  x1 < x2,  f(x1)*f(x2)<0 Ç≈Ç†ÇÈÇ±Ç∆  ***********************
 *  fÇ™É[ÉçÇ…Ç»ÇÈxÇï‘Ç∑   **********************************************************/
template<typename T> double Bisection<T>::search(double x1, double x2){
    double xLimit = std::abs(x1)<std::abs(x2) ? 1.0E-12*std::abs(x2) : 1.0E-12*std::abs(x1);
    double f1 {p->func(x1)};
    double f2 {p->func(x2)};
	
    // int count {0};////////////////////////////
    if(f1<0.0 && f2>0.0){
	while(true){
	    double xmid = (x1+x2)*0.5;
	    double fmid = p->func(xmid);
	    if(fmid <= 0.0) x1 = xmid;
	    else x2 = xmid;
	    // if( std::abs(fmid)<1.0E-8  || x2-x1 < xLimit ) return xmid;
	    if( x2-x1 < xLimit ) return xmid;
	    // if(count%10==0) printf("count=%d, x2-x1=%g, xLimit=%g\n", count, x2-x1, xLimit);///////////////////
	    // count++;////////////////////////////////
	}
    }else if(f1>0.0 && f2<0.0){
	while(true){
	    double xmid = (x1+x2)*0.5;
	    double fmid = p->func(xmid);
	    if(fmid >= 0.0) x1 = xmid;
	    else x2 = xmid;
	    // if( std::abs(fmid)<1.0E-8 || x2-x1 < xLimit ) return xmid;
	    if( x2-x1 < xLimit ) return xmid;
	    // if(count%10==0) printf("count=%d, x2-x1=%g, xLimit=%g\n", count, x2-x1, xLimit);///////////////////
	    // count++;
	}
    }else if(f1==0.0){
	return x1;
    }else if(f2==0.0){
	return x2;
    }else{
	return std::numeric_limits<double>::quiet_NaN();
    }
}

/***************** func2 *************************************************/
template<typename T> double Bisection<T>::search2(double x1, double x2){
    double xLimit = std::abs(x1)<std::abs(x2) ? 1.0E-12*std::abs(x2) : 1.0E-12*std::abs(x1);
    double f1 {p->func2(x1)};
    double f2 {p->func2(x2)};
	
    if(f1<0.0 && f2>0.0){
	while(true){
	    double xmid = (x1+x2)*0.5;
	    double fmid = p->func2(xmid);
	    if(fmid <= 0.0) x1 = xmid;
	    else x2 = xmid;
	    // if( std::abs(fmid)<1.0E-8  || x2-x1 < xLimit ) return xmid;
	    if( x2-x1 < xLimit ) return xmid;
	}
    }else if(f1>0.0 && f2<0.0){
	while(true){
	    double xmid = (x1+x2)*0.5;
	    double fmid = p->func2(xmid);
	    if(fmid >= 0.0) x1 = xmid;
	    else x2 = xmid;
	    // if( std::abs(fmid)<1.0E-8 || x2-x1 < xLimit ) return xmid;
	    if( x2-x1 < xLimit ) return xmid;
	}
    }else if(f1==0.0){
	return x1;
    }else if(f2==0.0){
	return x2;
    }else{
	return std::numeric_limits<double>::quiet_NaN();
    }
}

#endif
