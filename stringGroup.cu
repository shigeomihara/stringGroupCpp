// -*- C++ -*-
/************************************************************************************
 *                                  stringGroup.cu
 ***********************************************************************************/
#include "globalConstants.hh"
#include "stringGroup.hh"
#include "bisection.hh"
#include <iostream>  /////////////t
#include <iomanip>  /////////////t
#include <random>
#include <string>

/** Constructor     *****************************************************************
 * GG (W/m^2), TT (degC)     *******************************************************/
StringGroup::StringGroup(GlobalConstants &gcc, double GG, double TT)
    : m_gc{gcc}, m_G{GG}, m_T{TT}, m_smInit{gcc, GG, TT}{
    m_Ib0 = m_smInit.getIb0();
    m_Ib0sg = gcc.S*m_Ib0;

    m_Vmin = m_smInit.getVmin();
    m_Vmax = m_smInit.getVmax();
    
    m_Imin = getI(m_Vmax);
    m_Imax = getI(m_Vmin);
}

/** Calculate stringGroup current I (A) for voltage V (V)    ************************
 *  入力： m_Vmin, m_Vmax, m_Ib0sg, m_vsm
 *  関数呼び出し： m_smInit,getI, m_vsm[i].getI      *******************************/
double StringGroup::getI(double V){
    if(V < m_Vmin) return std::numeric_limits<double>::infinity();
    else if(m_Vmax < V) return -m_Ib0sg;
    else{
	double I = (m_gc.S-m_vsmNumSum) * m_smInit.getI(V);
	for(int i=0; i<m_vsm.size(); i++){
	    I += m_vsm[i].getI(V)*m_vsmNum[i];
	}
	return I;
    }
}

/** m_aVsg[K]に対するIkを計算して、Ik[K]に書き込む    *******************************/
void StringGroup::getI(double Ik[GlobalConstants::K]){
    for(int k=0; k<m_gc.K; k++){
	Ik[k] = 0.0;
	for(int s=0; s<m_vsm.size(); s++) Ik[k] += m_vsm[s].get_aIk(k);
	Ik[k] += m_smInit.get_aIk(k)*(m_gc.S - m_vsm.size());
    }
}

/** Calculate stringGroup voltage V (V) for current I (A)    ************************
 *  入力： m_Imax, m_Vmin        ***************************************************/
double StringGroup::getV(double I){
    if(m_Imax < I) return m_Vmin;
    else if(-m_Ib0sg <= I && I < m_Imin) return std::numeric_limits<double>::infinity();
    else if(I < -m_Ib0sg) return  std::numeric_limits<double>::quiet_NaN();
    else{
	Bisection<StringGroup> bis {this};
	m_Ifunc = I;
	return bis.search(m_Vmin, m_Vmax);
    }
    return 0.0;
}

/***************     getIで使われる２分法用のメソッド      **************************/
double StringGroup::func(double V){
    return getI(V)-m_Ifunc;
}
/************************************************************************/
void StringGroup::addSM(StringModel &sm, int s){
    m_vsm.push_back(sm);
    m_vsmNum.push_back(s);
    m_vsmNumSum += s;

    // printf("Vmin=%g, Vmax=%g, Imin=%g, Imax=%g\n",
    // 	   m_Vmin, m_Vmax, m_Imin, m_Imax);////////////////
    double Vmax {sm.getVmax()};
    if(Vmax > m_Vmax){
	m_Vmax = Vmax;
	m_Imin = getI(m_Vmax);
    }
    double Vmin {sm.getVmin()};
    if(Vmin > m_Vmin){  // Vmin = max Vmin_i !!
	m_Vmin = Vmin;
	m_Imax = getI(m_Vmin);
    }
    // printf("Vmin=%g, Vmax=%g, Imin=%g, Imax=%g\n",
    // 	   m_Vmin, m_Vmax, m_Imin, m_Imax);////////////////
}

/******************************************************************/
void StringGroup::setGT(double G, double T){
    m_G = G;
    m_T = T;

    m_smInit.setGT(G, T);
    for(int i=0; i<m_vsm.size(); i++) m_vsm[i].setGT(G, T);

    m_Ib0 = m_smInit.getIb0();
    m_Ib0sg = m_gc.S*m_Ib0;

    m_Vmin = m_smInit.getVmin();
    m_Vmax = m_smInit.getVmax();
    for(int i=0; i<m_vsm.size(); i++){
	double Vmax {m_vsm[i].getVmax()};
	if(Vmax > m_Vmax) m_Vmax = Vmax;
	double Vmin {m_vsm[i].getVmin()};
	if(Vmin > m_Vmin) m_Vmin = Vmin;  // Vmin = max Vmin_i !!
    }	
    
    m_Imin = getI(m_Vmax);
    m_Imax = getI(m_Vmin);

}

/****** Suppose refreshment after set5params of adaptiveModel ***********/
void StringGroup::refreshMaxMin(){
    m_smInit.refreshMaxMin();
    for(int i=0; i<m_vsm.size(); i++) m_vsm[i].refreshMaxMin();

    m_Vmin = m_smInit.getVmin();
    m_Vmax = m_smInit.getVmax();
    for(int i=0; i<m_vsm.size(); i++){
	double Vmax {m_vsm[i].getVmax()};
	if(Vmax > m_Vmax) m_Vmax = Vmax;
	double Vmin {m_vsm[i].getVmin()};
	if(Vmin > m_Vmin) m_Vmin = Vmin;  // Vmin = max Vmin_i !!
    }	
    
    m_Imin = getI(m_Vmax);
    m_Imax = getI(m_Vmin);
}

/********* Return -P for searching maximum P ***************************/
double StringGroup::operator()(const double V){
    return -getI(V)*V;
}
/***************************************************************/
double StringGroup::getPmax(double &Vmp, double &Imp){
    double Voc {getV(0.0)};
    // printf("Voc=%g\n", Voc);//////////////////////////////

    // Golden golden;
    // golden.ax = 0.0;
    // golden.cx = Voc;
    // golden.bx = 0.38197*Voc;  // Golden section
    // Vmp = golden.minimize(*this);
    // double Pmax = -golden.fmin;

    // Brent's method
    Brent brent;
    brent.ax = 0.0;
    brent.cx = Voc;
    brent.bx = 0.38197*Voc;  // Golden section
    Vmp = brent.minimize(*this);
    double Pmax = -brent.fmin;

    Imp = getI(Vmp);
    return Pmax;
}

/** 測定データG, T, I, V に最もフィットするパラメータを探索して、それに設定する   ***
 *  出力： m_aG[], m_aT[], m_aIsg[], m_aVsg[]     **********************************/
void StringGroup::paramEstim(double GG[], double TT[], double II[], double VV[], std::string Time[]){
    for(int k=0; k<m_gc.K; k++){
	m_aG[k] = GG[k];
	m_aT[k] = TT[k];
	m_aIsg[k] = II[k];
	m_aVsg[k] = VV[k];
    }
    
    double Eret = m_smInit.setVars(GG, TT);
    if(Eret==std::numeric_limits<double>::infinity()){
	std::cout << "StringGroup::paramEstim, smInit.setVars()=INFINITY, error\n";
	exit(-1);
    }

    double IkStr[GlobalConstants::K];

    for(int k=0; k<m_gc.K; k++){
	m_smInit.getIjunbi(k);
	IkStr[k] = m_smInit.getI(VV[k]);
	// printf("Imeas[%d]=%g, I=%g\n", k, II[k], IkStr[k]*m_gc.S);//////////////////t 
    }

    m_smInit.set_aIk(IkStr);
    
    double bunsi {0.0};
    double bunbo {0.0};
    for(int k=0; k<m_gc.K; k++){
	double x = II[k]-IkStr[k]*m_gc.S;
	bunsi += x*x;
	bunbo += II[k]*II[k];
    }
    m_bunbo = bunbo;
    double E = sqrt(bunsi/bunbo);
    printf("Initial E=%g\n", E);//////////////////t
    
    // m_smInit.print(-1);//////////////////////t

    double EminPre = E;
    double EminPreString = E;

    double delta = 0.002;
    double deltaString = 0.001;
    bool finish = false; 
    bool finishString = false;

    for(int s=0; s < m_gc.S; s++){ // 各ストリングsに対して
	double Emin = std::numeric_limits<double>::infinity();
	
	StringModel sm_s {m_gc, 1000.0, 25.0};
	Eret = sm_s.setVars(GG, TT);
	if(Eret==std::numeric_limits<double>::infinity()){
	    std::cout << "StringGroup::paramEstim, sm_s.setVars()=INFINITY, error\n";
	    exit(-1);
	}
	
	m_vsm.push_back(sm_s);

	for(int l=0; l<m_gc.L; l++){ // 各サブストリングlに対して
	    std::cout << "(s, l)=(" << s << ',' << l << ")" << std::endl;
	    m_vsm[s].addSS();
	    // m_vsm[s].m_vss[l].setVars(GG, TT);////////// 2024/6/8
	    m_vsm[s].setVarsVss(GG, TT, l);  // 2024/6/9

	    double X[4];
	    Emin = ELPSOsearch(X, s, l);

	    m_vsm[s].setSS_STCparams(X, l);
	    m_vsm[s].setVarsVssGetE(l);

	    printf("Emin=%g, EminPre=%g, EminPreString=%g, delta_now=%g, delta=%g, deltaStr_now=%g, deltaString=%g, IscSTC=%g, VocSTC=%g, ImpSTC=%g, VmpSTC=%g\n",
		   Emin, EminPre, EminPreString, 1.0-Emin/EminPre, delta, 1.0-Emin/EminPreString, deltaString, X[0], X[1], X[2]*X[0], X[3]*X[1]);
	    
	    // FiveParamsPlus fpp = m_vsm[s].get5paramsPlus(l);
	    // fpp.print();

	    if(Emin >= (1.0-delta)*EminPre){
		if(Emin >= (1.0-deltaString)*EminPreString){
		    finish = true;
		}else{
		    EminPreString = Emin;
		    finishString = true;
		}
	    }

	    if(finish) break;

	    EminPre = Emin;
	    
	    if(finishString){
		finishString = false;
		break;
	    }
	}

	m_vsm[s].set_aIk_for(VV);
	
	double Ik[GlobalConstants::K];
	getI(Ik);
	
	for(int k=0; k<m_gc.K; k++){
	    m_vsm[s].getIjunbi(k);
	    IkStr[k] = m_vsm[s].getI(VV[k]);
	    std::cout << "s=" << s << " k=" << k << " IkStr=" << IkStr[k]
		      << " Ik=" << Ik[k] << " Imeas,k=" << m_aIsg[k] << std::endl;///////t 
	}

	m_vsm[s].set_aIk(IkStr);

	EminPreString = Emin;
	
	if(finish) break;
    }
}

/** ストリングsのサブストリングlのSTCパラメータIscSTC, VocSTC, ImpSTC, VmpSTCを   ***
 *  Xで決めて、誤差関数Eの値を計算して返す
 *  入力： m_aG[], m_aT[], m_aVsg[], m_gc, m_vsm, m_aIsg, m_smInit
 *  関数呼び出し： StringModel::setSS_STCparams, StringModel::setVars,
 *                 StringModel::getIjunbi, StringModel::getI,
 *                 StringModel.set_aIk, StringModel::get_aIk   **********************/
double StringGroup::getE(double X[4], int s, int l){
    m_vsm[s].setSS_STCparams(X, l);
    // m_vsm[s].m_vss[l].m_am.setSTCparams(X);  // 一行上よりかえって遅くなった 2024/6/8
    // if(m_vsm[s].setVars(m_aG, m_aT)==INFINITY) return INFINITY;
    // if(m_vsm[s].setVarsVss(m_aG, m_aT, l)==INFINITY) return INFINITY;
    // m_vsm[s].m_vss[l].m_am.calcPCS(m_aG, m_aT);
    if(m_vsm[s].setVarsVssGetE(l)==std::numeric_limits<double>::infinity()) return std::numeric_limits<double>::infinity();
    
    double IkStr[GlobalConstants::K];
    
    for(int k=0; k<m_gc.K; k++){
	// m_vsm[s].getIjunbi(k, l);   // これ、ダメ 2024.6.22
	m_vsm[s].getIjunbi(k);
	IkStr[k] = m_vsm[s].getI(m_aVsg[k]);
	// printf("Isg[%d]/S=%g, Istr=%g\n", k, m_aIsg[k]/m_gc.S, IkStr[k]);//////////////////t 
    }
    // m_vsm[s].getI(l, m_aVsg, IkStr);  // かえって遅くなった、不採用

    // m_vsm[s].set_aIk(IkStr);
    
    double bunsi {0.0};
    // double bunbo {0.0};
    for(int k=0; k<m_gc.K; k++){
	double x = m_aIsg[k]-IkStr[k];
	// bunbo += x*x;
	// for(int s1=0; s1<=s; s1++) x -= m_vsm[s1].get_aIk(k);
	for(int s1=0; s1<s; s1++) x -= m_vsm[s1].get_aIk(k);
	x -= m_smInit.get_aIk(k)*(m_gc.S-m_vsm.size());
	bunsi += x*x;
    }
    double E = sqrt(bunsi/m_bunbo);
    // std::cout << "E=" << E << " bunsi=" << bunsi << " bunbo=" << bunbo << std::endl;//////t
    //    printf("E=%g\n", E);//////////////////t

    return E;
}

/*************************     ******************************************************/
double StringGroup::ELPSOsearch(double XX[4], int ss, int ll){
    const int varNum {4}; 
    // const int Np {100};
    const int Np {70};

    double VocSTCss = m_gc.getVocSTCss();
    const double Xmin[] = {0.5*m_gc.IscSTC, 0.5*VocSTCss, 0.4, 0.5};
    const double Xmax[] = {1.2*m_gc.IscSTC, 1.2*VocSTCss, 0.98, 1.0};
    const double Xinit[] = {m_gc.IscSTC, VocSTCss, m_gc.ImpSTC/m_gc.IscSTC, m_gc.VmpSTC/m_gc.VocSTC};

    // printf("Xmin: %g %g %g %g\n", Xmin[0],  Xmin[1],  Xmin[2],  Xmin[3]); 
    // printf("Xmax: %g %g %g %g\n", Xmax[0],  Xmax[1],  Xmax[2],  Xmax[3]);
    // printf("Xinit: %g %g %g %g\n", Xinit[0],  Xinit[1],  Xinit[2]*8.61,  Xinit[3]*37.05);/////////t
    // exit(0);/////////////
    // return;
    
    double X[Np][varNum];
    double P[Np][varNum];  // Personal best
    double Pg[varNum];  // Global best

    double Eval_P[Np];
    for(int i=0; i<Np; i++) Eval_P[i] = std::numeric_limits<double>::infinity();
    double Eval_Pg {std::numeric_limits<double>::infinity()};
    
    double V[Np][varNum];
    for(int i=0; i<Np; i++)
	for(int j=0; j<varNum; j++) V[i][j] = 0.0;

    // int tMax = 1000;
    int tMax = 6;   ////////////////////////////t

    double w {0.9};
    double del_w {(0.9-0.4)/tMax};
    double C1 {2.0};
    double C2 {2.0};

    double a {0.5};
    double del_a = {(0.5-0.1)/tMax};
    double s {0.9};
    double del_s {(0.9-0.4)/tMax};

    std::random_device rd;
    std::default_random_engine engine {rd()};
    std::uniform_real_distribution<double> uniform {0.0, 1.0}; 

    // for(int i=0; i<10; i++){
    // 	std::cout << distribution(engine) << " ";
    // }
    // printf("\n");

    for(int k=0; k<varNum; k++) X[0][k] = Xinit[k];
    
    for(int i=1; i<Np; i++){
	for(int k=0; k<varNum; k++){
	    X[i][k] = Xmin[k] +(Xmax[k]-Xmin[k])*uniform(engine);
	}
    }

    // for(int i=0; i<Np; i++){
    // 	printf("##### X[%d][]\n", i);
    // 	for(int k=0; k<varNum; k++){
    // 	    printf("%g ", X[i][k]);
    // 	}
    // 	printf("\n");
    // }

    // for(int i=1; i<Np; i++)
    // 	for(int k=0; k<varNum; k++) P[i][k] = X[i][k];

    // for(int k=0; k<varNum; k++) Pg[k] = X[0][k];

    std::normal_distribution<double> normal {0.0, 1.0};
    std::cauchy_distribution<double> cauchy {0.0, 1.0};

    for(int t=1; ; t++){
	for(int i=0; i<Np; i++){
	    double Eval_Xi = getE(X[i], ss, ll);
	    // if(i==0){//////////////////////////t
	    // 	printf("t=%d, ss=%d, ll=%d, Eval_X0=%g\n", t, ss, ll, Eval_Xi);
	    // }

	    // double IkStr[GlobalConstants::K];////////////////////t
	    // for(int k=0; k<GlobalConstants::K; k++){
	    // 	// m_vsm[ss].getIjunbi(k, ll);
	    // 	m_vsm[ss].getIjunbi(k);
	    // 	IkStr[k] = m_vsm[ss].getI(m_aVsg[k]);
	    // }
	    // m_vsm[ss].set_aIk(IkStr);
	    // // m_vsm[ss].get5paramsPlus(ll).print();
	    // m_vsm[ss].print(ll);/////////////////////t
	    // exit(0);////////////////////////t

	    if(Eval_Xi < Eval_P[i]){
		Eval_P[i] = Eval_Xi;
		for(int k=0; k<varNum; k++) P[i][k] = X[i][k];
	    }
	    if(Eval_Xi < Eval_Pg){
		Eval_Pg = Eval_Xi;
		for(int k=0; k<varNum; k++) Pg[k] = X[i][k];
	    }
	}

	// // Show Eval_Pg
	// printf("# t=%d: Eval_Pg= %.6f,", t, Eval_Pg);
	// for(int k=0; k<varNum; k++)
	//     printf(" Pg[%d]= %.4f", k, Pg[k]);
	// printf("\n");

	if(t == tMax || Eval_Pg < 1.0E-5) break; // 収束条件

	for(int i=0; i<Np; i++){
	    for(int k=0; k<varNum; k++){
		V[i][k] = w*V[i][k] +C1*uniform(engine)*(P[i][k]-X[i][k])
		    +C2*uniform(engine)*(Pg[k]-X[i][k]);
	    }
	}

	for(int i=0; i<Np; i++){
	    for(int k=0; k<varNum; k++){
		X[i][k] += V[i][k];
		if(X[i][k] < Xmin[k]) X[i][k] = Xmin[k];
		if(X[i][k] > Xmax[k]) X[i][k] = Xmax[k];
	    }
	}

	// mutation of Pg
	double Pg1[varNum];
	for(int k=0; k<varNum; k++){
	    Pg1[k] = Pg[k] +(Xmax[k]-Xmin[k])*normal(engine)*a;
	    if(Pg1[k] < Xmin[k]) Pg1[k] = Xmin[k];
	    if(Pg1[k] > Xmax[k]) Pg1[k] = Xmax[k];
	}
	
	double Eval_Pg1 = getE(Pg1, ss, ll);
	if(Eval_Pg1 < Eval_Pg){
	    Eval_Pg = Eval_Pg1;
	    for(int k=0; k<varNum; k++) Pg[k] = Pg1[k];
	}

	for(int k=0; k<varNum; k++){
	    Pg1[k] = Pg[k] +(Xmax[k]-Xmin[k])*cauchy(engine)*s;
	    if(Pg1[k] < Xmin[k]) Pg1[k] = Xmin[k];
	    if(Pg1[k] > Xmax[k]) Pg1[k] = Xmax[k];
	}

	Eval_Pg1 = getE(Pg1, ss, ll);
	if(Eval_Pg1 < Eval_Pg){
	    Eval_Pg = Eval_Pg1;
	    for(int k=0; k<varNum; k++) Pg[k] = Pg1[k];
	}

	for(int k=0; k<varNum; k++){
	    for(int kk=0; kk<varNum; kk++){
		if(kk==k){
		    Pg1[kk] = Xmin[kk]+Xmax[kk]-Pg[kk];
		    if(Pg1[kk] < Xmin[kk]) Pg1[kk] = Xmin[kk];
		    if(Pg1[kk] > Xmax[kk]) Pg1[kk] = Xmax[kk];
		}
		else Pg1[kk] = Pg[kk];
	    }

	    Eval_Pg1 = getE(Pg1, ss, ll);
	    if(Eval_Pg1 < Eval_Pg){
		Eval_Pg = Eval_Pg1;
		for(int j=0; j<varNum; j++) Pg[j] = Pg1[j];
	    }
	}

	int e = uniform(engine)*Np;
	int q = uniform(engine)*Np;
	for(int k=0; k<varNum; k++){
	    Pg1[k] = Pg[k]+0.5*(X[e][k]-X[q][k]);
	    if(Pg1[k] < Xmin[k]) Pg1[k] = Xmin[k];
	    if(Pg1[k] > Xmax[k]) Pg1[k] = Xmax[k];
	}

	Eval_Pg1 = getE(Pg1, ss, ll);
	if(Eval_Pg1 < Eval_Pg){
	    Eval_Pg = Eval_Pg1;
	    for(int k=0; k<varNum; k++) Pg[k] = Pg1[k];
	}

	w -= del_w;
	a -= del_a;
	s -= del_s;

    }

    for(int i=0; i<4; i++)
	XX[i] = Pg[i];

    return Eval_Pg;
}
