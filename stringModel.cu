// -*- C++ -*-
/******************************************************************************
 *            stringModel.cu
 *****************************************************************************/
#include "globalConstants.hh"
#include "stringModel.hh"
#include "bisection.hh"

/** Constructor    *****************************************************************
 *  GG (W/m^2), TT (degC)
 *  出力： m_gc, m_G, m_T, m_ssInit, m_Ib0, m_etabVt, m_IscSub, m_Vmin, m_Vmax
 *  関数呼び出し： Substring::getIb0, Substring::getEtabVt, Substring::getIsc
 *                 getV     ********************************************************/
StringModel::StringModel(GlobalConstants &gcc, double GG, double TT): m_gc{gcc}, m_G{GG}, m_T{TT},
								   m_ssInit{gcc, GG, TT}{
    m_Ib0 = m_ssInit.getIb0();
    m_etabVt = m_ssInit.getEtabVt();

    m_IscSub = m_ssInit.getI(0.0);
    m_mIb099 = -0.99999*m_Ib0;
    m_IscSub2 = 2.0*m_IscSub;
    
    try{
	m_Vmin = getV(2.0*m_IscSub);
	m_Vmax = getV(m_mIb099);
    }catch(std::exception e){
	std::cout << e.what() << "-->StringModel::StringModel\n";
	exit(0);//////////////t
    }

    // // printf("vss_size=%d, L=%d\n", m_vss.size(), m_gc.L);////////////////////////////////
    // std::cout << "StringModel Vmin=" << m_Vmin << " Vmax=" << m_Vmax << std::endl;////////////
}

/** m_aVmin[], m_aVmax[], m_aIb0[], m_aIscSub[], m_aEtabVt[]をセットする   **********
 *  入力： 引き数だけ
 *  出力： m_ssInit, m_vss, m_aIb0, m_aEtabVt, m_aIscSub, m_aVmin, m_aVmax
 *  関数呼び出し： Substring::setVars, Substring::getIb0(k), Substring::getEtabVt(k)
 *                 Substring::getIsc, getVIjunbi
 *  戻り値： Substring::setVarsの戻り値     *****************************************/
double StringModel::setVars(double GG[], double TT[]){
    double Eret = m_ssInit.setVars(GG, TT);
    if(Eret==INFINITY) return Eret;
    
    for(int i=0; i<m_vss.size(); i++){
	double E = m_vss[i].setVars(GG, TT);
	if(E==INFINITY) return E;
	else if(Eret < E) Eret= E;
    }
    
    for(int k=0; k<GlobalConstants::K; k++){
	m_aIb0[k] = m_ssInit.getIb0(k);      // 逆流防止ダイオードのもの
	m_aEtabVt[k] = m_ssInit.getEtabVt(k);

	m_aIscSub[k] = m_ssInit.getIsc(k);
	for(int i=0; i<m_vss.size(); i++){
	    double IscSub_i = m_vss[i].getIsc(k);
	    if(IscSub_i > m_IscSub) m_aIscSub[k] = IscSub_i;
	}

	getVjunbi(k);
	
	try{
	    m_aVmin[k] = getV(2.0*m_aIscSub[k]);
	    m_aVmax[k] = getV(-0.99*m_aIb0[k]);
	}catch(std::exception e){
	    std::cout << e.what() << "-->StringModel::setVars\n";
	    exit(0);//////////////t
	}
    }

	// printf("StringModel::setVars_k$$ k=%d, Vmin=%g, Vmax=%g\n", k, m_aVmin[k], m_aVmax[k]);////////////t
    return Eret;
}

/** m_vss[l]のsetVarsを実行する   **********
 *  入力： 引き数だけ
 *  出力： m_vss, m_aIscSub, m_aVmin, m_aVmax
 *  関数呼び出し： Substring::setVars, Substring::getIsc, getVIjunbi, getV
 *  戻り値： Substring::setVarsの戻り値     *****************************************/
double StringModel::setVarsVss(double GG[], double TT[], int l){
    double E =  m_vss[l].setVars(GG, TT);
    if(E==INFINITY) return E;

    for(int k=0; k<GlobalConstants::K; k++){
	double IscSub_l = m_vss[l].getIsc(k);
	// printf("IscSub_l=%g\n", IscSub_l);/////////////////t
	if(IscSub_l > m_aIscSub[k]){
	    m_aIscSub[k] = IscSub_l;
	    getVjunbi(k);
	    
	    try{
		m_aVmin[k] = getV(2.0*m_aIscSub[k]);
		// m_aVmax[k] = getV(-0.99*m_aIb0[k]);
	    }catch(std::exception e){
		std::cout << e.what() << "-->StringModel::setVarsVss" << " k="<< k
			  << " l=" << l << "\n";
		exit(0);//////////////t
	    }
	}
    }

    return E;
}

/** m_vss[l]のsetVarsを実行する, StringGroup::getE用   **********
 *  入力： 引き数だけ
 *  出力： m_vss, m_aIscSub, m_aVmin, m_aVmax
 *  関数呼び出し： Substring::setVars, Substring::getIsc, getVIjunbi, getV
 *  戻り値： Substring::setVarsの戻り値     *****************************************/
double StringModel::setVarsVssGetE(int l){
    double E =  m_vss[l].setVarsGetE();
    // return m_vss[l].setVarsGetE();
    if(E==INFINITY) return E;

    for(int k=0; k<GlobalConstants::K; k++){
	double IscSub_l = m_vss[l].getIsc(k);
	// printf("IscSub_l=%g\n", IscSub_l);/////////////////t
	if(IscSub_l > m_aIscSub[k]){
	    m_aIscSub[k] = IscSub_l;
	    getVjunbi(k);
	    
	    try{
		m_aVmin[k] = getV(2.0*m_aIscSub[k]);
		// m_aVmax[k] = getV(-0.99*m_aIb0[k]);
	    }catch(std::exception e){
		std::cout << e.what() << "-->StringModel::setVarsVss" << " k="<< k
			  << " l=" << l << "\n";
		exit(0);//////////////t
	    }
	}
    }

    return E;
}

/** 時点kに対するgetVの準備をする    ************************************************
 *  入力： m_aIb0[k], m_aEtabVt[k]
 *  出力： m_Ib0, m_etabVt
 *  関数呼び出し： m_ssInit.getVIjunbi, m_vss[i].getVIjunbi   ***********************/
void StringModel::getVjunbi(int k){
    m_Ib0 = m_aIb0[k];
    m_etabVt = m_aEtabVt[k];

    m_ssInit.getVIjunbi(k);
    for(int i=0; i<m_vss.size(); i++){
	m_vss[i].getVIjunbi(k);
    }
}

/** 時点kに対するgetVの準備をする    ************************************************
 *  入力： m_aIb0[k], m_aEtabVt[k]
 *  出力： m_Ib0, m_etabVt
 *  関数呼び出し： m_ssInit.getVIjunbi, m_vss[i].getVIjunbi   ***********************/
void StringModel::getVjunbiVss(int k, int l){
    m_Ib0 = m_aIb0[k];
    m_etabVt = m_aEtabVt[k];

    m_vss[l].getVIjunbi(k);
}

/** Calculate substring voltage V (V) for current I (A)   ***************************
 *  入力：  m_Ib0, m_etabVt
 *  関数呼び出し： m_vss.size, m_ssInit.getV, m_vss[i].getV     *********************/
double StringModel::getV(double I){
    if(I <= -m_Ib0) return std::numeric_limits<double>::quiet_NaN();

    double V0 = -m_etabVt*std::log(I/m_Ib0 +1.0);
    double V;
    try{
	V = V0 + (m_gc.L-m_vssNumSum)*m_ssInit.getV(I);
	for(int i=0; i<m_vss.size(); i++) V += m_vss[i].getV(I)*m_vssNum[i];
    }catch(std::exception e){
	std::string str {e.what()};
	throw std::runtime_error(str+"-->StringModel::getV");
    }
    
    return V;
}

/** 時点kに対するgetIの準備をする     ***********************************************
 *  入力： m_aIb0[k], m_aEtabVt[k]
 *  出力： m_Ib0, m_etabVt, m_Vmin, m_Vmax, m_IscSub
 *  関数呼び出し： m_ssInit.getVIjunbi, m_vss[i].getVIjunbi     *********************/
void StringModel::getIjunbi(int k){
    m_Ib0 = m_aIb0[k];
    m_etabVt = m_aEtabVt[k];

    m_ssInit.getVIjunbi(k);
    for(int i=0; i<m_vss.size(); i++){
	m_vss[i].getVIjunbi(k);
    }

    m_Vmin = m_aVmin[k];
    m_Vmax = m_aVmax[k];
    m_IscSub = m_aIscSub[k];
}

/** substring lのパラメータが変わったときの、時点kに対するgetIの準備をする    ********
 *  ssInitとlより前のsubstringのgetIjunbiは終わっているとする
 *  入力： m_aIb0[k], m_aEtabVt[k]
 *  出力： m_Ib0, m_etabVt, m_Vmin, m_Vmax, m_IscSub
 *  関数呼び出し： m_vss[l].getVIjunbi     *********************/
void StringModel::getIjunbi(int k, int l){
    m_Ib0 = m_aIb0[k];
    m_etabVt = m_aEtabVt[k];

    m_vss[l].getVIjunbi(k);
    
    m_Vmin = m_aVmin[k];
    m_Vmax = m_aVmax[k];
    m_IscSub = m_aIscSub[k];
}

/** Calculate substring current I (A) for voltage V (V)    **************************
 *  入力： m_Vmin, m_Vmax, m_Ib0, m_IscSub
 *  出力： m_Vfunc
 *  関数呼び出し： Bisection::search     ********************************************/
double StringModel::getI(double V){
    if(V < m_Vmin) return std::numeric_limits<double>::infinity();
    else if(m_Vmax < V) return -m_Ib0;
    else{
	Bisection<StringModel> bis {this};
	m_Vfunc = V;
	// return bis.search(-0.99*m_Ib0, 2.0*m_IscSub);
	return bis.search(m_mIb099, m_IscSub2);
    }
}

/** Calculate substring current Ik (A) for voltage Vk (V), (k=0,1,..K-1)    *********
 *  l はXによってVocSTCなどのパラメータが変わっているサブストリング番号
 *  getIjunbiの事前呼び出しは不要
 *  上のgetIをK回繰り返すより高速になると思ったが、かえって遅くなった --> 不採用
 *  入力： m_Vmin, m_Vmax, m_Ib0, m_IscSub
 *  出力： m_Vfunc
 *  関数呼び出し： Bisection::search     ********************************************/
void StringModel::getI(int l, double Vk[GlobalConstants::K], double Ik[GlobalConstants::K]){
    for(int k=0; k<GlobalConstants::K; k++){
	m_Ib0 = m_aIb0[k];
	m_etabVt = m_aEtabVt[k];

	m_vss[l].getVIjunbi(k);

	if(Vk[k] < m_aVmin[k]) Ik[k] = std::numeric_limits<double>::infinity();
	else if(m_aVmax[k] < Vk[k]) Ik[k] = -m_Ib0;
	else{
	    Bisection<StringModel> bis {this};
	    m_Vfunc = Vk[k];
	    Ik[k] = bis.search(-0.99*m_Ib0, 2.0*m_aIscSub[k]);
	}
    }
}

/** getIで使われる２分法用のメソッド   **********************************************
 *  入力： m_Vfunc
 *  関数呼び出し： getV     *********************************************************/
double StringModel::func(double I){
    return getV(I)-m_Vfunc;
}

/****************   Add a substring     ***********************************************/
void StringModel::addSS(){
    // Substring ss {m_ssInit};   // default copy constructor
    Substring ss {m_gc, 1000.0, 25.0};
    m_vss.push_back(ss);
}

/****************   Add substrings     ***********************************************/
void StringModel::addSS(Substring &ss, int l){
    m_vss.push_back(ss);
    m_vssNum.push_back(l);
    m_vssNumSum += l;

    double Isc_ss {ss.getI(0.0)};
    if(Isc_ss > m_IscSub) m_IscSub = Isc_ss;
			      
    m_Vmin = getV(2.0*m_IscSub);
    m_Vmax = getV(m_mIb099);
}

/******************************************************************************/
void StringModel::setGT(double G, double T){
    m_G = G;
    m_T = T;
    m_ssInit.setGT(G, T);
    for(int i=0; i<m_vss.size(); i++) m_vss[i].setGT(G, T);

    m_Ib0 = m_ssInit.getIb0();
    m_etabVt = m_ssInit.getEtabVt();

    m_IscSub = m_ssInit.getI(0.0);
    m_mIb099 = -0.99999*m_Ib0;
    m_IscSub2 = 2.0*m_IscSub;

    for(int i=0; i<m_vss.size(); i++){
	double Isc_ss {m_vss[i].getI(0.0)};
	if(Isc_ss > m_IscSub) m_IscSub = Isc_ss;
    }
    
    m_Vmin = getV(2.0*m_IscSub);
    m_Vmax = getV(m_mIb099);
}

/******* Suppose refreshment after set5params of adaptiveModel ************/
void StringModel::refreshMaxMin(){
    m_ssInit.refreshMaxMin();
    for(int i=0; i<m_vss.size(); i++) m_vss[i].refreshMaxMin();
    
    m_IscSub = m_ssInit.getI(0.0);
    m_IscSub2 = 2.0*m_IscSub;

    for(int i=0; i<m_vss.size(); i++){
	double Isc_ss {m_vss[i].getI(0.0)};
	if(Isc_ss > m_IscSub) m_IscSub = Isc_ss;
    }
    
    m_Vmin = getV(2.0*m_IscSub);
    m_Vmax = getV(m_mIb099);
}

/** サブストリングlのSTCパラメータをXでセットする    *********************************
 *  関数呼び出し： Substring::setAM_STCparams     ***********************************/
void StringModel::setSS_STCparams(double X[4], int l){
    m_vss[l].setAM_STCparams(X);
}

/** 時点kに対するストリングの電流m_aIkをセットする    *******************************
 *  出力： m_aIk[]     *************************************************************/
void StringModel::set_aIk(double Ik[GlobalConstants::K]){
    for(int k=0; k<GlobalConstants::K;k++){
	m_aIk[k] = Ik[k];
    }
}

/** 電圧Vに対する現在のパラメータでのストリング電流 Ikをm_aIkにセットする    *******/
void StringModel::set_aIk_for(double V[GlobalConstants::K]){
    for(int k=0; k<GlobalConstants::K;k++){
	getIjunbi(k);
	m_aIk[k] = getI(V[k]);
    }
}

/**   ****************************************************************************/
FiveParams StringModel::get5params(int l) const{
    if(l >= 0) return m_vss[l].get5params();
    else return m_ssInit.get5params();
}

/**   ****************************************************************************/
FiveParamsPlus StringModel::get5paramsPlus(int l) const{
    if(l >= 0) return m_vss[l].get5paramsPlus();
    else return m_ssInit.get5paramsPlus();
}

/**     **************************************************************************/
void StringModel::print(int l) const{
    if(l>=0) m_vss[l].print();
    else m_ssInit.print();
    for(int k=0; k<GlobalConstants::K; k++){
	printf("StringModel: k=%d, aVmin=%g, aVmax=%g, aIb0=%g, aIscSub=%g, aEtabVt=%g, aIk=%g\n",
	       k, m_aVmin[k], m_aVmax[k], m_aIb0[k], m_aIscSub[k], m_aEtabVt[k], m_aIk[k]);
    }
}

