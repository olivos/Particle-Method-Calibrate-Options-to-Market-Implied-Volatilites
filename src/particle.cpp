//============================================================================
// Name        : test.cpp
// Author      : Olivier Ort-Snep
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================


#include <iostream>
#include <random>
#include <armadillo>
#include <thread>

#include "topov2.h"
#include "Sdepaths.h"
#include "Mcoption.h"
#include "Mknpaths.h"


using namespace std;
using namespace arma;
using namespace vSpace;

const double sigmaR = 0.3; /* Realized vol */
const double mu = 0.;

//Drift term
inline double a(const double& t,const arma::vec& Path, const int& end){
	return (0.004-0.0249)*Path[end];
}
//sigma term

inline double b2(const double& t,const arma::vec& Path, const int& end){
	return Path(end)*0.15;
}

//inline double b(const double& t,const arma::vec& Path, const int& end){
//	int i = max(0,end - 20);
//	if(Path(end) > mean(Path.subvec(i,end))  ){
//		return Path(end)*0.20; /* Low vol */
//	}
//	else{
//		return Path(end)*0.10; /* high vol */
//	}
//}

double b(const double& t,const arma::vec& Path, const int& end,const realSpace & T,const double & sigmaL, const double & sigmaH, const double & days){
	int delta = days/252*(T.getNx()+1)/T.getXf( );
	int i = max(0,end - delta);
//	cout <<mean(Path.subvec(i,end))<<endl;
	if(Path(end) >= mean(Path.subvec(i,end))  ){
//		cout<< "low vol "<<endl;

		return Path(end)*sigmaL; /* Low vol */
	}
	else{
//		cout<< "high vol "<<endl;

		return Path(end)*sigmaH; /* high vol */

	}
}

inline double price(const double & t){
	return 100.;
}

inline double pay(const double& x, const double & K){
	return (x-K)>0 ? (x-K):0;
}

double err (const double S0,const arma::mat DKP, realSpace & T,const arma::mat & v){
//	 v_0 = sigmaL , v_1 = sigmaH , v_2 = days for mvavg
	const int npaths = 1000;
	const int nprices = DKP.n_rows;
	double err = 0;
	double pred_price;
	double K;
	double maturity;
	std::function<double  (const double & t, const arma::vec & Path,const int & end)> aa = [](const double & t, const arma::vec & Path,const int & end){return a(t,Path,end);};
	std::function<double  (const double & t, const arma::vec & Path,const int & end)> bb = [v,T](const double & t, const arma::vec & Path,const int & end){return b(t, Path, end, T,v(0), v(1), v(2));};
	Sdepaths P{T,S0,aa,bb,npaths};
	for(int i = 0 ; i != nprices ; ++i){
		maturity = DKP(i,0);
		K = DKP(i,1);
		std::function<double(const double & x)> payoff = [K](const double & x){return (x-K)>0?(x-K):0;};
		pred_price = P.E(  maturity, payoff );
		err += pow(  pred_price - DKP(i,2)  ,2);
	}
//	cout << pred_price;
	return err/nprices;
}

void err_threaded (double * res,const double S0,const arma::mat * DKP,const realSpace T,const arma::mat v){
//	 v_0 = sigmaL , v_1 = sigmaH , v_2 = days for mvavg
	const int npaths = 1000;
	const int nprices = (*DKP).n_rows;
	double err = 0;
	double pred_price;
	double K;
	double maturity;
	std::function<double  (const double & t, const arma::vec & Path,const int & end)> aa = [](const double & t, const arma::vec & Path,const int & end){return a(t,Path,end);};
	std::function<double  (const double & t, const arma::vec & Path,const int & end)> bb = [v,T](const double & t, const arma::vec & Path,const int & end){return b(t, Path, end, T,v(0), v(1), v(2));};
	Sdepaths P{T,S0,aa,bb,npaths};
	for(int i = 0 ; i != nprices ; ++i){
		maturity = (*DKP)(i,0);
		K = (*DKP)(i,1);
		std::function<double(const double & x)> payoff = [K](const double & x){return (x-K)>0?(x-K):0;};
		pred_price = P.E(  maturity, payoff );
		err += pow(  pred_price - (*DKP)(i,2)  ,2);
	}
//	cout << pred_price;
	*res = err/nprices;
}

arma::mat grad(const double S0,const arma::mat DKP, realSpace & T,const arma::mat & v){
	vec one{1,0,0};
	vec two{0,1,0};
	vec three{0,0,1};
	vec res(3,fill::zeros);
	const double dv= 1e-5;
	double oneL,oneR,twoL,twoR;

	thread th1{err_threaded,&oneR,S0,&DKP,T,v+dv*one};
	thread th2{err_threaded,&oneL,S0,&DKP,T,v-dv*one};
	thread th3{err_threaded,&twoR,S0,&DKP,T,v+dv*two};
	thread th4{err_threaded,&twoL,S0,&DKP,T,v-dv*two};

	th1.join();
	th2.join();
	th3.join();
	th4.join();


	res(0) = (oneR - oneL)/(2*dv);
	res(1) = (twoR - twoL)/(2*dv);


	return res;
}

arma::mat search(const double S0, const arma::mat DKP, realSpace & T){
	const int grid_len = 100;
	const double start = 0.01;
	const double end = 1.;
	const double dv = (end-start)/(grid_len-1);
	vec v{start,start,20};
	mat ERR(grid_len,grid_len,fill::zeros);

	int i_min = 0;
	int j_min = 0;
	double min = err(S0,DKP,T,v);


	for(int i = 0 ; i != grid_len ; ++i){
		v(0) = start+dv*i;
		for(int j = i; j != grid_len ; ++j){
			v(1) = start+dv*j;
			ERR(i,j) = err(S0,DKP,T,v);
			if( min > ERR(i,j)){
				i_min = i;
				j_min = j;
				min = ERR(i,j);
			}
		}
	}
	cout << i_min <<' '<< j_min<<'\n';
	cout << start+dv*i_min << ' '<<start+dv*j_min<<'\n';
	return ERR;
}

arma::vec predict (const double S0,const arma::mat DKP, realSpace & T,const arma::mat & v){
//	 v_0 = sigmaL , v_1 = sigmaH , v_2 = days for mvavg
	const int npaths = 1000000;
	const int nprices = DKP.n_rows;
	double K;
	double maturity;
	std::function<double  (const double & t, const arma::vec & Path,const int & end)> aa = [](const double & t, const arma::vec & Path,const int & end){return a(t,Path,end);};
	std::function<double  (const double & t, const arma::vec & Path,const int & end)> bb = [v,T](const double & t, const arma::vec & Path,const int & end){return b(t, Path, end, T,v(0), v(1), v(2));};
	Sdepaths P{T,S0,aa,bb,npaths};

	vec Preds(nprices,fill::zeros);

	for(int i = 0 ; i != nprices ; ++i){
		maturity = DKP(i,0);
		K = DKP(i,1);
		std::function<double(const double & x)> payoff = [K](const double & x){return (x-K)>0?(x-K):0;};
		Preds(i) = P.E(  maturity, payoff );
	}
	return Preds;
}


arma::mat grad_descent(const arma::mat& v0, std::function<arma::mat(const arma::mat &)> & grad,const double & eth,const double & eps){
	mat v{v0};
	mat g = grad(v);
	const int maxIt{100};
	for(int i = 0 ; i != maxIt ; ++i ){
		if(norm(g) < eps){
			cout << "numit" << i <<' ' << norm(g) <<'\n';
			return v;
		}
		else{
			v = v-eth*g;
			v(0) = (v(0)>0)?v(0):0.01;
			v(1) = (v(1)>v(0))?v(1):(v(0)+0.001);
			g = grad(v);
		}
	}
	cout << "reached max it" << norm(g) << '\n';
	return v;

}


int main() {


	const double maturity = 2.290411; /* maturity in years : here 3 months */
	const double r = 0.;
	const double t = 0.;
	const double S0 = 1403.44;
	realSpace T{0,maturity,100}; /* Time space */

	std::function<double  (const double & t, const arma::vec & Path,const int & end)> aa = [](const double & t, const arma::vec & Path,const int & end){return a(t,Path,end);};
	std::function<double  (const double & t, const arma::vec & Path,const int & end)> bb = [T](const double & t, const arma::vec & Path,const int & end){return b(t,Path,end,T,0.2,0.2,20);};

//	Sdepaths P{T,S0,aa,bb,1000000};
	std::function<double(const double & x)> payoff = [S0](const double & x){return (x-S0)>0?(x-S0):0;};

//	cout <<  P.E(1.,payoff)<<' ';

	dataframe dkp{222,3,"dkp.csv"};
	mat DKP =  dkp.getData();
	vec v{0.05,0.34,20};

	auto start = std::chrono::high_resolution_clock::now();

	vec Preds = predict(S0,DKP,T,v);

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsedC = finish - start;


	dataframe preds(Preds);
	preds.write_csv("preds.csv");
	cout << elapsedC.count();


//	auto start = std::chrono::high_resolution_clock::now();
//
//
//	mat ERR = search(S0,DKP,T);
//
//	auto finish = std::chrono::high_resolution_clock::now();
//	std::chrono::duration<double> elapsedC = finish - start;
//
//	cout <<"time"<< elapsedC.count()<<"\n\nPred\n";
//	print(S0,DKP,T,v);
//	cout << "DKP\n"<<DKP.col(2) <<'\n';
//	dataframe err(ERR);
//	err.write_csv("err.csv");




}

//	cout << DKP;

//	std::function<arma::mat(const arma::mat &)> g = [&](const arma::mat & v){return grad(S0,DKP,T,v);};
//
//
//		auto start = std::chrono::high_resolution_clock::now();
//
//		vec vinf = grad_descent(v,g,1e-3,1e-6);
//
//
//		auto finish = std::chrono::high_resolution_clock::now();
//		std::chrono::duration<double> elapsedC = finish - start;
//
//		cout << elapsedC.count()<<"\n\nPred\n";
//		print(S0,DKP,T,v);
//		cout << "DKP\n"<<DKP.col(2) <<'\n';
//
////	cout << P.E(1.,payoff);
//
//		cout << "vinf\n"<<vinf<<'\n';
//		double d;
//		err_threaded(&d,S0,&DKP,T,vinf);
//		cout << "err\n"<< d;

