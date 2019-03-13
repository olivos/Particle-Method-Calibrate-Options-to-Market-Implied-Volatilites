//============================================================================
// Name        : test.cpp
// Author      : Olivier Ort-Snep
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================


#include <iostream>
#include <random>

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
	return mu;
}
//sigma term
inline double b2(const double& t,const arma::vec& Path, const int& end){
	return Path(end)*0.15;
}

inline double b(const double& t,const arma::vec& Path, const int& end){
	int i = max(0,end - 20);
	if(Path(end) > mean(Path.subvec(i,end))  ){
		return Path(end)*0.20; /* Low vol */
	}
	else{
		return Path(end)*0.10; /* high vol */
	}
}
inline double price(const double & t){
	return 100.;
}

inline double pay(const double& x, const double & K){
	return (x-K)>0 ? (x-K):0;
}



int main() {


	const double maturity = 1; /* maturity in years : here 3 months */
	const double r = 0;
	const double t = 0.;
	realSpace T{0,maturity,100}; /* Time space */

	Sdepaths P{T,100,a,b,10000};


// Defining the option
	fun S{T,price};
	Mcoption C{S,80.,r,P,pay};
//	cout << C(t) << ' ';
//	cout << func::findCallIV(C(t),100.,80., 0., maturity-t,0., 1.,100);


	const int nK = 100;
	realSpace K{70,120,nK-1};
	const int nT = 30;
	realSpace Times{0.,0.9,nT-1};

	auto start = std::chrono::high_resolution_clock::now();

	mat Z(nK,nT,fill::zeros);
	for(int i = 0 ; i != nK ; ++i){
		Mcoption C(S,K(i),r,P,pay);
	for(int j = 0 ; j!= nT ; ++j  ){
//		cout <<C(maturity-Times(j))<< ' ';
		Z(i,j) = func::findCallIV(C(Times(j)),100.,K(i), r, maturity-Times(j),0., 1.,100);
		}
	}

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	cout << elapsed.count()<<endl;

	dataframe d{Z};
	d.write_csv("IVc.csv");

	cout << Times << '\n';
	cout << K << '\n';
}
