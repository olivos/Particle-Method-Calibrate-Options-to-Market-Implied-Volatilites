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
double a(const double& t,const arma::vec& Path, const int& end){
	return mu;
}
//sigma term
double b2(const double& t,const arma::vec& Path, const int& end){
	return Path(end)*0.15;
}

double b(const double& t,const arma::vec& Path, const int& end){
	int i = max(0,end - 20);
	if(Path(end) > mean(Path.subvec(i,end))  ){
		return Path(end)*0.10; /* Low vol */
	}
	else{
		return Path(end)*0.20; /* high vol */
	}
}
double price(const double & t){
	return 100.;
}

double pay(const double& x, const double & K){
	return (x-K)>0 ? (x-K):0;
}

double f(const double & t){
	return 100;
}

int main() {


	const double maturity = 1; /* maturity in years : here 3 months */
	const double r = 0;
	const double t = 0.;
	realSpace T{0,maturity,250}; /* Time space */



	Sdepaths P{T,100,a,b2,10000};

	fun S{T,f};

	Mcoption C(S,100.,r,P,pay);

	cout << C(t) <<'\n';

	cout <<func::callprice(100., 100., 0., maturity,  t,0.15);


	vfun dist = P.pdf(maturity,20);

//// Defining the option
//	fun S{T,price};
//	const int nK = 100;
//	realSpace K{100,1000,nK-1};
//	const int nT = 30;
//	realSpace Times{0.,0.9,nT-1};
//
//	mat Z(nK,nT,fill::zeros);
//	for(int i = 0 ; i != nK ; ++i){
//		Mcoption C(S,K(i),r,P,pay);
//	for(int j = 0 ; j!= nT ; ++j  ){
////		cout <<C(maturity-Times(j))<< ' ';
//		Z(i,j) = C(maturity-Times(j));
//		}
//	}

//	dataframe d{Z};
//	d.write_csv("Zpricesf.csv");

//	cout << Times << '\n';
//	cout << K << '\n';
}
