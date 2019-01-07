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


using namespace std;
using namespace arma;


double fu(const double& x){
	return pow(x,2);
}
using namespace vSpace;

double a(const double& t,const arma::vec& Path, const int& end){
	return 0;
}
double b(const double& t,const arma::vec& Path, const int& end){
	int i = max(0,end - 20);
	if(Path(end) > mean(Path.subvec(i,end))  ){
		return Path(end)*0.10;
	}
	else{
		return Path(end)*0.20;
	}
//	return  Path(end)*0.2;
}
//double b(const double& t,const arma::vec& Path, const int& end){
//	return 0.2*Path(end);
//}
double Ide(const double& x){
	return x;
}
double payoff(const double& x){
	double K{0.2};
	return (x-K)>0 ? (x-K):0;
}
double pay(const double& x, const double & K){
	return (x-K)>0 ? (x-K):0;
}

double price( const double& t){
return 1;
}

int f(int a, int b){
	return a + b;
}


int main() {

	double maturity = 0.5; /* maturity in years : here 3 months */
	double t = 0; /* current time in years */
	double sigma = 0.2; /* annualized volatility */
	double r = 0; /* risk free interest rate */
	double q = 0; /* continuous dividend */
	double K = 1; /* strike */


	realSpace T{0,maturity,125}; /* Time space */
	Sdepaths P{T,1,a,b,10000};
//	cout << P.pdf(0.2,20);

	fun S{T,price};

	realSpace Strikes{0.87,1.07,25};
	vec IV{26,fill::zeros};

	for(int i = 0 ; i != 26 ; ++i ){
		Mcoption C(S,Strikes(i),r,P,pay);
		IV(i) = C.vol(0);
//		cout << s<< ' '<<C(0) << ' '<< func::callprice(1,Strike,r,maturity,0,s,q)<<'\n';
	}

	vfun smile = vfun(Strikes, IV);
	outputC::plot(Strikes,smile);



}
