// Simulating PDV model and printing the distribution
// Sdepaths simulates the paths of this stochastic diffrential equation :
// dSt = a(t,X,St)dt + b(t,X,St)dBt


double a(const double& t,const arma::vec& Path, const int& end){
	return 0;
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



int main() {

	double maturity = 1; /* maturity in years : here 3 months */

	realSpace T{0,maturity,125}; /* Time space */
	Sdepaths P{T,1,a,b,10000}; /* Montecarlo simulation */
	vfun dist = P.pdf(maturity,20);
	cout << dist.getX()<< '\n';
	cout << dist;


}


// Using Sde Paths to price an option

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
}
//double b(const double& t,const arma::vec& Path, const int& end){
//	return 0.2*Path(end);
//}
//double Ide(const double& x){
//	return x;
//}
//double payoff(const double& x){
//	double K{0.2};
//	return (x-K)>0 ? (x-K):0;
//}
double pay(const double& x, const double & K){
	return (x-K)>0 ? (x-K):0;
}

double price( const double& t){
return 1;
}
//
//int f(int a, int b){
//	return a + b;
//}


int main() {

	double maturity = 1; /* maturity in years : here 3 months */
	double t = 0; /* current time in years */
	double sigma = 0.2; /* annualized volatility */
	double r = 0; /* risk free interest rate */
	double q = 0; /* continuous dividend */
	double K = 1; /* strike */


	realSpace T{0,maturity,125}; /* Time space */
	Sdepaths P{T,1,a,b,10000};
	vfun dist = P.pdf(maturity,20);

// Defining the option
	fun S{T,price};
	Mcoption C(S,1,r,P,pay);

	cout << C(0); /* option price */
	cout << C.vol(0); /* option vol */


}