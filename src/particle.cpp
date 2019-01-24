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

double a(const double& t,const arma::vec& Path, const int& end){
	return 0;
}
double b(const double& t,const arma::vec& Path, const int& end){
	int i = max(0,end - 20);
	if(Path(end) > mean(Path.subvec(i,end))  ){
		return 0.10;
	}
	else{
		return 0.20;
	}
}
//double b(const double& t,const arma::vec& Path, const int& end){
//	return 0.2;
//}

double sdup(const double& t, const double& St ){
	return 0.3;
}
double K = 1; /* strike */
double d{0.01};


double pay2( const double & x){
	return (x>K+d)?(x-K-d):0;
}
double pay1( const double & x){
	return (x>K)?(x-K):0;
}
double pay0( const double & x){
	return (x>K-d)?(x-K+d):0;
}
int main() {

	double maturity = 1; /* maturity in years : here 3 months */

	double r = 0; /* risk free interest rate */
	double q = 0; /* continuous dividend */


	realSpace T{0,maturity,125}; /* Time space */
//	Mknpaths P{T,1,a,b,sdup,1000};


//	double CT0 =  P.E(0.8-d,pay1);
//	double CT2 =  P.E(0.8+d,pay1);
//
//
//	double CK0 =  P.E(0.8,pay0);
//	double CK1 =  P.E(0.8,pay1);
//	double CK2 =  P.E(0.8,pay2);

	double CT0 =  func::callprice(1,K,r,0.8-d,0,0.3,q);
	double CT2 =  func::callprice(1,K,r,0.8+d,0,0.3,q);

	double CK0 =  func::callprice(1,K-d,r,0.8,0,0.3,q);
	double CK1 =  func::callprice(1,K,r,0.8,0,0.3,q);
	double CK2 =  func::callprice(1,K+d,r,0.8,0,0.3,q);

//	0.0863425


	double dT = (CT2 - CT0)/(2*d);
	double ddK = (CK2 - 2*CK1 + CK0)/pow(d,2);

	cout << sqrt( 2*dT/(pow(K,2)*ddK) );

}
