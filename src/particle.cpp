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
	return 0.1;
}


int main() {

	double maturity = 1; /* maturity in years : here 3 months */
	double t = 0; /* current time in years */
	double sigma = 0.2; /* annualized volatility */
	double r = 0; /* risk free interest rate */
	double q = 0; /* continuous dividend */
	double K = 1; /* strike */


	realSpace T{0,maturity,125}; /* Time space */
	Mknpaths P{T,1,a,b,sdup,100};

	cout << '\n'<<P.getPath(99);
}
