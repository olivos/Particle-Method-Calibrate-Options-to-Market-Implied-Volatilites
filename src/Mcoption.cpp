/*
 * Mcoption.cpp
 *
 *  Created on: Dec 31, 2018
 *      Author: oliv
 */

#include "Mcoption.h"

namespace vSpace {


Mcoption::Mcoption(const fonction& S, const double& K, const double& r,const Sdepaths& P, double (&payoff) (const double & x, const double & strike)):Option(S,K,r),P(P),payoff(payoff) {
	T = final(0);
}

Mcoption::~Mcoption() {
	// TODO Auto-generated destructor stub
}

double Mcoption::operator ()(const double& t) const {
	double res{0};
	int N = P.getN();
	int j = round( (T) /deltas(0) );
	for(int i = 0 ; i != N ; ++i){
		res += exp(-r*(T-t))* payoff( P.getPath(i)(j),K );
	}
	return res/N;
}

double Mcoption::vol(const double& t) const {
	double c = this->operator ()(t);
	return func::findCallIV(c,1,K,r,T,0,1,100,0.00001);

}

double Mcoption::delta(const double& t) const {
}

} /* namespace vSpace */
