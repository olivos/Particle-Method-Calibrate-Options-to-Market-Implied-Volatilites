/*
 * Sdepaths.h
 *
 *  Created on: Dec 30, 2018
 *      Author: oliv
 *
 *      N Paths for a SDE : dSt = a(t,X_t) dt + b(t,X_t) dBt where X_t = (S_u)_{u<=t}
 *
 */

#ifndef SDEPATHS_H_
#define SDEPATHS_H_

#include "topov2.h"


namespace vSpace {

class Sdepaths {
public:
	Sdepaths(const realSpace & T,const double & S0,
			std::function<double  (const double & t, const arma::vec & Path,const int & end)> &a,
			std::function<double  (const double & t, const arma::vec & Path,const int & end)> &b,
			const int & N);
	const arma::vec & getPath (const int & i) const;
	const arma::mat getPathsMat () const;
	double E(const double & t,std::function< double  (const double &)> &f ) const;
	double V(const double & t) const;
	vfun pdf(const double & t, int n = 100) ;
	virtual ~Sdepaths();

	void setPaths(const std::vector<arma::vec>& paths) {
		Paths = paths;
	}

	const int getN() const {
		return N;
	}

private:
	const realSpace& T;
	std::function<double (const double & t, const arma::vec & Path,const int & end)> &a;
	std::function<double (const double & t, const arma::vec & Path,const int & end)> &b;
	const int N;
	std::vector<arma::vec> Paths;
	/* To spare computation time the two following are not implemented unless pdf is called !! */
	arma::vec proba;
	realSpace logprice;

//	static std::function< double (const double & x)> Ide= [](const double & x) {return x;};
	static double square(const double & x) {
		return pow(x,2);
	}


};

} /* namespace vSpace */

#endif /* SDEPATHS_H_ */
