/*
 * Mknpaths.h
 *
 *  Created on: Jan 20, 2019
 *      Author: oliv
 *
 *      N Paths for a SDE : dSt = a(t,X_t) dt + b(t,X_t,P) St sigma_loc(t,St) dBt where X_t = (S_u)_{u<=t}
 *
 */

#ifndef MKNPATHS_H_
#define MKNPATHS_H_

#include "topov2.h"


namespace vSpace {

class Mknpaths {
public:
	Mknpaths(const realSpace & T,const double & S0,
			double (&a) (const double & t, const arma::vec & Path,const int & end),
			double (&b) (const double & t, const arma::vec & Path,const int & end),
			double (&sDup) (const double & t, const double & St),
			const int & N);
	virtual ~Mknpaths();

	const arma::vec & getPath (const int & i) const{
		return Paths[i];
	}


private:
	const realSpace& T;
	double (&a) (const double & t, const arma::vec & Path,const int & end);
	double (&b) (const double & t, const arma::vec & Path,const int & end);
	double (&sDup) (const double & t, const double & St);
	const int N;
	std::vector<arma::vec> Paths;
	arma::vec minIndexes;
	/* To spare computation time the two following are not implemented unless pdf is called !! */
	arma::vec proba;
	realSpace logprice;
	double leverage(const int& i, const int & j,const double & h);
	double kernel( const double & x,const double& h){
		if( abs(x) < 1  ){
			return 15./16/h *pow(1-pow(x/h,2),2) ;}
		return 0;
	}




};

} /* namespace vSpace */

#endif /* MKNPATHS_H_ */
