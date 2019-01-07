/*
 * Mcoption.h
 *
 *  Created on: Dec 31, 2018
 *      Author: oliv
 */

#ifndef MCOPTION_H_
#define MCOPTION_H_

#include <Option.h>
#include "Sdepaths.h"


namespace vSpace {

class Mcoption: public Option {
public:
	Mcoption(const fonction& S, const double& K, const double& r,const Sdepaths& P, double (&payoff) (const double & x, const double & strike) );
	double operator() (const double& t) const;
	double vol(const double& t)const; /* returns the implied vol */
	double delta(const double& t) const;
	virtual ~Mcoption();

private:
	double T;
	const Sdepaths & P;
	double (&payoff) (const double & x, const double & strike) ;
};

} /* namespace vSpace */

#endif /* MCOPTION_H_ */
