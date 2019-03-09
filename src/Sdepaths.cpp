/*
 * Sdepaths.cpp
 *
 *  Created on: Dec 30, 2018
 *      Author: oliv
 */

#include "Sdepaths.h"
#include <random>
using namespace std;
using namespace arma;

namespace vSpace {


Sdepaths::Sdepaths(const realSpace& T,const double & S0,
		double (&a)(const double& t,const arma::vec& Path, const int& end),
		double (&b)(const double& t, const arma::vec& Path, const int& end),
		const int& N):T(T),a(a),b(b),N(N),proba(),logprice()
{
	if(N > 10000){
		throw "maximum 10 000 paths";
	}
	Paths = vector<arma::vec>(N,vec(T.getNx()+1));
//	t= 0
	vector<mt19937> generators;
	generators.reserve(N);

	string temp{};
	ifstream ip("10kseeds.txt");

	for(int i = 0 ; i != N ; ++i){
		Paths[i](0) = S0; /* Intialize */
		getline(ip,temp,'\n');
		generators.push_back(mt19937( stoi(temp)  )); /* seed the mersenne twisters */
	}
	ip.close();

	normal_distribution<double> Z;

	for(int j = 1 ; j != T.getNx()+1 ; ++j){
		for(int i = 0 ; i != N ; ++i){
			Paths[i](j) = Paths[i](j-1) + a(T(j),Paths[i],j-1)*T.getHx() + b(T(j),Paths[i],j-1)*sqrt(T.getHx())*Z(generators[i]);
//			Paths[i](j) = (generators[i])();

		}
	}


//	for(vector<arma::vec>::iterator it = Paths.begin(); it != Paths.end(); ++it){
//		cout << *it << '\n';
//	}



}

Sdepaths::~Sdepaths() {
	// TODO Auto-generated destructor stub
}

const arma::vec& Sdepaths::getPath(const int& i) const {
	return Paths[i];
}

double Sdepaths::E(const double& t,
		double (&f)(const double&)) const {
	double res{0};
	int j = round( (t-T.getXi()) /T.getHx() );

	for(int i = 0 ; i != N ; ++i){
		res += f( Paths[i](j) );
	}
	return res/N;
}

double Sdepaths::V(const double & t) const {
	return this->E(t,square)- this->E(t);
}

vfun Sdepaths::pdf(const double& t, int n) {
//	Find the max
	double maxi{0};

	int j = round( (t-T.getXi()) /T.getHx() );

	double mini{(Paths[0])(j)};
	double temp{0};
	for(int i = 0; i != N ; ++i){
		temp = (Paths[i])(j) ;
		maxi = temp < maxi ? maxi:temp;
		mini = temp > mini ? mini:temp;
	}
	double semirange = max(abs(log(mini)),abs(log(maxi)));
	if(semirange == 0 ){
		semirange = 0.1;
	}

	proba.reshape(n+1,1);


	double hx = 2*semirange/n;


	int c{0};
	for(int i = 0 ; i != N ; ++i){
		c = round( ( log( (Paths[i])(j) ) + semirange) / hx ) ;
		proba(c)++;
	}
	proba = proba/10000;
	return vfun(realSpace(-semirange,semirange,n),proba);
}

const arma::mat Sdepaths::getPathsMat() const {
	int p = T.getNx()+1;
	mat Paths_mat(N,p);
	for(int i = 0 ; i != N ; ++i){
		for(int j = 0 ; j != p ; ++j){
			Paths_mat(i,j) = Paths[i](j);
		}
	}
	return Paths_mat;
}

} /* namespace vSpace */
