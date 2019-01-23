/*
 * Mknpaths.cpp
 *
 *  Created on: Jan 20, 2019
 *      Author: oliv
 */

#include "Mknpaths.h"
#include <random>
#include <algorithm>
using namespace std;
using namespace arma;

namespace vSpace {

Mknpaths::Mknpaths(const realSpace & T,const double & S0,
		double (&a) (const double & t, const arma::vec & Path,const int & end),
		double (&b) (const double & t, const arma::vec & Path,const int & end),
		double (&sDup) (const double & t, const double & St),
		const int & N):T(T),a(a),b(b),sDup(sDup),N(N),proba(),logprice(){

		if(N > 10000){
			throw "maximum 10 000 paths";
		}
		const int M(T.getNx()+2);
		Paths = vector<arma::vec>(N,vec(M)); /* Last column stores the most recently
		calculated value for each row in order to enable a sort on spot val*/
		minIndexes = vec(N,fill::zeros);
	//	t= 0
		vector<mt19937> generators;
		generators.reserve(N);

		string temp{};
		ifstream ip("10kseeds.txt");


		for(int i = 0 ; i != N ; ++i){
			Paths[i](0) = S0; /* Intialize */
			Paths[i](M-1) = S0;
			getline(ip,temp,'\n');
			generators.push_back(mt19937( stoi(temp)  )); /* seed the mersenne twisters */
		}
		ip.close();

		normal_distribution<double> Z;

//      j = 1
		{
			int j{1};
				cout << "before sorting\n";

				for(int i = 0  ; i != N ; ++i ){
					cout << Paths[i]<<'\n';
				}



				cout << "after sorting\n";
				for(int i = 0  ; i != N ; ++i ){
					cout << Paths[i]<<'\n';
				}


				for(int i = 0 ; i != N ; ++i){
//					Make sure it is correct to remove both b and l in stochastic term and replace it by sDup for this stage
					Paths[i](j) = Paths[i](j-1) + a(T(j),Paths[i],j-1)*T.getHx() + sDup(T(j),S0)*sqrt(T.getHx())*Z(generators[i]);
					Paths[i](M-1) = Paths[i](j);
				}

		}


		for(int j = 2 ; j != T.getNx()+1 ; ++j){

			cout << "before sorting\n";

			double h = 1.5*S0*sqrt(max(T(j),0.25))*pow(N,-1./5);


			for(int i = 0  ; i != N ; ++i ){
				cout << Paths[i]<<'\n';
			}

			sort(Paths.begin(),Paths.end(),[M](arma::vec a, arma::vec b){if(a(M-1) < b(M-1)){return 1;}return 0;});
			/* Is mixing the seeds an issue ? */

			cout << "after sorting\n";
			for(int i = 0  ; i != N ; ++i ){
				cout << Paths[i]<<'\n';
			}

			for(int i = 0 ; i != N ; ++i){
				Paths[i](j) = Paths[i](j-1) + a(T(j),Paths[i],j-1)*T.getHx() + Paths[i](j-1)*b(T(j),Paths[i],j-1)*leverage(i,j,h)*sqrt(T.getHx())*Z(generators[i]);
				Paths[i](M-1) = Paths[i](j);
//				cout << 1./b(T(j),Paths[i],j)<<'\n';
			}


		}


}


Mknpaths::~Mknpaths() {
	// TODO Auto-generated destructor stub
}

double Mknpaths::leverage(const int& i, const int& j,const double & h) {
	double eps{0.01};
	double num{0};
	double den{0};
	int imin = minIndexes(i);
	double temp{0};
	double spot{Paths[i](j-1)};
	int k{imin};


	temp = kernel(Paths[k](j-1)-spot,h);

	while ( k>= 0 and temp > eps ){ /* We look left of imin */

//		cout << temp;

		num += temp;
		den += temp * pow(b(T(j),Paths[k],j-1),2);
		--k;
		temp = (k>=0)?kernel(Paths[k](j-1)-spot,h):0;
	}



//	The 0 might be missing if think not a big deal
	minIndexes(i) = max(k,0);
	k = imin + 1;


	temp = (k<N)?(Paths[k](j-1)-spot):0; /* beware temp is now kernless */


	while(k < N and temp < 0){ /* we take between imin and spot */


		temp = kernel(temp,h);
//		cout << temp;

		num += temp;
		den += temp * pow(b(T(j),Paths[k],j-1),2);
		++k;
		temp = (k<N)?(Paths[k](j-1)-spot):0;
	}


	temp = (k<N)?kernel(Paths[k](j-1)-spot,h):0;
	while(k<N and temp > eps){ /* we look right of spot */

//		cout << temp;


		num += temp;
		den += temp* pow(b(T(j),Paths[k],j-1),2);
		++k;
		temp = (k<N)?(kernel(Paths[k](j-1)-spot,h)):0;
	}
//	cout << sqrt(num/den)<<' ';
	return sDup(T(j), spot)*sqrt(num/den);
}

} /* namespace vSpace */
