//lognormal dist tester
#include <iostream>
#include <random>
#include "dist_fit.hpp"
#include <vector>

int main()
{
  DistFit distTest;
  const int nrolls= 100;  // number of experiments

  //std::default_random_engine generator;
  std::random_device rd;
  std::mt19937 gen(rd());
  double m, v, mu_log, sigma_log;
  m = 1.0;
  v =2.0;
  mu_log = log((pow(m,2))/sqrt(v+pow(m,2)));
  sigma_log = sqrt(log(v/(pow(m,2))+1));
  std::lognormal_distribution<double> distribution(mu_log,sigma_log);
  std::cout << "lognormal_distribution (" << m << "," << v << "):" << std::endl;
  int num_steps = 10000;
  vector< double > edges;
  vector < vector < double > > logMV, MV;
  for (int i = 0; i < num_steps; i++){
	  for (int i=0; i<nrolls; ++i) {
		double number = distribution(gen);
		edges.push_back(number);
		//cout << "NUMBER: " << number << endl;
	  }
	  distTest.add_vector(edges);
	  logMV =  distTest.get_logMV();
	  MV = distTest.get_normMV();
	  edges.clear();
  }
  for(int i = 0; i < logMV.size(); i++){
	  //std::cout << "NORMAL MEAN: " << MV[i][0] << std::endl;
	  //std::cout << "NORMAL VARIANCE: " << MV[i][1] << std::endl;
	  std::cout << "LOG MEAN: " << logMV[i][0] << std::endl;
	  std::cout << "LOG VARIANCE: " << logMV[i][1] << std::endl;
  }

  return 0;
}
