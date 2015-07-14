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
  double m, u;
  m = 1.0;
  u =2.0;
  std::lognormal_distribution<double> distribution(u,m);
  std::cout << "lognormal_distribution (" << m << "," << u << "):" << std::endl;
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
