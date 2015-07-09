//lognormal dist tester
#include <iostream>
#include <random>
#include "dist_fit.hpp"
#include <vector>

int main()
{
  DistFit distTest;
  const int nrolls= 100;  // number of experiments
  vector< double > edges;

  //std::default_random_engine generator;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::lognormal_distribution<double> distribution(0.0,1.0);

  for (int i=0; i<nrolls; ++i) {
    double number = distribution(gen);
    edges.push_back(number);
  }

  distTest.add_vector(edges);
  std::cout << "lognormal_distribution (0.0,1.0):" << std::endl;
  std::cout << "MEAN: " << distTest.get_logMV()[0][0] << std::endl;
  std::cout << "VARIANCE: " << distTest.get_logMV()[0][1] << std::endl;



  return 0;
}
