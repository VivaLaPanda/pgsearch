//lognormal dist tester
#include <iostream>
#include <random>
#include "dist_fit.hpp"

int main()
{
  DistFit distTest;
  const int nrolls= 100;  // number of experiments

  //std::default_random_engine generator;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::lognormal_distribution<double> distribution(0.0,1.0);

  for (int i=0; i<nrolls; ++i) {
    double number = distribution(gen);
    distTest.add_point(number);
  }

  std::cout << "lognormal_distribution (0.0,1.0):" << std::endl;
  std::cout << "MEAN: " << distTest.get_mean() << std::endl;
  std::cout << "VARIANCE: " << distTest.get_variance() << std::endl;



  return 0;
}
