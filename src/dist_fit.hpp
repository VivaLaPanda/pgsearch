#include <vector>
#include <iostream> 
#include <numeric>
#include <math.h>
//#include <boost/accumulators/numeric/functional/vector.hpp>
//#include <boost/accumulators/statistics/variance.hpp>
//#include <boost/fusion/algorithm/iteration/for_each.hpp>
//#include <boost/fusion/include/for_each.hpp>
//vector< vector< double > > myvec ;
using namespace std;

class DistFit{
		vector< double > points;
		vector< double > log_points;
		double mean, variance;
	public:
		DistFit();
		void add_point(double new_point);
		double get_mean() const {return mean; };
		double get_variance() const {return variance; };
};

DistFit::DistFit(){
	mean = 0;
	variance = 0;
}


void DistFit::add_point(double new_point){
	points.push_back(new_point);
	log_points.push_back(log(new_point));

	double sum = accumulate(log_points.begin(), log_points.end(), 0.0);
	double mu = sum / log_points.size();
	double sumSquares = inner_product(log_points.begin(), log_points.end(), log_points.begin(), 0.0);
	double sigma = sumSquares/log_points.size() - mu*mu;

	mean= exp(mu + pow(sigma,2)/2);
	variance = exp(2*mu + pow(sigma,2))*(exp(pow(sigma,2)-1));

	cout << "Mean: " << mean << endl;
	cout << "Variance: " << variance << endl;
}


/*
void gd_slope(float &slope_mu, float &slope_sigma) {
	int nmax = residual.size();
	float tmp_sigma = 0.f, tmp_mu = 0.f;
	for(int i = 0; i < nmax; ++i) {
		tmp_sigma += 2.f*residual[i]*exp(-pow(x_mu_sigma[i],2))/sqrt(M_PI)*(-x_mu_sigma[i]/sigma);
		tmp_mu += 2.f*residual[i]*exp(-pow(x_mu_sigma[i],2))/sqrt(M_PI)*(-1.f/sqrt(2.f*pow(sigma,2)));
	}
	slope_sigma = tmp_sigma;
	slope_mu = tmp_mu;
	return;
}

*/
