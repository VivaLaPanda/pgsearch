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
	private:
		double mu, sigma, mean, variance;
		vector< vector < double > > edges;
		vector< vector < double > > log_edges;
		vector< vector< double > > mean_and_vars;
		vector< vector< double > > log_mean_and_vars;
	public:
		DistFit();
		vector< double > calc_mean_var(vector< double >);
		void add_edge(double, int);
		void add_vector( vector< double >);
		vector< vector< double> > get_normMV() const {return mean_and_vars; };
		vector< vector< double> > get_logMV() const {return log_mean_and_vars;};
		double get_variance() const {return variance; };
};

DistFit::DistFit(){
	mean = 0;
	variance = 0;
	mu = 0;
	sigma = 0;
	cout << "TEST3";
}

vector< double > DistFit::calc_mean_var(vector< double > points){

	double sum = accumulate(points.begin(), points.end(), 0.0);
	double mu = sum / points.size();
	double sumSquares = inner_product(points.begin(), points.end(), points.begin(), 0.0);
	double sigma = sumSquares/points.size() - mu*mu;

	//cout << "Mean: " << mean << endl;
	//cout << "Variance: " << variance << endl;

	vector< double > m_and_v;
	m_and_v.push_back(mu);
	m_and_v.push_back(sigma);

	cout << "TEST4";
	return m_and_v;
}


void  DistFit::add_vector(vector<double> edge_vector){
	cout << "TEST5";
	//log_new_edges = for_each(edge_vector.begin(), edge_vector.end(), &myfunction);
	for(int i; i < edge_vector.size(); i++){
		log_edges[i].push_back(log(edge_vector[i]));
		edges[i].push_back(edge_vector[i]);
		mean_and_vars.push_back(calc_mean_var(edges[i]));
		log_mean_and_vars.push_back(calc_mean_var(log_edges[i]));
		mu = log_mean_and_vars[i][0];
		sigma = log_mean_and_vars[i][1];
		mean= exp(mu + pow(sigma,2)/2);
		variance = exp(2*mu + pow(sigma,2))*(exp(pow(sigma,2)-1));
		vector< double> tmp;
		log_mean_and_vars[i][0] = mean;
		log_mean_and_vars[i][1] = variance;
	}
}

void DistFit::add_edge(double new_edge, int loc){
	edges[loc].push_back(new_edge);
	log_edges[loc].push_back(log(new_edge));
	mean_and_vars.push_back(calc_mean_var(edges[loc]));
	log_mean_and_vars.push_back(calc_mean_var(log_edges[loc]));
	mu = log_mean_and_vars[loc][0];
	sigma = log_mean_and_vars[loc][1];
	mean= exp(mu + pow(sigma,2)/2);
	variance = exp(2*mu + pow(sigma,2))*(exp(pow(sigma,2)-1));
	vector< double> tmp;
	tmp.push_back(mean);
	tmp.push_back(variance);
	log_mean_and_vars[loc] = tmp;
}
