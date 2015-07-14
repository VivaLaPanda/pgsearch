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
};

DistFit::DistFit(){
}

vector< double > DistFit::calc_mean_var(vector< double > points){

	double sum  = 0, mu = 0, sigma_sq = 0, sdev = 0, dev = 0;
	sum = accumulate(points.begin(), points.end(), 0.0);
	mu = sum / points.size();
	sdev = inner_product(points.begin(), points.end(), points.begin(), 0.0);
	sigma_sq = sdev/points.size() - mu*mu;

	//cout << "Mean: " << mean << endl;
	//cout << "Variance: " << variance << endl;

//	for(int i = 0; i < points.size(); i++){
//		sum += points[i];
//	}
//	mu = sum/points.size();
//        for(int i = 1; i <= points.size(); ++i){
//		dev = (points[i] - mu)*(points[i] - mu);
//		sdev = sdev + dev;
//	}
//	sigma_sq = sdev / (points.size());

	vector< double > m_and_v;
	m_and_v.push_back(mu);
	m_and_v.push_back(sigma_sq);

	return m_and_v;
}


void DistFit::add_vector(vector<double> edge_vector){
	//log_new_edges = for_each(edge_vector.begin(), edge_vector.end(), &myfunction);

	double mu1, sigma1, mean1, variance1;

	if(log_edges.empty() == true){
		for(int i =0; i < edge_vector.size(); i++){
			vector< double > tmp, tmplog;
			tmp.push_back(edge_vector[i]);
			tmplog.push_back(log(edge_vector[i]));
			log_edges.push_back(tmplog);
			edges.push_back(tmp);
			mean_and_vars.push_back(calc_mean_var(edges[i]));
			log_mean_and_vars.push_back(calc_mean_var(log_edges[i]));
			mu1 = log_mean_and_vars[i][0];
			sigma1 = log_mean_and_vars[i][1];
			mean1 = exp(mu1 + pow(sigma1,2)/2);
			variance1 = exp(2*mu1 + pow(sigma1,2))*(exp(pow(sigma1,2))-1);
			log_mean_and_vars[i][0] = mean1;
			log_mean_and_vars[i][1] = variance1;
		}
	}

	else{
		for(int i = 0; i < edge_vector.size(); i++){
			log_edges[i].push_back(log(edge_vector[i]));
			edges[i].push_back(edge_vector[i]);
			mean_and_vars[i] = calc_mean_var(edges[i]);
			log_mean_and_vars[i] = calc_mean_var(log_edges[i]);
			mu1 = log_mean_and_vars[i][0];
			sigma1 = log_mean_and_vars[i][1];
			mean1 = exp(mu1 + pow(sigma1,2)/2);
			variance1 = exp(2*mu1 + pow(sigma1,2))*(exp(pow(sigma1,2))-1);
			log_mean_and_vars[i][0] = mean1;
			log_mean_and_vars[i][1] = variance1;
		}
	}
}

void DistFit::add_edge(double new_edge, int loc){
	double mu1, sigma1, mean1, variance1;
	if(loc == 0 && edges.empty() == true){
		vector< double > tmp, tmplog;
		tmp.push_back(new_edge);
		tmplog.push_back(log(new_edge));
		edges.push_back(tmp);
		log_edges.push_back(tmplog);
		mean_and_vars.push_back(calc_mean_var(edges[loc]));
		log_mean_and_vars.push_back(calc_mean_var(log_edges[loc]));
		mu1 = log_mean_and_vars[loc][0];
		sigma1 = log_mean_and_vars[loc][1];
		mean1 = exp(mu1 + pow(sigma1,2)/2);
		variance1 = exp(2*mu1 + pow(sigma1,2))*(exp(pow(sigma1,2))-1);
		log_mean_and_vars[loc][0] = mean1;
		log_mean_and_vars[loc][1] = variance1;
	}
	else{
		edges[loc].push_back(new_edge);
		log_edges[loc].push_back(log(new_edge));
		mean_and_vars[loc] = calc_mean_var(edges[loc]);
		log_mean_and_vars[loc] = calc_mean_var(log_edges[loc]);
		mu1 = log_mean_and_vars[loc][0];
		sigma1 = log_mean_and_vars[loc][1];
		mean1 = exp(mu1 + pow(sigma1,2)/2);
		variance1 = exp(2*mu1 + pow(sigma1,2))*(exp(pow(sigma1,2))-1);
		log_mean_and_vars[loc][0] = mean1;
		log_mean_and_vars[loc][1] = variance1;
	}
}





