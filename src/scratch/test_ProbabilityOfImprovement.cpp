#include <iostream>
//#include <algorithm>
//#include <math.h>
#include <float.h>
#include <random>
//#include <time.h>
using namespace std ;

const double pi = 3.14159265358979323846264338328 ;

vector<double> linspace(double a, double b, int n)
{
	vector<double> array ;
	double step = (b-a)/(n-1) ;
	while (a<=b)
	{
		array.push_back(a) ;
		a += step ;
	}
	return array ;
}

/*double IntNormDist(double mu, double sig, double a)
{
	int n = 1000 ;
	vector<double> x = linspace(a,mu+3*sig,n) ;
	double dx = (x[1]-x[0]) ;
	
	double integral = 0.0 ;
	for (int i = 0; i < x.size(); i++)
	{
		double fac = 1/(sig*sqrt(2*pi)) ;
		double ind = -(pow(x[i]-mu,2))/(2*pow(sig,2)) ;
		integral += fac*exp(ind)*dx ;
	}
	return integral ;
}*/

double ComputeImprovementProbability(double c_A0, double c_B0, vector<double> mu_A, vector<double> sig_A, vector<double> mu_B, vector<double> sig_B)
{
	double max_3sig = mu_A[0] + 3*sig_A[0] ;
	double min_3sig = mu_A[0] - 3*sig_A[0] ;
	for (int i = 0; i < mu_A.size(); i++)
	{
		if (max_3sig < mu_A[i]+3*sig_A[i])
			max_3sig = mu_A[i]+3*sig_A[i] ;
		if (min_3sig > mu_A[i]-3*sig_A[i])
			min_3sig = mu_A[i]-3*sig_A[i] ;
	}
	for (int i = 0; i < mu_B.size(); i++)
	{
		if (max_3sig < mu_B[i]+3*sig_B[i])
			max_3sig = mu_B[i]+3*sig_B[i] ;
		if (min_3sig > mu_B[i]-3*sig_B[i])
			min_3sig = mu_B[i]-3*sig_B[i] ;
	}
	
	int n = 10000 ;
	vector<double> x = linspace(min_3sig,max_3sig,n) ;
	double dx = x[1]-x[0] ;
	double pImprove = 0.0 ;
	for (int k = 0; k < x.size(); k++)
	{
		double p_cAi = 0.0 ;
		for (int i = 0; i < mu_A.size(); i++)
		{
			double p_cA1 = (1/(sig_A[i]*sqrt(2*pi)))*exp(-(pow(x[k]-mu_A[i],2))/(2*pow(sig_A[i],2))) ;
			double p_cA2 = 1.0 ;
			for (int j = 0; j < mu_A.size(); j++)
			{
				if (j != i)
					p_cA2 *= 0.5*erfc((x[k]-mu_A[j])/(sig_A[j]*sqrt(2))) ;
			}
			p_cAi += p_cA1*p_cA2 ;
		}
		double p_cBi = 1.0 ;
		for (int i = 0; i < mu_B.size(); i++)
			p_cBi *= 0.5*erfc((x[k]-(c_B0-c_A0)-mu_B[i])/(sig_B[i]*sqrt(2))) ;
		pImprove += (p_cAi)*(1-p_cBi)*dx ;
	}
	return pImprove;
}

int main()
{
	
	random_device generator ;
	
	uniform_real_distribution<double> mu_dist(0.0,10.0) ;
	uniform_real_distribution<double> sig_dist(0.0,10.0) ;
	
	vector<double> mu_A ;
	vector<double> sig_A ;

	int m = 3 ;
	for (int i = 0; i < m; i++)
	{
		mu_A.push_back(mu_dist(generator)) ;
		sig_A.push_back(sig_dist(generator)) ;
	}

	vector<double> mu_B ;
	vector<double> sig_B ;
	
	int n = 10 ;
	for (int i = 0; i < n; i++)
	{
		mu_B.push_back(mu_dist(generator)) ;
		sig_B.push_back(sig_dist(generator)) ;
	}
	
	double c_A0 = mu_dist(generator) ;
	double c_B0 = mu_dist(generator) ;

	clock_t t_fun1 = clock() ;
	double pImproveA_B = ComputeImprovementProbability(c_B0,c_A0,mu_B,sig_B,mu_A,sig_A) ;
	t_fun1 = clock() - t_fun1 ;

	cout << "Probability of improvement: " << pImproveA_B << endl ;
	cout << "Tics elapsed: " << t_fun1 << ", " << ((float)t_fun1)/CLOCKS_PER_SEC << "s" << endl ;

	clock_t t_fun0 ;
	t_fun0 = clock() ;
	double pImproveB_A = ComputeImprovementProbability(c_A0,c_B0,mu_A,sig_A,mu_B,sig_B) ;
	t_fun0 = clock() - t_fun0 ;

	cout << "Probability of improvement: " << pImproveB_A << endl ;
	cout << "Tics elapsed: " << t_fun0 << ", " << ((float)t_fun0)/CLOCKS_PER_SEC << "s" << endl ;

	n = 10000 ;
	vector<double> A_cost(n,c_A0) ;
	vector<double> B_cost(n,c_B0) ;
	
	double pBoverA = 0.0 ;
	clock_t t_MC ;
	t_MC = clock() ;
	for (int k = 0; k < n; k++)
	{
		double best_A = DBL_MAX ;
		for (int i = 0; i < mu_A.size(); i++)
		{
			normal_distribution<double> distribution(mu_A[i],sig_A[i]) ;
			double c_Ai = distribution(generator) ;
				if (c_Ai < best_A)
					best_A = c_Ai ;
		}
		A_cost[k] += best_A ;

		double best_B = DBL_MAX ;
		for (int i = 0; i < mu_B.size(); i++)
		{
			normal_distribution<double> distribution(mu_B[i],sig_B[i]) ;
			double c_Bi = distribution(generator) ;
				if (c_Bi < best_B)
					best_B = c_Bi ;
		}
		B_cost[k] += best_B ;

		if (B_cost[k]<=A_cost[k])
			pBoverA++ ;
	}
	t_MC = clock() - t_MC ;

	double pBA = pBoverA/n ;
	cout << "Monte Carlo probability: " << pBA << endl ;
	cout << "Tics elapsed: " << t_MC << ", " << ((float)t_MC)/CLOCKS_PER_SEC << "s" << endl ;

	return 0 ;
}
