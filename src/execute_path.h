#ifndef EXECUTE_PATH
#define EXECUTE_PATH

#include <vector>
//#include <cmath>
#include <math.h>
#include <limits.h>
//#include <set>
//#include "node.h"
using namespace std;
const double pi = 3.14159265358979323846264338328 ;

double generateGaussianNoise(double mu, double sigma)
{
		const double epsilon = std::numeric_limits<double>::min();
		const double two_pi = 2.0*3.14159265358979323846;
		static double z0, z1;
		static bool generate;
		generate = !generate;
		if (!generate)
		   return z1 * sigma + mu;

		double u1, u2;
		do
		{
			u1 = rand() * (1.0 / RAND_MAX);
			u2 = rand() * (1.0 / RAND_MAX);
		}
		while ( u1 <= epsilon );
			 
			z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
			z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
			return z0 * sigma + mu;
}

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
bool ComputeImprovementProbability(Vertex* A, Vertex* B) //  double c_A0, double c_B0, vector<double> mu_A, vector<double> sig_A, vector<double> mu_B, vector<double> sig_B)
{ // The longer vectors should be set to B
	//rewrite this to work as a  comparitor
	double c_A0, c_B0;
	bool isBetter;
	vector<Node*> ANodes, BNodes;
	c_A0 = A->GetActualCost();
	c_B0 = B->GetActualCost();
	ANodes = A->GetNodes();
	BNodes = B->GetNodes();
	double max_3sig = ANodes[0]->GetMeanCost() + 3*ANodes[0]->GetVarCost() ;
	double min_3sig = ANodes[0]->GetMeanCost() - 3*ANodes[0]->GetVarCost() ;
	for (int i = 0; i < ANodes.size(); i++)
	{
		if (max_3sig < ANodes[i]->GetMeanCost()+3*ANodes[i]->GetVarCost())
			max_3sig = ANodes[i]->GetMeanCost()+3*ANodes[i]->GetVarCost() ;
		if (min_3sig > ANodes[i]->GetMeanCost()-3*ANodes[i]->GetVarCost())
			min_3sig = ANodes[i]->GetMeanCost()-3*ANodes[i]->GetVarCost() ;
	}
	for (int i = 0; i < BNodes.size(); i++)
	{
		if (max_3sig < BNodes[i]->GetMeanCost()+3*BNodes[i]->GetVarCost())
			max_3sig = BNodes[i]->GetMeanCost()+3*BNodes[i]->GetVarCost() ;
		if (min_3sig > BNodes[i]->GetMeanCost()-3*BNodes[i]->GetVarCost())
			min_3sig = BNodes[i]->GetMeanCost()-3*BNodes[i]->GetVarCost() ;
	}
	
	int n = 10000 ;
	vector<double> x = linspace(min_3sig,max_3sig,n) ;
	double dx = x[1]-x[0] ;
	double pImprove = 0.0 ;
	for (int k = 0; k < x.size(); k++)
	{
		double p_cAi = 0.0 ;
		for (int i = 0; i < ANodes.size(); i++)
		{
			double p_cA1 = (1/(ANodes[i]->GetVarCost()*sqrt(2*pi)))*exp(-(pow(x[k]-ANodes[i]->GetMeanCost(),2))/(2*pow(ANodes[i]->GetVarCost(),2))) ;
			double p_cA2 = 1.0 ;
			for (int j = 0; j < ANodes.size(); j++)
			{
				if (j != i)
					p_cA2 *= 0.5*erfc((x[k]-ANodes[j]->GetMeanCost())/(ANodes[j]->GetVarCost()*sqrt(2))) ;
			}
			p_cAi += p_cA1*p_cA2 ;
		}
		double p_cBi = 1.0 ;
		for (int i = 0; i < BNodes.size(); i++)
			p_cBi *= 0.5*erfc((x[k]-(c_B0-c_A0)-BNodes[i]->GetMeanCost())/(BNodes[i]->GetVarCost()*sqrt(2))) ;
		pImprove += (p_cAi)*(1-p_cBi)*dx ;
	}
	if(pImprove > .5){
		isBetter = true;
	}
	return isBetter;
}

bool TestAstarComparitor(Vertex* A, Vertex* B){
	if(A->GetActualCost() > B->GetActualCost()){
		return true;
	}
	else{
		return false;
	}
}

vector<double>  executePath(vector< Node*> GSPaths){
	//Initialize Memory
	cout << endl << "EXECUTING PATH" << endl;
	Vertex* goal;
	Vertex* cur_loc;
	vector < Node*> SGPaths, NewNodes;
	vector <Vertex*> vertices, nextVertices;
	Vertex* TmpVertex;
	
	//Reverse Paths because they are first given from goal to start
	//Outputs the reversed paths using the display paths function
	for(int i = 0; i < GSPaths.size(); i++){
		SGPaths.push_back(GSPaths[i]->ReverseList(0));
		cout << endl << "Path" << i << endl;
		SGPaths[i]->DisplayPath();
		cout << endl;
	}

	//Begin Iterating through paths, set start & goal and initialize the NewNodes vector to the path vector containing all possible start nodes.
	cout << "Iterating Through Paths" << endl;
	cur_loc = SGPaths[0]->GetVertex();
	goal = GSPaths[0]->GetVertex();
	NewNodes = SGPaths;
	ofstream executedCost;
	executedCost.open("executedPath.txt") ;
	cout << "Current Location: " << cur_loc->GetX() << ", " <<  cur_loc->GetY() << endl;
	double totalCost;
	//Begin execution
	while(cur_loc != goal){
		cout << "step" << endl;

		//For each node available next, parse into the vertices vector and assign that vertice an actual cost to it
		nextVertices.clear();
		for(int i = 0; i < NewNodes.size(); i++){
			cout << "NODE:" << i << " (" << NewNodes[i]->GetParent()->GetVertex()->GetX() << "," 
				<< NewNodes[i]->GetParent()->GetVertex()->GetY() << "), cost: " 
				<< NewNodes[i]->GetParent()->GetMeanCost() << ", var: " << NewNodes[i]->GetParent()->GetVarCost() << endl;
			nextVertices.push_back(NewNodes[i]->GetParent()->GetVertex());
			// Assign the parent nodes to the next vertex
			if(find(vertices.begin(), vertices.end(), NewNodes[i]->GetParent()->GetVertex()) == vertices.end()){ // this is rechecking old vertices that shouldn't need to be checked.
				cout << "found new vertex" << endl;
				NewNodes[i]->GetParent()->GetVertex()->SetNodes(NewNodes[i]->GetParent());
				vertices.push_back(NewNodes[i]->GetParent()->GetVertex());
				NewNodes[i]->GetParent()->GetVertex()->SetActualCost(generateGaussianNoise((NewNodes[i]->GetParent()->GetMeanCost()-NewNodes[i]->GetMeanCost()),(NewNodes[i]->GetParent()->GetVarCost()-NewNodes[i]->GetVarCost() )));
			}
		}

		//Sort the vertices by the improvement probability function
		//Set the cur location to the first vertice in the sorted vector
		//Get the nodes available at the new location
		cout << "Nodes List Size: " << NewNodes.size() << endl;
		sort(nextVertices.begin(), nextVertices.end(), ComputeImprovementProbability);
		for(int j = 0; j < vertices.size(); j++){
			cout << "Vertex list: " << nextVertices[j]->GetX() << ", "<< nextVertices[j]->GetY() << endl;
		}
		cur_loc = nextVertices[0];// this is where the bug is!
		cout << goal->GetX() << ", " << goal->GetY() << endl;
		NewNodes = cur_loc->GetNodes();
		//Status Prints
		cout << "Current Location: " << cur_loc->GetX() << ", " <<  cur_loc->GetY() << endl;
		totalCost +=  cur_loc->GetActualCost() ;
		cout << "Current Cost: " << cur_loc->GetActualCost() << endl;
	}
	cout << "TOTAL COST = " << totalCost << endl;

	executedCost.close();
	
	double min_path;
	int t;
	min_path = 100000; // make this max double	
	for(int i = 0; i < GSPaths.size(); i++){
		if(GSPaths[i]->GetMeanCost() < min_path){
			min_path = GSPaths[i]->GetMeanCost();
			t = i;
		}
	}
	cout << "Min Mean Path Cost " << GSPaths[t]->GetMeanCost() << endl;
	
	Node* loc;
	loc = SGPaths[t];
	double astarCost, cost;
	astarCost = 0;
	cost = 0;
	while(loc->GetVertex() != goal){
		if(find(vertices.begin(), vertices.end(), loc->GetParent()->GetVertex()) == vertices.end() ){ // make this true if vertex hasn't had cost set yet
			cout << "Found new Vertice" << endl;
			cost = generateGaussianNoise((loc->GetParent()->GetMeanCost()-loc->GetMeanCost()), (loc->GetParent()->GetVarCost()-loc->GetVarCost()));// rewrite to set the vertex cost based on edge to get there.
			loc->GetParent()->GetVertex()->SetActualCost(cost); //need to only generate cost if hasn't been generated in previous search
		}
		else{
			cout << "Already Defined Vertex" << endl;
			cout << "Current Cost: " << loc->GetParent()->GetVertex()->GetActualCost() << endl;
			cost = loc->GetParent()->GetVertex()->GetActualCost();
		}
		astarCost += cost;
		loc = loc->GetParent();
	}
	cout << "TOTAL ASTAR COST = " << astarCost << endl; // Need to print the total cost to a file.
	vector< double> costs;
	costs.push_back(totalCost);
	costs.push_back(astarCost);
	return costs;
}

#endif
