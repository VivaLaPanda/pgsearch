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


void executePath(vector< Node*> GSPaths){
	//Initialize Values
	cout << endl << "EXECUTING PATH" << endl;
	Vertex* goal_loc;
	Vertex* goal;
	Vertex* cur_loc;
	vector < Node*> SGPaths, NewNodes;
	vector <Vertex*> vertices;
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
	goal_loc = GSPaths[0]->GetVertex();
	NewNodes = SGPaths;

	//Begin execution
	while(cur_loc != goal){
		cout << "step" << endl;
		//For each node available next, parse into the vertices vector and assign that vertice an actual cost to it
		for(int i = 0; i < NewNodes.size(); i++){
			cout << "NODE:" << i << NewNodes[i]->GetParent() << endl;
			// Assign the parent nodes to the next vertex
			if(find(vertices.begin(), vertices.end(), NewNodes[i]->GetParent()->GetVertex()) == vertices.end()){
				cout << "found new vertex" << endl;
				NewNodes[i]->GetParent()->GetVertex()->SetNodes(NewNodes[i]->GetParent());
				vertices.push_back(NewNodes[i]->GetParent()->GetVertex());
				NewNodes[i]->GetParent()->GetVertex()->SetActualCost(generateGaussianNoise(NewNodes[i]->GetParent()->GetMeanCost(), NewNodes[i]->GetParent()->GetVarCost())); // I think this is right now? --- this line is wrong and needs to generate the actual path cost for each vertice based on the chosen path to get there.
			}
		}

		//Sort the vertices by the improvement probability function
		//Set the cur location to the first vertice in the sorted vector
		//Get the nodes available at the new location
		sort(vertices.begin(), vertices.end(), ComputeImprovementProbability);
		cur_loc = vertices[0];
		NewNodes = cur_loc->GetNodes();
		//Debug prints
		cout << "Current Location: " << cur_loc->GetX() << ", " <<  cur_loc->GetY() << endl;
		cout << "NewNodes Vertices: " << NewNodes[0]->GetVertex() << endl; //this is for current debugging
		cout << "NewNodes Vertices: " << NewNodes[1]->GetVertex() << endl;
		 
	}
}




//Old stuff saving for reference... PLEASE IGNORE
/*
	for(int i = 0; i < SGPaths.size(); i++){


		cout << "Path" <<  i  << endl;
		goal = GSPaths[i];
		goal_loc = goal->GetVertex();
		cur_loc = SGPaths[i];
		double gu = goal->GetMeanCost();
		double gv = goal->GetVarCost();

		//NextVertices.insert(SGPaths[i]->GetParent());
		//SGPaths = NextVertices[0]->GetNodes();


		cout << "Before While" << endl;
		while(cur_loc->GetVertex() != goal_loc){
			cout << "Step" << endl;
			cur_loc = cur_loc->GetParent();
			cout << "Current Location: " << cur_loc->GetVertex()->GetX() << ", " <<  cur_loc->GetVertex()->GetY() << endl;
			//Track the score
			cout << gu - cur_loc->GetMeanCost() << endl;
			cout << gv - cur_loc->GetVarCost() << endl << endl;
		}


	}
*/


/*
void executePath(vector<Node*> bestPaths)
{
	vector< Node *> reversedPaths(bestPaths.size());
	vector< Node *> vertices;
	cout << endl << "EXECUTING PATH" << endl;
	Vertex * cur_loc;
	vector< Node *> parents;

	while(cur_loc != bestPaths.end()){
		parents = []
		for(int i = 0; i < bestPaths.size(); i++){
			reversedPaths[i] =  bestPaths[i]->ReverseList(0);
			cout << endl << "Path" << i << endl;
			reversedPaths[i]->DisplayPath();
			cout << endl;
		}

		for(int i = 0; i < reversedPaths.size(); i++){
			cout << endl << "Path" << i << endl;
			reversedPaths[i]->DisplayPath();
			cout << endl;
			// store parents for next iter
			parents[i] = reversedPaths[i]->GetParent();// this is wrong and needs to be dependent on next loc vertices

			// find the mean cost of the first & last node

			reversedPaths[i]->GetMeanCost();
			bestPaths[i]->GetMeanCost();

			// Here add in while loop to iterate through each node along path

			if(find(vertices.begin(), vertices.end(), reversedPaths[i]->GetVertex())== vertices.end()){ //need to rewrite this correctly
				vertices.push_back(reveresedPaths[i]->GetVertex());
			}

		}

		//New paths/location
		reversedPaths = [];
		reversedPaths = parents;

		cur_loc = bestPaths[i]->GetVertex; // need to change this later
		cout << "CURRENT LOCATION: " << cur_loc->GetVertex();


		/* Get list of paths
		 * find all paths to get to different vertices, I think this needs to be only the next level because otherwise we will need to trim as we go becaues if we take a round about path to a vertex that doesn't get trimmed.
		 * create a vector of paths to vertices (aka vector of nodes that are at a vertice)
		 * Call *probability of best path*
		 * Take next action based on sorted (store in final path taken!)
		 i* Repeat
		 */

		/*
		 * set attributes of each vertice with score from them to the goal
		
	}
}

*/

#endif
