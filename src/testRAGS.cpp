#include <iostream>
#include <utility>
#include <fstream>
#include <sstream>
#include <random>
#include "RAGS.h"

using namespace std ;
using namespace easymath ;

typedef unsigned long int ULONG ;
typedef pair<int, int> edge ;

// Read graph from configuration files
void DefineGraph(vector< XY > & locations, vector< edge > & edge_array, vector< vector<double> > & weights) 
{
	ifstream verticesFile("config_files/vertices0.txt") ;

	cout << "Reading vertices from file..." ;
	vector<double> v(2) ;
	string line ;
	while (getline(verticesFile,line))
	{
		stringstream lineStream(line) ;
		string cell ;
		int i = 0 ;
		while (getline(lineStream,cell,','))
		{
			v[i++] = atof(cell.c_str()) ;
		}
		locations.push_back(XY(v[0],v[1])) ;
	}
	cout << "complete.\n" ;
	
	ifstream edgesFile("config_files/edges0.txt") ;
	
	cout << "Reading edges from file..." ;
	vector<double> e(4) ;
	while (getline(edgesFile,line))
	{
		stringstream lineStream(line) ;
		string cell ;
		int i = 0 ;
		while (getline(lineStream,cell,','))
		{
			e[i++] = atof(cell.c_str()) ;
		}
		edge e_ind((int)e[0],(int)e[1]) ;
		vector<double> w(2) ;
		w[0] = e[2] ;
		w[1] = e[3] ;
		edge_array.push_back(e_ind) ;
		weights.push_back(w) ;
	}
	cout << "...complete.\n" ;
}

int main()
{
  // Create and store graph definition
  vector< XY > locations ;
  vector< edge > edge_array ;
  vector< vector<double> > weights ;
  
  DefineGraph(locations, edge_array, weights) ;
  
  // Create RAGS object
  RAGS * testRAGS = new RAGS(locations, edge_array, weights) ;
  
  // Store and set true edge costs
  default_random_engine generator(0) ;
  ULONG nEdges = testRAGS->GetGraph()->GetNumEdges() ;
  Edge ** eRAGS = testRAGS->GetGraph()->GetEdges() ;
  vector<double> trueWeights ;
  for (ULONG i = 0; i < nEdges; i++)
  {
    eRAGS[i]->SetTrueCost(generator) ;
    trueWeights.push_back(eRAGS[i]->GetTrueCost()) ;
  }
  
  // Set start and goal locations
  XY curLoc = locations[0] ;
  XY goal = locations[locations.size()-1] ;
  
  // Loop through execution weights
  XY nextLoc ;
  while (!(curLoc == goal))
  {
    cout << "Current traversal: (" << curLoc.x << "," << curLoc.y << ")" ;
    
    nextLoc = testRAGS->SearchGraph(curLoc, goal, trueWeights) ;
    
    cout << "-> (" << nextLoc.x << "," << nextLoc.y << ")" << endl ;
    
    curLoc = nextLoc ;
  }
  
  return 0 ;
}
