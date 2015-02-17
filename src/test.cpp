#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <queue>
#include <vector>
#include <functional>
#include <cmath>
#include <math.h>
#include <float.h>
#include <limits.h>
using namespace std ;

typedef unsigned long int ULONG ;

#include "config.h"
#include "vertex.h"
#include "edge.h"
#include "node.h"
#include "graph.h"
#include "queue.h"
#include "search.h"
#include "path.h"

void DefineGraph(vector< vector<double> > & vertVec, vector< vector<double> > & edgeVec) 
{
	ifstream verticesFile("sector_vertices.txt") ;

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
		vertVec.push_back(v) ;
	}
	cout << "complete.\n" ;
	
	//cout << "Size of vertices matrix: " << vertVec.size() << " x " << vertVec[0].size() << endl ;
	/*for (ULONG i = 0; i < vertVec.size(); i++)
		cout << "Vertex [" << i << "]: (" << vertVec[i][0] << "," << vertVec[i][1] << ")\n" ;*/
	
	ifstream edgesFile("sector_edges.txt") ;
	
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
		edgeVec.push_back(e) ;
	}
	cout << "...complete.\n" ;
	
	//cout << "Size of edges matrix: " << edgeVec.size() << " x " << edgeVec[0].size() << endl ;
	/*for (ULONG i = 0; i < edgeVec.size(); i++)
	{
		cout << "Edge [" << i << "]: (" << vertVec[edgeVec[i][0]][0] << "," << vertVec[edgeVec[i][0]][1] 
		<< ") to (" << vertVec[edgeVec[i][1]][0] << "," << vertVec[edgeVec[i][1]][1]
		<< "), cost: " << edgeVec[i][2] << ", var: " << edgeVec[i][3] << endl ;
	}*/
		
}

void DefineWorld(vector< vector<bool> > & obstacles, vector< vector<int> > & membership)
{
	ifstream obstaclesFile("obstacles.txt") ;
	
	cout << "Reading obstacle map from file..." ;
	vector<bool> obs(25) ;
	string line ;
	int lineCount = 0 ;
	while (getline(obstaclesFile,line))
	{
		stringstream lineStream(line) ;
		string cell ;
		int i = 0 ;
		while (getline(lineStream,cell,','))
		{
			if (atoi(cell.c_str()) == 1)
				obs[i++] = true ;
			else
				obs[i++] = false ;
		}
		obstacles.push_back(obs) ;
	}
	cout << "complete.\n" ;
	
	ifstream sectorsFile("membership.txt") ;
	
	cout << "Reading sector map from file..." ;
	while (getline(sectorsFile,line))
	{
		stringstream lineStream(line) ;
		string cell ;
		vector<int> mem ;
		while (getline(lineStream,cell,','))
		{
			mem.push_back(atoi(cell.c_str())) ;
		}
		membership.push_back(mem) ;
	}
	cout << "complete.\n" ;
}

int main()
{
	cout << "Test program...\n" ;
	
	/*// Testing on a 4 or 8 connected grid
	double xMin = 0.0 ;
	double xInc = 1.0 ;
	int lenX = 100 ;
	double yMin = 0.0 ;
	double yInc = 1.0 ;
	int lenY = 100 ;
	vector<double> xGrid(lenX) ;
	vector<double> yGrid(lenY) ;
	for (int i = 0; i < lenX; i++)
		xGrid[i] = xMin + i*xInc ;
	
	for (int i = 0; i < lenY; i++)
		yGrid[i] = yMin + i*yInc ;
		
	Graph * testGraph = new Graph(xGrid,yGrid,8) ;*/
	
	// Create vectors and edges in graph from text files
	vector< vector<double> > vertVec ;
	vector< vector<double> > edgeVec ;
	DefineGraph(vertVec, edgeVec) ;
	
	// Create graph
	Graph * testGraph = new Graph(vertVec, edgeVec) ;
	Vertex * sourceSec = testGraph->GetVertices()[0] ; // top-left sector
	Vertex * goalSec = testGraph->GetVertices()[testGraph->GetNumVertices()-1] ; // bottom-right sector
	
	//Create search object and perform path search from source to goal
	cout << "Creating search object..." ;
	Search * testSearch = new Search(testGraph, sourceSec, goalSec) ;
	cout << "complete.\n" ;
	
	cout << "Performing path search from (" <<  sourceSec->GetX() << "," << sourceSec->GetY() << ") to (" ;
	cout << goalSec->GetX() << "," << goalSec->GetY() << ")...\n" ;
	pathOut pType = BEST ;
	vector<Node *> bestPaths = testSearch->PathSearch(pType) ;
	cout << "Path search complete.\n" ;
	
	if (bestPaths.size() != 0)
	{
		for (ULONG i = 0; i < (ULONG)bestPaths.size(); i++)
		{
			cout << "Path " << i << endl ;
			bestPaths[i]->DisplayPath() ;
		}
	}
	
	//Read in membership and obstacle text files
	vector< vector<bool> > obstacles ;
	vector< vector<int> > membership ;
	DefineWorld(obstacles, membership) ;
	
	//Create path object to plan low level path
	Vertex * sourcePath = new Vertex(2.0,3.0) ;
	Vertex * goalPath = new Vertex(48.0,17.0) ;
	cout << "Performing low level path search from (" << sourcePath->GetX() << "," ;
	cout << sourcePath->GetY() << ") to (" << goalPath->GetX() << "," ;
	cout << goalPath->GetY() << ")...\n" ;
	Path * testPath = new Path(obstacles, membership, bestPaths[0], sourcePath, goalPath) ;
	Node * lowLevelPath = testPath->ComputePath(8) ;
	cout << "Low level path search complete.\n" ;
	
	delete testGraph ;
	testGraph = 0 ;
	delete testSearch ;
	testSearch = 0 ;
	delete sourcePath ;
	sourcePath = 0 ;
	delete goalPath ;
	goalPath = 0 ;
	delete testPath ;
	testPath = 0 ;
	
	return 0 ;
}
