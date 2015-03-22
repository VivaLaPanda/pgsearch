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
#include <stdlib.h>
#include <typeinfo>

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
#include "execute_path.h"

void DefineGraph(vector< vector<double> > & vertVec, vector< vector<double> > & edgeVec) 
{
	ifstream verticesFile("sector_vertices1.txt") ;

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

	ifstream edgesFile("sector_edges1.txt") ;

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

vector< vector< double > > makeVertices(double x, double y){
	vector< vector< double > > vertices;
	vector< double > tmp;
	srand (time(NULL));
	int numVerts;
	double vertx, verty;
	int xx = x;
	int yy = y;
	double testx, testy;
	numVerts = 50; // number of vertices to generate wihin specified x, y area
	tmp.push_back(0);
	tmp.push_back(0);

	vertices.push_back(tmp);
	for(int i = 0; i < numVerts; i++){
		vertx = rand() % xx;
		verty = rand() % yy;
		testx = vertx;
		testy = verty;
		vector< double > row;
		row.push_back(testx);
		row.push_back(testy);
		vertices.push_back(row);
	}

//	for(int i = 0; i < vertices.size(); i++){
//		for(int j = 0; vertices[0].size(); j++){
//			cout << vertices[i][j] << endl;
//		}
//	}
	return vertices;
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

	// Uncomment to create vectors and edges in graph from text files
	//vector< vector<double> > vertVec ;
	//vector< vector<double> > edgeVec ;
	//DefineGraph(vertVec, edgeVec) ;

	// Create graph
	// Need to  make a vector of vertices or adapt the graph.h file to generate them automatically given a x,y area.
	vector< vector< double > > vertVec2;
	double x,y, radius;
	x = 50;
	y = 50;
	cout << "Generating Random Vertices in " << x << " by " << y << endl;
	vertVec2 = makeVertices(x,y);
	radius= 10.0;// need to rewrite into the formula that determines radius automatically based on max vertex distance
	cout << "Connecting with radius " << radius << endl;
	//Graph * testGraph = new Graph(vertVec, edgeVec) ;
	Graph * testGraph = new Graph(vertVec2, radius);
	Vertex * sourceSec = testGraph->GetVertices()[0] ; // top-left sector
	Vertex * goalSec = testGraph->GetVertices()[testGraph->GetNumVertices()-1] ; // bottom-right sector

	//Create search object and perform path search from source to goal
	cout << "Creating search object..." ;
	Search * testSearch = new Search(testGraph, sourceSec, goalSec) ;
	cout << "complete.\n" ;

	cout << "Performing path search from (" <<  sourceSec->GetX() << "," << sourceSec->GetY() << ") to (" ;
	cout << goalSec->GetX() << "," << goalSec->GetY() << ")...\n" ;
	pathOut pType = ALL ;
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

	executePath(bestPaths);




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
	return 0 ;
}
