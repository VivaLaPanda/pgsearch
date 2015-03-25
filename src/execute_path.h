#include <vector>
#include <math.h>
#include <limits.h>

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

bool ComputeImprovementProbability(Vertex * A, Vertex * B)
{
	vector<Node *> ANodes = A->GetNodes();
	vector<Node *> BNodes = B->GetNodes();
	
	double c_A0 = A->GetCTC() ;
	double c_B0 = B->GetCTC() ;
	double max_3sig = ANodes[0]->GetMeanCTG() + 3*ANodes[0]->GetVarCTG() ;
	double min_3sig = ANodes[0]->GetMeanCTG() - 3*ANodes[0]->GetVarCTG() ;
	for (int i = 0; i < ANodes.size(); i++)
	{
		if (max_3sig < ANodes[i]->GetMeanCTG()+3*ANodes[i]->GetVarCTG())
			max_3sig = ANodes[i]->GetMeanCTG()+3*ANodes[i]->GetVarCTG() ;
		if (min_3sig > ANodes[i]->GetMeanCTG()-3*ANodes[i]->GetVarCTG())
			min_3sig = ANodes[i]->GetMeanCTG()-3*ANodes[i]->GetVarCTG() ;
	}
	for (int i = 0; i < BNodes.size(); i++)
	{
		if (max_3sig < BNodes[i]->GetMeanCTG()+3*BNodes[i]->GetVarCTG())
			max_3sig = BNodes[i]->GetMeanCTG()+3*BNodes[i]->GetVarCTG() ;
		if (min_3sig > BNodes[i]->GetMeanCTG()-3*BNodes[i]->GetVarCTG())
			min_3sig = BNodes[i]->GetMeanCTG()-3*BNodes[i]->GetVarCTG() ;
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
			double mu_Ai = ANodes[i]->GetMeanCTG() ;
			double sig_Ai = ANodes[i]->GetVarCTG() ;
			double p_cA1 = (1/(sig_Ai*sqrt(2*pi)))*exp(-(pow(x[k]-mu_Ai,2))/(2*pow(sig_Ai,2))) ;
			double p_cA2 = 1.0 ;
			for (int j = 0; j < ANodes.size(); j++)
			{
				double mu_Aj = ANodes[j]->GetMeanCTG() ;
				double sig_Aj = ANodes[j]->GetVarCTG() ;
				if (j != i)
					p_cA2 *= 0.5*erfc((x[k]-mu_Aj)/(sig_Aj*sqrt(2))) ;
			}
			p_cAi += p_cA1*p_cA2 ;
		}
		double p_cBi = 1.0 ;
		for (int i = 0; i < BNodes.size(); i++)
		{
			double mu_Bi = BNodes[i]->GetMeanCTG() ;
			double sig_Bi = BNodes[i]->GetVarCTG() ;
			p_cBi *= 0.5*erfc((x[k]-(c_B0-c_A0)-mu_Bi)/(sig_Bi*sqrt(2))) ;
		}
		pImprove += (p_cAi)*(1-p_cBi)*dx ;
	}
	
	return (pImprove<=0.5);
}

bool GreedyComparison(Vertex* A, Vertex* B)
{
	return (A->GetCTC() < B->GetCTC()) ;
}

void AssignTrueEdgeCosts(Graph * searchGraph, int seed)
{
	Edge ** allEdges = searchGraph->GetEdges() ;
	ULONG numEdges = searchGraph->GetNumEdges() ;
	default_random_engine generator(seed) ;
	cout << "Using generator(" << generator() << ")...\n" ;
	for (ULONG i = 0; i < numEdges; i++)
		allEdges[i]->SetTrueCost(generator) ;
}

void SetTrueEdgeCosts(Vertex * v1, vector<Vertex *> v2, Graph * searchGraph)
{
	// Find all outgoing edges
	Edge ** allEdges = searchGraph->GetEdges() ;
	ULONG numEdges = searchGraph->GetNumEdges() ;
	vector<Edge *> connectedEdges ;
	for (int j = 0; j < numEdges; j++)
	{
		if (allEdges[j]->GetVertex1()->GetX() == v1->GetX() && 
			allEdges[j]->GetVertex1()->GetY() == v1->GetY())
			connectedEdges.push_back(allEdges[j]) ;
	}
	
	// Set cost-to-come of v2 vertices as true cost of edge traversal
	for (int i = 0; i < v2.size(); i++)
	{
		for (int j = 0; j < connectedEdges.size(); j++)
		{
			if (connectedEdges[j]->GetVertex2()->GetX() == v2[i]->GetX() && 
				connectedEdges[j]->GetVertex2()->GetY() == v2[i]->GetY())
			{
				v2[i]->SetCTC(connectedEdges[j]->GetTrueCost()) ;
				break ;
			}
		}
	}
}

vector<double> executePath(vector< Node*> GSPaths, Graph * searchGraph)
{
	vector <Node *> SGPaths ; // store all paths, start to goal
	vector <Node *> newNodes, tmpNodes ; // store remaining paths, start to goal
	vector <Vertex *> nextVerts ; // store connected vertices
	vector<double> allCosts ; // store costs of PGsearch, A*, D* searches
	double totalCost = 0 ; // store total cost
	
	// Set linked list from start node to goal node
	for(int i = 0; i < GSPaths.size(); i++)
		SGPaths.push_back(GSPaths[i]->ReverseList(0));
	
	// Compute cost-to-go for all path nodes
	for(int i = 0; i < SGPaths.size(); i++)
		SGPaths[i]->SetCTG(GSPaths[i]->GetMeanCost(),GSPaths[i]->GetVarCost()) ;
	
	// Step through paths
	cout << "Stepping through paths..." << endl ;
	Vertex * curLoc = SGPaths[0]->GetVertex() ;
	Vertex * goal = GSPaths[0]->GetVertex() ;
	newNodes = SGPaths ;

	// Assign nodes to first vertex
	curLoc->SetNodes(newNodes) ;
	
	while (curLoc->GetX() != goal->GetX() || curLoc->GetY() != goal->GetY())
	{
		cout << "Current Location: (" << curLoc->GetX() << "," <<  curLoc->GetY() << ")\n" ;
		
		// Extract nodes of current vertex
		newNodes = curLoc->GetNodes() ;
		
		// Identify next vertices
		for (int i = 0; i < newNodes.size(); i++)
		{
			bool newVert = true ;
			for (int j = 0; j < nextVerts.size(); j++)
			{
				if (nextVerts[j]->GetX() == newNodes[i]->GetParent()->GetVertex()->GetX() &&
					nextVerts[j]->GetY() == newNodes[i]->GetParent()->GetVertex()->GetY())
				{
					newVert = false ;
					break ;
				}
			}
			if (newVert)
				nextVerts.push_back(newNodes[i]->GetParent()->GetVertex()) ;
		}
		
		// Identify next vertex path nodes
		tmpNodes.clear() ;
		for (int i = 0; i < nextVerts.size(); i++)
		{
			for (int j = 0; j < newNodes.size(); j++)
			{
				if (nextVerts[i]->GetX() == newNodes[j]->GetParent()->GetVertex()->GetX() &&
					nextVerts[i]->GetY() == newNodes[j]->GetParent()->GetVertex()->GetY())
					tmpNodes.push_back(newNodes[j]->GetParent()) ;
			}
			nextVerts[i]->SetNodes(tmpNodes) ;
		}
		
		// Set cost-to-come for next vertices
		SetTrueEdgeCosts(curLoc, nextVerts, searchGraph) ;
		
		// Rank next vertices according to probability of improvement
		sort(nextVerts.begin(),nextVerts.end(),ComputeImprovementProbability) ;
		
		// Move to best vertex and log the cost
		curLoc = nextVerts[0] ;
		totalCost += nextVerts[0]->GetCTC() ;
		
		// Clear vectors for next step
		newNodes.clear() ;
		nextVerts.clear() ;
	}
	
	cout << "Goal vertex (" << curLoc->GetX() << "," <<  curLoc->GetY() << ") reached! "
		<< "Total cost: " << totalCost << endl ;
	allCosts.push_back(totalCost) ;
	
	newNodes.clear() ;
	nextVerts.clear() ;
	
	/**********************************************************************************************/
	// Traverse A* path (lowest mean cost path, no adaptivity)
	double lowestCost = GSPaths[0]->GetMeanCost() ;
	double lowestCostInd = 0 ;
	for (int i = 0; i < GSPaths.size(); i++)
	{
		if (GSPaths[i]->GetMeanCost() < lowestCost)
		{
			lowestCost = GSPaths[i]->GetMeanCost() ;
			lowestCostInd = i ;
		}
	}
	
	cout << "Traversing A* path...\n" ;
	totalCost = 0 ; // reset path cost
	Node * curNode = SGPaths[lowestCostInd] ;
	curLoc = curNode->GetVertex() ;
	while (curNode->GetParent())
	{
		// Identify edge to traverse
		nextVerts.push_back(curNode->GetParent()->GetVertex()) ;
		
		cout << "Current Location: (" << curLoc->GetX() << "," <<  curLoc->GetY() << ")\n" ;
		
		// Set cost-to-come for next vertex
		SetTrueEdgeCosts(curLoc, nextVerts, searchGraph) ;
		
		// Move to next vertex and log the cost
		curNode = curNode->GetParent() ;
		curLoc = curNode->GetVertex() ;
		totalCost += nextVerts[0]->GetCTC() ;
		
		// Clear vectors for next step
		newNodes.clear() ;
		nextVerts.clear() ;
	}
	
	cout << "Goal vertex (" << curLoc->GetX() << "," <<  curLoc->GetY() << ") reached! "
		<< "Total cost: " << totalCost << endl ;
	allCosts.push_back(totalCost) ;
	
	/**********************************************************************************************/
	// Traverse greedy path (lowest immediate cost, searching only non-dominated paths)
	cout << "Traversing greedy path...\n" ;
	curLoc = SGPaths[0]->GetVertex() ; // reset path search
	newNodes = SGPaths ; // reset path search
	totalCost = 0 ; // reset path cost
	
	// Assign nodes to first vertex
	curLoc->SetNodes(newNodes) ;
	
	while (curLoc->GetX() != goal->GetX() || curLoc->GetY() != goal->GetY())
	{
		cout << "Current Location: (" << curLoc->GetX() << "," <<  curLoc->GetY() << ")\n" ;
		
		// Extract nodes of current vertex
		newNodes = curLoc->GetNodes() ;
		
		// Identify next vertices
		for (int i = 0; i < newNodes.size(); i++)
		{
			if (!newNodes[i]->GetParent())
				continue ;
			bool newVert = true ;
			for (int j = 0; j < nextVerts.size(); j++)
			{
				
				if (nextVerts[j]->GetX() == newNodes[i]->GetParent()->GetVertex()->GetX() &&
					nextVerts[j]->GetY() == newNodes[i]->GetParent()->GetVertex()->GetY())
				{
					newVert = false ;
					break ;
				}
			}
			if (newVert)
				nextVerts.push_back(newNodes[i]->GetParent()->GetVertex()) ;
		}
		
		// Identify next vertex path nodes
		tmpNodes.clear() ;
		for (int i = 0; i < nextVerts.size(); i++)
		{
			for (int j = 0; j < newNodes.size(); j++)
			{
				if (!newNodes[j]->GetParent())
					continue ;
				if (nextVerts[i]->GetX() == newNodes[j]->GetParent()->GetVertex()->GetX() &&
					nextVerts[i]->GetY() == newNodes[j]->GetParent()->GetVertex()->GetY())
					tmpNodes.push_back(newNodes[j]->GetParent()) ;
			}
			nextVerts[i]->SetNodes(tmpNodes) ;
		}
		
		// Set cost-to-come for next vertices
		SetTrueEdgeCosts(curLoc, nextVerts, searchGraph) ;
		
		// Rank next vertices according to probability of improvement
		sort(nextVerts.begin(),nextVerts.end(),GreedyComparison) ;
		
		// Move to best vertex and log the cost
		curLoc = nextVerts[0] ;
		totalCost += nextVerts[0]->GetCTC() ;
		
		// Clear vectors for next step
		newNodes.clear() ;
		nextVerts.clear() ;
	}
	
	cout << "Goal vertex (" << curLoc->GetX() << "," <<  curLoc->GetY() << ") reached! "
		<< "Total cost: " << totalCost << endl ;
	allCosts.push_back(totalCost) ;
	
	return allCosts ;
	
}
