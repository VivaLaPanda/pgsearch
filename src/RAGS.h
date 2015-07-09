#include <random>
#include <utility>
#include <algorithm>
#include <time.h>
#include <float.h>
#include <limits.h>

#include "../Math/easymath.h"

using namespace easymath ;
using namespace std ;

typedef unsigned long int ULONG ;
enum searchType {ASTAR, DIJKSTRA} ; // BREADTH, DEPTH
enum heuristic {ZERO, MANHATTAN, EUCLIDEAN} ;
enum pathOut {BEST,ALL} ;
enum nodeType {SOURCE, OTHER} ;

const double pi = 3.14159265358979323846264338328 ;

// Search configuration
searchType SEARCH_TYPE = ASTAR ;
heuristic HEURISTIC = ZERO ;
pathOut PSET = ALL ;

// Vertex class to contain location of vertices
class Node ;
class Vertex
{
	public:
		Vertex(double x, double y):
		  itsX(x), itsY(y) {}
		~Vertex() {}
		
		double GetX() const {return itsX ;}
		void SetX(double x) {itsX = x ;}
		double GetY() const {return itsY ;}
		void SetY(double y) {itsY = y ;}
		double GetCTC()const{return itsCTC ;}
		void SetCTC(double cost_to_come) {itsCTC = cost_to_come ;}
		void SetNodes(vector<Node *> NewNodes) {itsNodes = NewNodes ;}
		vector<Node *> GetNodes() {return itsNodes ;}
		double SetActualCost(double ac) {itsActualCost = ac;}
		double GetActualCost() const {return itsActualCost ;}
		
		friend bool operator==(const Vertex &lhs, const Vertex &rhs){
			return lhs.GetX()==rhs.GetX() && lhs.GetY()==rhs.GetY();
		}

	private:
		double itsActualCost ;
		double itsX ;
		double itsY ;
		double itsCTC ;
		vector<Node *> itsNodes ;
};

// Edge class to contain mean and variance of cost along an edge
class Edge
{
	public:
		Edge(Vertex * v1, Vertex * v2, double cost, double var):
		  itsVertex1(v1), itsVertex2(v2), itsMeanCost(cost), itsVarCost(var), itsMeanSearch(cost), itsVarSearch(var) {}
		~Edge()
		{
	    delete itsVertex1 ;
	    itsVertex1 = 0 ;
	    delete itsVertex2 ;
	    itsVertex2 = 0 ;
    }
		
		Vertex * GetVertex1() const {return itsVertex1 ;}
		Vertex * GetVertex2() const {return itsVertex2 ;}
		double GetMeanCost() const {return itsMeanCost ;}
		void SetMeanCost(double cost) {itsMeanCost = cost ;}
		double GetVarCost() const {return itsVarCost ;}
		void SetVarCost(double var) {itsVarCost = var ;}
		double GetTrueCost() {return itsTrueCost ;}
		void SetTrueCost(double cost) {itsTrueCost = cost ;}
		double GetMeanSearch() const {return itsMeanSearch ;}
		void SetMeanSearch(double cost) {itsMeanSearch = cost ;}
		double GetVarSearch() const {return itsVarSearch ;}
		void SetVarSearch(double var) {itsVarSearch = var ;}
		void SetTrueCost(default_random_engine generator)
		{
	    double diff_x = itsVertex1->GetX() - itsVertex2->GetX() ;
	    double diff_y = itsVertex1->GetY() - itsVertex2->GetY() ;
	    double diff = sqrt(pow(diff_x,2) + pow(diff_y,2)) ;
	    normal_distribution<double> distribution(itsMeanCost,itsVarCost) ;
	    itsTrueCost = distribution(generator) ;
	    if (itsTrueCost < diff)
		    itsTrueCost = diff ;
    }
	private:
		Vertex * itsVertex1 ;
		Vertex * itsVertex2 ;
		double itsMeanCost ;
		double itsVarCost ;
		double itsTrueCost ;
		double itsMeanSearch ; // Actual value used in search
		double itsVarSearch ; // Actual value used in search
} ;

// Graph class to create and store graph structure
class Graph
{
	public:
	  typedef pair<int, int> edge ;
	  Graph(vector<XY> &locations, vector<edge> &edge_array, vector< vector<double> > &weights)
	  {
	    itsVertices = GenerateVertices(locations) ;
	    itsEdges = GenerateEdges(edge_array, weights) ;
    }
    
		~Graph()
		{
	    delete [] itsVertices ;
	    itsVertices = 0;
	    delete [] itsEdges ;
	    itsEdges = 0 ;
    }
		
		Vertex ** GetVertices() const {return itsVertices ;}
		Edge ** GetEdges() const {return itsEdges ;}
		ULONG GetNumVertices() const {return numVertices ;}
		ULONG GetNumEdges() const {return numEdges ;}
		
		vector<Edge *> GetNeighbours(XY v)
		{
      Vertex * pV ;
      for (ULONG i = 0; i < numVertices; i++)
      {
        if ( v.x == itsVertices[i]->GetX() && v.y == itsVertices[i]->GetY() )
          pV = itsVertices[i] ;
      }
      return GetNeighbours(pV) ;
		}
		
		vector<Edge *> GetNeighbours(Vertex * v)
		{
	    vector<Edge *> neighbours(numEdges) ;
	    ULONG k = 0 ;
	
	    for (ULONG i = 0; i < numEdges; i++)
	    {
		    double x1 = itsEdges[i]->GetVertex1()->GetX() ;
		    double y1 = itsEdges[i]->GetVertex1()->GetY() ;
		    if (x1 == v->GetX() && y1 == v->GetY())
		    {
			    neighbours[k] = itsEdges[i] ;
			    k++ ;
		    }
	    }
	
	    neighbours.resize(k) ;
	    return neighbours ;
    }
    
	private:
		Vertex ** itsVertices ;
		Edge ** itsEdges ;
		ULONG numVertices ;
		ULONG numEdges ;
		
		Vertex ** GenerateVertices(vector<XY> &vertices)
		{
      numVertices = (ULONG)vertices.size() ;
	    Vertex ** allVertices = new Vertex * [numVertices] ;
	
	    for (ULONG i = 0; i < numVertices; i++)
	    {
		    double x = vertices[i].x ;
		    double y = vertices[i].y ;
		
		    allVertices[i] = new Vertex(x,y) ;
	    }
	
	    return allVertices ;
    }
    
		Edge ** GenerateEdges(vector<edge> &edges, vector< vector<double> > &weights)
		{
      numEdges = (ULONG)edges.size() ;
	    Edge ** allEdges = new Edge * [numEdges] ;
	
	    for (ULONG i = 0; i < numEdges; i++)
	    {
		    Vertex * v1 = itsVertices[(ULONG)edges[i].first] ;
		    Vertex * v2 = itsVertices[(ULONG)edges[i].second] ;
		    double cost = weights[i][0] ;
		    double var = weights[i][1] ;
		
		    allEdges[i] = new Edge(v1, v2, cost, var) ;
	    }
	
	    return allEdges ;
    }
} ;

// Node class to maintain path information up to a vertex
// contains parent node, (mu, sigma) of cost-to-come
class Node
{
	public:
		Node(Vertex * vertex):
		  itsVertex(vertex), itsParent(0), itsMeanCost(0.0), itsVarCost(0.0), itsDepth(0),
		  itsHeuristic(0.0), itsMeanCTG(0.0), itsVarCTG(0.0) {}
		
		Node(Vertex * vertex, nodeType n):
		itsVertex(vertex), itsParent(0), itsHeuristic(0.0)
		{
	    switch (n)
	    {
		    case SOURCE:
			    itsMeanCost = 0.0 ;
			    itsVarCost = 0.0 ;
			    itsDepth = 0.0 ;
			    itsMeanCTG = 0.0 ;
			    break ;
		    default:
			    itsMeanCost = DBL_MAX ;
			    itsVarCost = DBL_MAX ;
			    itsDepth = ULONG_MAX ;
			    itsMeanCTG = DBL_MAX ;
			    itsVarCTG = DBL_MAX ;
	    }
    }
    
		Node(Node * parent, Edge * edge):
		itsParent(parent), itsHeuristic(0.0), itsMeanCTG(0.0), itsVarCTG(0.0)
		{
	    itsVertex = edge->GetVertex2() ;
	    itsMeanCost = itsParent->GetMeanCost() + edge->GetMeanSearch() ;
	    itsVarCost = itsParent->GetVarCost() + edge->GetVarSearch() ;
	    itsDepth = itsParent->GetDepth() + 1 ;
    }
    
		~Node(){} ;
		
		Node * GetParent() const {return itsParent ;}
		void SetParent(Node * parent) {itsParent = parent ;}
		double GetMeanCost() const {return itsMeanCost ;}
		void SetMeanCost(double cost) {itsMeanCost = cost ;}
		double GetVarCost() const {return itsVarCost ;}
		void SetVarCost(double var) {itsVarCost = var ;}
		Vertex * GetVertex() const {return itsVertex ;}
		void SetVertex(Vertex * vertex) {itsVertex = vertex ;}
		ULONG GetDepth() const {return itsDepth ;}
		void SetDepth(ULONG depth) {itsDepth = depth ;}
		double GetHeuristic() const {return itsHeuristic ;}
		void SetHeuristic(double h) {itsHeuristic = h ;}
		double GetMeanCTG() const {return itsMeanCTG ;}
		double GetVarCTG() const {return itsVarCTG ;}
		
		void DisplayPath()
		{
	    cout << "Vertex: (" << itsVertex->GetX() << "," << itsVertex->GetY() << ")\n" ;
	    cout << "Mean cost-to-come: " << itsMeanCost << ", " << "variance: " << itsVarCost << endl ;
	    cout << "Mean cost-to-go: " << itsMeanCTG << ", " << "variance: " << itsVarCTG << endl ;
	
	    if (itsParent)
		    itsParent->DisplayPath() ;
    }
    
		Node * ReverseList(Node * itsChild)
		{
	    Node * itsParentR = new Node(GetVertex()) ;
	    itsParentR->SetMeanCost(GetMeanCost()) ;
	    itsParentR->SetVarCost(GetVarCost()) ;
	    itsParentR->SetParent(itsChild) ;
	
	    Node * itsReverse ;
	
	    if (GetParent())
	    {
		    itsReverse = GetParent()->ReverseList(itsParentR) ;
		    return itsReverse ;
	    }
	    else
		    return itsParentR ;
    }

		void SetCTG(double totalMean, double totalVar)
		{
	    itsMeanCTG = totalMean - itsMeanCost ;
	    itsVarCTG = totalVar - itsVarCost ;
	
	    if (itsParent)
		    itsParent->SetCTG(totalMean, totalVar) ;
    }
    
	private:
		Vertex * itsVertex ;
		Node * itsParent ;
		double itsMeanCost ;
		double itsVarCost ;
		ULONG itsDepth ;
		double itsHeuristic ;
		double itsMeanCTG ;
		double itsVarCTG ;
} ;

// Node comparison class
class CompareNode
{
	public:
		bool operator() (const Node * n1, const Node * n2) const
		{
//			switch (SEARCH_TYPE)
//				{
//					case BREADTH:
//						return (n1->GetDepth() > n2->GetDepth()) ;
//					default:
						double n1Cost = n1->GetMeanCost() + n1->GetHeuristic() ;
						double n2Cost = n2->GetMeanCost() + n2->GetHeuristic() ;
						return (n1Cost >= n2Cost && n1->GetVarCost() >= n2->GetVarCost()) ;
//						if (n1Cost > n2Cost && n1->GetVarCost() > n2->GetVarCost())
//							return true ;
//						else if (n2Cost > n1Cost && n2->GetVarCost() > n1->GetVarCost())
//							return false ;
//						else
//							return (n1Cost > n2Cost) ;
//				}
		}
} ;

// Custom queue type to perform priority queue updates
class Queue
{
	public:
  	typedef priority_queue<Node *, vector<Node *>, CompareNode> QUEUE ;
		Queue(Node * source)
		{
	    itsPQ = new QUEUE ;
	    itsPQ->push(source) ;
    }
    
		~Queue()
		{
	    delete itsPQ ;
	    itsPQ = 0 ;
	    for (ULONG i = 0; i < closed.size(); i ++)
	    {
		    delete closed[i] ;
		    closed[i] = 0 ;
	    }
    }
		
		vector<Node *> GetClosed() const {return closed ;}
		bool EmptyQueue() const {return itsPQ->empty() ;}
		ULONG SizeQueue() const {return (ULONG)itsPQ->size() ;}
		
		void UpdateQueue(Node * newNode)
		{
	    // Compare newNode to nodes in closed set
	    // if closed contains node with same vertex, compare their costs
	    // choose whether or not to create a new node
	    bool dom = false ;
	    for (ULONG i = 0; i < closed.size(); i++)
	    {
		    if (closed[i]->GetVertex() == newNode->GetVertex())
		    {
			    dom = CompareNodes(newNode, closed[i]) ;
			    if (dom)
				    return ;
		    }
	    }
	    itsPQ->push(newNode) ;
    }
    
		Node * PopQueue()
		{
      // Check if next node is already dominated by existing node in closed set
      Node * newNode = itsPQ->top() ;
      bool dom = false ;
      for (ULONG i = 0; i < closed.size(); i++)
      {
	      if (closed[i]->GetVertex() == newNode->GetVertex())
	      {
		      dom = CompareNodes(newNode, closed[i]) ;
		      if (dom)
		      {
			      itsPQ->pop() ;
			      return 0 ;
		      }
	      }
      }
      closed.push_back(itsPQ->top()) ;
      itsPQ->pop() ;
      return closed[(ULONG)closed.size()-1] ;
    }
    
	private:
		QUEUE * itsPQ ;
		vector<Node *> closed ;
		
		bool CompareNodes(const Node * n1, const Node * n2) const
		{
	    double n1Cost = n1->GetMeanCost() ;
	    double n2Cost = n2->GetMeanCost() ;
	    return (n1Cost >= n2Cost && n1->GetVarCost() >= n2->GetVarCost()) ;
      /*	if (n1Cost > n2Cost && n1->GetVarCost() > n2->GetVarCost())
		    return true ;
	    else if (n2Cost > n1Cost && n2->GetVarCost() > n1->GetVarCost())
		    return false ;
	    else
		    return (n1Cost > n2Cost) ;*/
    }
} ;

// Path search class to store and search a graph
// A* search with zero, Manhattan or Euclidean distance heuristics
class Search
{
	public:
		Search(Graph * graph, Vertex * source, Vertex * goal):
		  itsGraph(graph), itsSource(source), itsGoal(goal) {}

		~Search()
		{
	    delete itsQueue ;
	    itsQueue = 0 ;
    }
		
		Graph * GetGraph() const {return itsGraph ;}
		Queue * GetQueue() const {return itsQueue ;}
		void SetQueue(Queue * queue) {itsQueue = queue ;}
		Vertex * GetSource() const {return itsSource ;}
		Vertex * GetGoal() const {return itsGoal ;}
		
		vector<Node *> PathSearch(pathOut pType)
		{
	    ULONG sourceID = FindSourceID() ;
	    itsQueue = new Queue(new Node(itsGraph->GetVertices()[sourceID], SOURCE)) ;
	
	    clock_t t_start = clock() ;
	    double t_elapse = 0.0 ;
	
	    while (!itsQueue->EmptyQueue() && t_elapse < 5)
	    {
		    // Pop cheapest node from queue
		    Node * currentNode = itsQueue->PopQueue() ;
		    if (!currentNode) { continue ; }
		
		    if (pType == BEST)
		    {
			    // Terminate search once one path is found
			    if (currentNode->GetVertex() == itsGoal)
				    break ;
		    }
		
		    Node * currentNeighbour ;
		
		    // Find all neighbours
		    vector<Edge *> neighbours = itsGraph->GetNeighbours(currentNode->GetVertex()) ;

		    // Update neighbours
		    for (ULONG i = 0; i < (ULONG)neighbours.size(); i++)
		    {
			    // Check if neighbour vertex is already in closed set
			    bool newNeighbour = true ;
			    if (pType == BEST)
			    {
				    Vertex * vcheck = neighbours[i]->GetVertex2() ;
				    for (ULONG j = 0; j < itsQueue->GetClosed().size(); j++)
				    {
					    if (itsQueue->GetClosed()[j]->GetVertex() == vcheck)
					    {
						    newNeighbour = false ;
						    break ;
					    }
				    }
			    }
			
			    if (newNeighbour)
			    {
				    // Create neighbour node
				    currentNeighbour = new Node(currentNode, neighbours[i]) ;
				    UpdateNode(currentNeighbour) ;
				    itsQueue->UpdateQueue(currentNeighbour) ;
			    }
		    }
		
		    t_elapse = (float)(clock() - t_start)/CLOCKS_PER_SEC ;
	    }
	
	    // Check if a path is found
	    bool ClosedAll = false ;
	    for (ULONG i = 0; i < itsQueue->GetClosed().size(); i++)
	    {
		    if (itsQueue->GetClosed()[i]->GetVertex() == itsGoal)
		    {
			    ClosedAll = true ;
			    break ;
		    }
	    }
	
	    if (!ClosedAll)
	    {
		    cout << "No path found from source to goal.\n" ;
		    vector<Node *> bestPath ;
		    return bestPath ;
	    }
	    else
	    {
		    ULONG k = 0 ;
		    vector<Node *> bestPath((ULONG)itsQueue->GetClosed().size()) ;
		
		    for (ULONG i = 0; i < (ULONG)itsQueue->GetClosed().size(); i++)
		    {
			    if (itsGoal == itsQueue->GetClosed()[i]->GetVertex())
			    {
				    bestPath[k] = itsQueue->GetClosed()[i] ;
				    k++ ;
			    }
		    }
		
		    bestPath.resize(k) ;
		    return bestPath ;
	    }
    }

	private:
		Graph * itsGraph ;
		Queue * itsQueue ;
		Vertex * itsSource ;
		Vertex * itsGoal ;
		
		ULONG FindSourceID()
		{
	    for (ULONG i = 0; i < itsGraph->GetNumVertices(); i++)
	      if (itsSource == itsGraph->GetVertices()[i])
		      return i ;
    }
    
		double ManhattanDistance(Vertex * v1, Vertex * v2)
		{
	    double diffX = abs(v1->GetX() - v2->GetX()) ;
	    double diffY = abs(v1->GetY() - v2->GetY()) ;
	    double diff = diffX + diffY ;
	    return diff ;
    }
    
		double EuclideanDistance(Vertex * v1, Vertex *v2)
		{
	    double diffX = pow(v1->GetX() - v2->GetX(),2) ;
	    double diffY = pow(v1->GetY() - v2->GetY(),2) ;
	    double diff = sqrt(diffX+diffY) ;
	    return diff ;
    }
    
		void UpdateNode(Node * n)
		{
	    if (SEARCH_TYPE == ASTAR)
	    {
		    double diff ;
		    switch (HEURISTIC)
		    {
			    case ZERO:
				    diff = 0.0 ;
				    n->SetHeuristic(diff) ;
				    break ;
			    case MANHATTAN:
				    diff = ManhattanDistance(itsGoal, n->GetVertex()) ;
				    n->SetHeuristic(diff) ;
				    break ;
			    case EUCLIDEAN:
				    diff = EuclideanDistance(itsGoal, n->GetVertex()) ;
				    n->SetHeuristic(diff) ;
				    break ;
		    }
	    }
	    else if (SEARCH_TYPE == DIJKSTRA)
	    {
	      HEURISTIC = ZERO ;
	      SEARCH_TYPE = ASTAR ;
	      UpdateNode(n) ;
	    }
    }
} ;
    
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

// Function to compare vertices: returns TRUE if vertex A is better than vertex B
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

  // If code is taking too long, change n to a smaller value
  // Note that since this is the numerical integration discretisation
  // smaller n will provide coarser approximation to solution
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

// RAGS class for interfacing with sector agents
class RAGS
{
  public:
    typedef pair<int, int> edge ;
    RAGS(vector<XY> &locations, vector<edge> &edge_array, vector< vector<double> > &weights)
    { itsGraph = new Graph(locations, edge_array, weights) ; }
    
    ~RAGS()
    {
      delete itsGraph ;
      itsGraph = 0 ;
      delete itsSearch ;
      itsSearch = 0 ;
    }
    
    Graph * GetGraph() const {return itsGraph ;}
    Vertex * GetVert() const {return itsVert ;}
    
    void SetInitialVert(Vertex * start)
    {
      Vertex ** allVerts = itsGraph->GetVertices() ;
      for (ULONG i = 0; i < itsGraph->GetNumVertices(); i++)
      {
        if ( start == allVerts[i] )
          itsVert = allVerts[i] ;
      }
    }
    
    XY SearchGraph(XY start, XY goal, vector<double> &weights)
    // Create function that converts from XY to Vertex *
    {
      Vertex ** allVerts = itsGraph->GetVertices() ;
      Vertex * sVert ;
      Vertex * gVert ;
      for (ULONG i = 0; i < itsGraph->GetNumVertices(); i++)
      {
        if ( start.x == allVerts[i]->GetX() && start.y == allVerts[i]->GetY() )
          sVert = allVerts[i] ;
      }
      for (ULONG i = 0; i < itsGraph->GetNumVertices(); i++)
      {
        if ( goal.x == allVerts[i]->GetX() && goal.y == allVerts[i]->GetY() )
          gVert = allVerts[i] ;
      }
      return SearchGraph(sVert, gVert, weights) ;
    }
    
    XY SearchGraph(Vertex * start, Vertex * goal, vector<double> &weights)
    {
      AssignCurrentEdgeCosts(weights) ;
      
      // Initialise non-dominated path set
      if (itsNDSet.empty())
      {
        itsSearch = new Search(itsGraph, start, goal) ;
        vector<Node *> GSPaths = itsSearch->PathSearch(PSET) ;
        
        if (GSPaths.empty())
        {
          XY s = XY(start->GetX(),start->GetY()) ;
	        return s ; // no paths found, stay where you are
        }
	      
	      for (ULONG i = 0; i < (ULONG)GSPaths.size(); i++)
		      itsNDSet.push_back(GSPaths[i]->ReverseList(0)) ;
	      for (ULONG i = 0; i < (ULONG)itsNDSet.size(); i++)
	        itsNDSet[i]->SetCTG(GSPaths[i]->GetMeanCost(),GSPaths[i]->GetVarCost()) ;
        SetInitialVert(itsNDSet[0]->GetVertex()) ;
      }
	      
      // Flag and return if current vertex does not match start vertex
      if (!(itsVert == start))
      {
        cout << "\nError: input start vertex (" << start->GetX() << "," << start->GetY() <<
        ") does not match stored vertex (" << itsVert->GetX() << "," << itsVert->GetY() << ")." <<
        endl ; 
        XY s = XY(start->GetX(),start->GetY()) ;
        return s ;
      }
      
      vector<Node *> newNodes = itsNDSet ;
      itsVert->SetNodes(newNodes) ;
      vector<Node *> tmpNodes ;
      vector<Vertex *> nextVerts ;
      
      // Identify next vertices
		  for (int i = 0; i < newNodes.size(); i++)
		  {
			  bool newVert = true ;
			  for (int j = 0; j < nextVerts.size(); j++)
			  {
				  if (nextVerts[j] == newNodes[i]->GetParent()->GetVertex() || nextVerts[j] == itsVert)
				  {
					  newVert = false ;
					  break ;
				  }
			  }
			  if (newVert)
				  nextVerts.push_back(newNodes[i]->GetParent()->GetVertex()) ;
		  }
		  
		  // Identify next vertex path nodes
		  for (int i = 0; i < nextVerts.size(); i++)
		  {
			  tmpNodes.clear() ;
			  for (int j = 0; j < newNodes.size(); j++)
			  {
				  if (nextVerts[i] == newNodes[j]->GetParent()->GetVertex())
				  {
					  tmpNodes.push_back(newNodes[j]->GetParent()) ;
				  }
			  }
			  nextVerts[i]->SetNodes(tmpNodes) ;
		  }
		  
		  // Rank next vertices according to probability of improvement
  		sort(nextVerts.begin(),nextVerts.end(),ComputeImprovementProbability) ;
  		
  		itsVert = nextVerts[0] ;
  		itsNDSet = itsVert->GetNodes() ;
  		XY vertXY = XY(itsVert->GetX(),itsVert->GetY()) ;
  		return vertXY ;
    }
    
  private:
    Graph * itsGraph ;
    Search * itsSearch ;
    Vertex * itsVert ;
    vector<Node *> itsNDSet ;
    
    void AssignCurrentEdgeCosts(vector<double> &weights)
    {
      ULONG n = itsGraph->GetNumEdges() ;
      Edge ** e = itsGraph->GetEdges() ;
      
      for (ULONG i = 0; i < n; i++)
        e[i]->SetTrueCost(weights[i]) ;
    }
} ;
