// Node class to maintain path information up to a vertex
// contains parent node, (mu, sigma) of cost-to-come

enum nodeType {SOURCE, OTHER} ;

class Node
{
	public:
		Node(Vertex *) ;
		Node(Vertex *, nodeType) ;
		Node(Node *, Edge *) ;
		~Node() ;
		
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
		
		void DisplayPath() ;
	private:
		Vertex * itsVertex ;
		Node * itsParent ;
		double itsMeanCost ;
		double itsVarCost ;
		ULONG itsDepth ;
		double itsHeuristic ;
} ;

Node::Node(Vertex * vertex)
{
	itsVertex = vertex ;
	itsParent = 0 ;
	itsMeanCost = 0.0 ;
	itsVarCost = 0.0 ;
	itsDepth = 0 ;
	itsHeuristic = 0.0 ;
}

Node::Node(Vertex * vertex, nodeType n)
{
	switch (n)
	{
		case SOURCE:
			itsVertex = vertex ;
			itsParent = 0 ;
			itsMeanCost = 0.0 ;
			itsVarCost = 0.0 ;
			itsDepth = 0.0 ;
			itsHeuristic = 0.0 ;
			break ;
		default:
			itsVertex = vertex ;
			itsParent = 0 ;
			itsMeanCost = DBL_MAX ;
			itsVarCost = DBL_MAX ;
			itsDepth = ULONG_MAX ;
			itsHeuristic = 0.0 ;
	}
}

Node::Node(Node * parent, Edge * edge)
{
	itsVertex = edge->GetVertex2() ;
	itsParent = parent ;
	itsMeanCost = itsParent->GetMeanCost() + edge->GetMeanCost() ;
	itsVarCost = itsParent->GetVarCost() + edge->GetVarCost() ;
	itsDepth = itsParent->GetDepth() + 1 ;
	itsHeuristic = 0.0 ;
}

Node::~Node()
{
}

void Node::DisplayPath()
{
	cout << "Vertex: (" << itsVertex->GetX() << "," << itsVertex->GetY() << ")\n" ;
	cout << "Mean cost: " << itsMeanCost << endl ;
	cout << "Variance: " << itsVarCost << endl ;
	
	if (itsParent)
		itsParent->DisplayPath() ;
}
