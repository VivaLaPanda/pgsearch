// Vertex class to contain location of vertices
class Node;
class Vertex
{
	public:
		Vertex(double x, double y) ;
		~Vertex() {}
		
		double GetX() const {return itsX ;}
		void SetX(double x) {itsX = x ;}
		double GetY() const {return itsY ;}
		void SetY(double y) {itsY = y ;}
		double GetCTG()const{return itsCTG ;}
		void SetCTG(double cost_to_goal) {itsCTG = cost_to_goal ;}
		void SetNodes(Node* NewNode) {Nodes.push_back(NewNode) ;}
		vector<Node*> GetNodes() {return Nodes ; }
		double SetActualCost(double ac) {itsActualCost = ac;}
		double GetActualCost() const {return itsActualCost ;}


	private:
		double itsActualCost;
		double itsX ;
		double itsY ;
		double itsCTG;
		vector< Node* > Nodes;
};

Vertex::Vertex(double x, double y)
{
	itsX = x ;
	itsY = y ;
}
