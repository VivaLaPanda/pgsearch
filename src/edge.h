// Edge class to contain mean and variance of cost along an edge
#include <random>
using namespace std;
class Edge
{
	public:
		Edge(Vertex * v1, Vertex * v2, double cost, double var) ;
		~Edge() ;
		
		Vertex * GetVertex1() const {return itsVertex1 ;}
		void SetVertex1(Vertex * vertex1) {itsVertex1 = vertex1 ;}
		Vertex * GetVertex2() const {return itsVertex2 ;}
		void SetVertex2(Vertex * vertex2) {itsVertex2 = vertex2 ;}
		double GetMeanCost() const {return itsMeanCost ;}
		void SetMeanCost(double cost) {itsMeanCost = cost ;}
		double GetVarCost() const {return itsVarCost ;}
		void SetVarCost(double var) {itsVarCost = var ;}
		double GetTrueCost() {return itsTrueCost ;}
		void SetTrueCost() ;
	private:
		Vertex * itsVertex1 ;
		Vertex * itsVertex2 ;
		double itsMeanCost ;
		double itsVarCost ;
		double itsTrueCost ;
} ;

Edge::Edge(Vertex * v1, Vertex * v2, double cost, double var)
{
	itsVertex1 = v1 ;
	itsVertex2 = v2 ;
	itsMeanCost = cost ;
	itsVarCost = var ;
}

Edge::~Edge()
{
	delete itsVertex1 ;
	itsVertex1 = 0 ;
	delete itsVertex2 ;
	itsVertex2 = 0 ;
}

void Edge::SetTrueCost()
{
	default_random_engine generator ;
	normal_distribution<double> distribution(itsMeanCost,itsVarCost) ;
	itsTrueCost = distribution(generator) ;
	if (itsTrueCost < 0.0)
		itsTrueCost = 0.0 ; // May need to change this to Euclidean distance between vertices
}
