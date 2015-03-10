// Vertex class to contain location of vertices

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
		double SetCTG(double cost_to_goal) {itsCTG = cost_to_goal ;}
	private:
		double itsX ;
		double itsY ;
		double itsCTG;
} ;

Vertex::Vertex(double x, double y)
{
	itsX = x ;
	itsY = y ;
}
