class CompareNode
{
	public:
		bool operator() (const Node * n1, const Node * n2) const
		{
			switch (SEARCH_TYPE)
				{
					case BREADTH:
						return (n1->GetDepth() > n2->GetDepth()) ;
					default:
						double n1Cost = n1->GetMeanCost() + n1->GetHeuristic() ;
						double n2Cost = n2->GetMeanCost() + n2->GetHeuristic() ;
						return (n1Cost >= n2Cost && n1->GetVarCost() >= n2->GetVarCost()) ;
						/*if (n1Cost > n2Cost && n1->GetVarCost() > n2->GetVarCost())
							return true ;
						else if (n2Cost > n1Cost && n2->GetVarCost() > n1->GetVarCost())
							return false ;
						else
							return (n1Cost > n2Cost) ;*/
				}
		}
} ;
