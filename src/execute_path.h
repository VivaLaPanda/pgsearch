#ifndef EXECUTE_PATH
#define EXECUTE_PATH

#include <vector>
//#include "node.h"
void executePath(vector<Node*> bestPaths)
{
	vector< Node *> reversedPaths;
	cout << endl << "EXECUTING PATH" << endl;
	for(int i = 0; i < bestPaths.size(); i++){
		reversedPaths[i] = bestPaths[i]->ReverseList(bestPaths[i]);
		cout << endl << "Path" << i << endl;
		reversedPaths[i]->DisplayPath();
		cout << endl;
		//bestPaths[i]->GetVertex()->SetCTG(10);
	}

	for(int i = 0; i < reversedPaths.size(); i++){
		cout << endl << "Path" << i << endl;
		reversedPaths[i]->DisplayPath();  
		cout << endl;
		// find the mean cost of the first & last node
		// 
		reversedPaths[i]->GetMeanCost();
		Node * itsParent;
		Node * oldParent;
		itsParent = reversedPaths[i];
		while(itsParent != oldParent)
			oldParent = itsParent;
			itsParent = itsParent->GetParent();
		cout << "mean cost = " << reversedPaths[i]->GetMeanCost() - itsParent->GetMeanCost() << endl;
		
	}

	cout << endl << reversedPaths[2]->GetVertex()->GetCTG() << endl;

	/* Get list of paths
	 * Call *probability of best path*
	 * Take next action based on sorted (store in final path taken!)
	 i* Repeat
	 */

	/*set attributes of each vertice with score from them to the goal
	*/
}


#endif
