#include "vector"

void executePath(vector<Node*> bestPaths)
{
	cout << endl << "EXECUTING PATH" << endl;
	for(int i = 0; i < bestPaths.size(); i++){
		bestPaths[i]->ReverseList(bestPaths[i]);
		cout << endl << "Path" << i << endl;
		bestPaths[i]->DisplayPath();
		cout << endl;
		//bestPaths[i]->GetVertex()->SetCTG(10);
	}

	for(int i = 0; i < bestPaths.size(); i++){
		cout << endl << "Path" << i << endl;
		bestPaths[i]->DisplayPath();  
		cout << endl;
	}

	cout << endl << bestPaths[2]->GetVertex()->GetCTG() << endl;

	/* Get list of paths
	 * Call *probability of best path*
	 * Take next action based on sorted (store in final path taken!)
	 i* Repeat
	 */

	/*set attributes of each vertice with score from them to the goal
	*/
}



