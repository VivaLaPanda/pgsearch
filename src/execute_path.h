#include "vector"

void executePath(vector<Node*> bestPaths)
{
	//cout<<"PARENT" << bestPaths[0]->GetParent() << endl;
	//vector< Vertex *> Vertices;
	//for(int i =0; i < bestPaths.size(); i++){
	//	cout << "PARENT"<< i << bestPaths[i]->GetParent()->GetVertex()<< endl;
	//	Vertices.push_back(bestPaths[i]->GetParent()->GetVertex());
	//	cout << "VERTICES " <<Vertices[i]<<endl;
	//	Node* p;
	//	Node* cur;
	//	cur = bestPaths[i];
	//	while(cur){
	//		p = cur->GetParent();
	//		cout << p->GetVertex() << endl;
	//		cur = p;
	//	}

	//}
	//cout << endl << bestPaths[1] << endl;

	for(int i = 0; i < bestPaths.size(); i++){
		cout << endl << bestPaths[i] << endl;
		bestPaths[i]->GetVertex()->SetCTG(10);
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



