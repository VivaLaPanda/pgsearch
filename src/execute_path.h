#include "vector"
//#include "node.h"

void executePath(vector<Node*> bestPaths)
{
	cout<<"PARENT" << bestPaths[0]->GetParent() << endl;
	for(int i =0; i < bestPaths.size(); i++){
		cout << "PARENT"<< i << bestPaths[i]->GetParent()<< endl;
		//GetVertex();
		//stuff
	}
}

