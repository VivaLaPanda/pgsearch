#ifndef EXECUTE_PATH
#define EXECUTE_PATH

#include <vector>
//#include "node.h"
using namespace std;

void executePath(vector< Node*> GSPaths){
	cout << endl << "EXECUTING PATH" << endl;
	Vertex* goal_loc;
	Node* goal;
	Node* cur_loc;
	vector < Node*> SGPaths;

	for(int i = 0; i < GSPaths.size(); i++){
		SGPaths.push_back(GSPaths[i]->ReverseList(0));
		cout << endl << "Path" << i << endl;
		SGPaths[i]->DisplayPath();
		cout << endl;
	}

	cout << "Iterating Through Paths" << endl;
	cur_loc = SGPaths[0];
	for(int i = 0; i < SGPaths.size(); i++){

		cout << "Path" <<  i  << endl;
		goal = GSPaths[i];
		goal_loc = goal->GetVertex();
		cur_loc = SGPaths[i];
		double gu = goal->GetMeanCost();
		double gv = goal->GetVarCost();

		cout << "Before While" << endl;
		while(cur_loc->GetVertex() != goal_loc){
			cout << "Step" << endl;
			cur_loc = cur_loc->GetParent();
			cout << "Current Location: " << cur_loc->GetVertex()->GetX() << ", " <<  cur_loc->GetVertex()->GetY() << endl;
			//Track the score
			cout << gu - cur_loc->GetMeanCost() << endl;
			cout << gv - cur_loc->GetVarCost() << endl << endl;
			

		}
	}

}

/*
void executePath(vector<Node*> bestPaths)
{
	vector< Node *> reversedPaths(bestPaths.size());
	vector< Node *> vertices;
	cout << endl << "EXECUTING PATH" << endl;
	Vertex * cur_loc;
	vector< Node *> parents;

	while(cur_loc != bestPaths.end()){
		parents = []
		for(int i = 0; i < bestPaths.size(); i++){
			reversedPaths[i] =  bestPaths[i]->ReverseList(0);
			cout << endl << "Path" << i << endl;
			reversedPaths[i]->DisplayPath();
			cout << endl;
		}

		for(int i = 0; i < reversedPaths.size(); i++){
			cout << endl << "Path" << i << endl;
			reversedPaths[i]->DisplayPath();
			cout << endl;
			// store parents for next iter
			parents[i] = reversedPaths[i]->GetParent();// this is wrong and needs to be dependent on next loc vertices

			// find the mean cost of the first & last node

			reversedPaths[i]->GetMeanCost();
			bestPaths[i]->GetMeanCost();

			// Here add in while loop to iterate through each node along path

			if(find(vertices.begin(), vertices.end(), reversedPaths[i]->GetVertex())== vertices.end()){ //need to rewrite this correctly
				vertices.push_back(reveresedPaths[i]->GetVertex());
			}

		}

		//New paths/location
		reversedPaths = [];
		reversedPaths = parents;

		cur_loc = bestPaths[i]->GetVertex; // need to change this later
		cout << "CURRENT LOCATION: " << cur_loc->GetVertex();


		/* Get list of paths
		 * find all paths to get to different vertices, I think this needs to be only the next level because otherwise we will need to trim as we go becaues if we take a round about path to a vertex that doesn't get trimmed.
		 * create a vector of paths to vertices (aka vector of nodes that are at a vertice)
		 * Call *probability of best path*
		 * Take next action based on sorted (store in final path taken!)
		 i* Repeat
		 */

		/*
		 * set attributes of each vertice with score from them to the goal
		
	}
}

*/

#endif
