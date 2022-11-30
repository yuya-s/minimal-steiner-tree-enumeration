#ifndef FRONTIER_NETWORK_ZDD_H
#define FRONTIER_NETWORK_ZDD_H

#include "utility.h"
#include "graph.h"
#include "state.h"


//******************************************************************************
// FrontierNode 構造体
// Frontier計算用のノードを表す。
class FrontierNode {


public:

	vector<int> deg;  // deg 配列（フロンティア法アルゴリズムの文献参照）
	vector<int> comp; // comp 配列（フロンティア法アルゴリズムの文献参照）
    vector<int> terminal; // terminal 配列（フロンティア法アルゴリズムの文献参照）
	vector<int> addedge;
    double hueristicValue;
    double number_of_existence;
    int nodeId;
	bool add_non_terminal;
    double weight;
	// BDD, FrontierNode を使用するときはこの関数を必ず呼ぶこと。
	// 終端ノードの初期化を行う
    //。
	static void Initialize()
	{

	}

	~FrontierNode()
	{

	}
//    bool operator<(const FrontierNode* another)
//    {
//        return hueristicValue < another->hueristicValue ;
//    }



//    bool operator>(const FrontierNode* right)
//    {
//    return hueristicValue < right->hueristicValue ;
//    }



    void HueristicValue_terminal(){

        hueristicValue=0;
        for(int i =0;i<terminal.size();i++){
            hueristicValue+=terminal[i];
        }
        hueristicValue-=weight*0.001;
        hueristicValue=number_of_existence;
    }

    void HueristicValue_num_of_existence(){

        hueristicValue=number_of_existence;

    }



	static FrontierNode* CreateRootNode(State *state)
	{

		FrontierNode* node = new FrontierNode();
		//node->SetNextId();

        //int numberVertices = state->graph->GetNumberOfVertices();

		node->deg.resize(0);// = new int[0];
		node->comp.resize(0);// = new int[0];
        node->terminal.resize(0);// =  new int[0];
		node->addedge.resize(0);
        node->hueristicValue=0;
        node->number_of_existence=0;
        node->weight=0;


		//node->hueristicValue=1;

		// deg, comp を初期化

		return node;
	}


	FrontierNode* MakeCopy()
	{
		FrontierNode* node = new FrontierNode();
		node->deg.resize(deg.size());// = new int[0];
		node->comp.resize(deg.size());// = new int[0];
        node->terminal.resize(deg.size());// =  new int[0];
		node->addedge.resize(deg.size());

		//node->hueristicValue = hueristicValue;

		for(int i=0;i<deg.size();i++){
			node->deg[i]=deg[i];
			node->comp[i]=comp[i];
			node->terminal[i]=terminal[i];
			node->addedge[i]=addedge[i];
		}
		return node;
	}



    int RelaxEditDistance(FrontierNode* _zdd, int _threshold){

        int distance=0;
        for(unsigned int i=0; i < _zdd->comp.size();i++){

			if(terminal[i] > 0 && comp[i] != _zdd->comp[i])return _threshold;

			//if(comp[i] < _zdd->comp[i])return _threshold;
			if(comp[i] != _zdd->comp[i])distance++;

			if(distance >= _threshold)return _threshold;
        }
        return distance;
    }



//    bool operator>(const FrontierNode &zdd_node) const
//    {
//		return existenceProbability > zdd_node.existenceProbability;
//    }
//
//	bool operator<(const FrontierNode &zdd_node) const
//    {
//		return !(this->operator>(zdd_node));
//		//existenceProbability < zdd_node.existenceProbability;
//    }


	//bool operator()(const FrontierNode *zdd1, const FrontierNode *zdd2) const{
	//	return zdd1->hueristicValue < zdd2->hueristicValue;
	//}


};

static int BDDSize(int frontierSize){

	return sizeof(int)*frontierSize*3*2;
}





//******************************************************************************
// FrontierNode
class BDDNode {

public:

	int existencePointer; // -1 is to true, -2 is to false
	int nonexistencePointer; // -1 is to true, -2 is to false


	static void Initialize()
	{

	}

	BDDNode()
	{
		existencePointer=-2;
		nonexistencePointer=-2;
	}


	~BDDNode()
	{

	}

};




//******************************************************************************
// BDD 構造体
// BDD のノードが node_list_array に格納される。
// レベル i のノードは (*node_list_array)[i] に格納される。
// i は1始まり。0は使わない。i = m + 1 はダミー。
// レベル i の j 番目のノードは (*node_list_array)[i][j] で参照できる。
class BDD
{
public:
	vector<vector<BDDNode>> node_list_array;
	double minimumCost;
	double averageCost;
	int treenum;

	BDD()
	{
		node_list_array.clear();
	}

	void Initialize(int n){
		node_list_array.resize(n);
	}

	~BDD()
	{
		for (unsigned int i = 0; i < node_list_array.size(); ++i) {
//			for (unsigned int j = 0; j < node_list_array[i].size(); ++j) {
				node_list_array[i].clear();
//			}
		}
//		delete node_list_array;
	}
    void Reduce();
	OriginalGraph Traversal(State* state,double threshold, int);
    OriginalGraph MinimumTraversal(State* state);
	int GetSize();
};


class BDDNodeonPriorityQueue{
public:
    int edgeIdinBDD; // edge ID;
    int nodeId; // node ID in BDD;
	double cost;

	BDDNodeonPriorityQueue(){

	};

	BDDNodeonPriorityQueue(int arrayId, int nodeId, double cost){ // コンストラクタ
        this->edgeIdinBDD = arrayId;
        this->nodeId = nodeId;
		this->cost = cost;
    };
	bool operator<(const BDDNodeonPriorityQueue &b) const{
	    return cost < b.cost;
	};
};

class BDDNodeTraversal{
public:
    //int edgeIdinBDD; // edge ID;
    int nodeId; // node ID in BDD;
//	vector<double> cost;

	BDDNodeTraversal(){

	};

	BDDNodeTraversal(int nodeId){ // コンストラクタ
        //this->edgeIdinBDD = arrayId;
        this->nodeId = nodeId;
//		this->cost = cost;
    };

};


#endif //FRONTIER_NETWORK_ZDD_H
