

#ifndef FRONTIER_NETWORK_STATE_H
#define FRONTIER_NETWORK_STATE_H

#include "utility.h"
#include "graph.h"

//******************************************************************************
// State structure

class State
{
public:
	BDDinputGraph* graph; // graph
    vector<int> terminals;
	vector<int> Fc; // current Frontier
	vector<int> Fp; // previous Frontier
	vector<int> **F;

public:
	State(){};

	State(SeedSteinerTreeSet* g, vector<int> terminals_,string ordering_)
	{
        terminals = terminals_;
    	if(ordering_=="BFS")graph = CopyGraphBFSbased(g, terminals);
		else cerr;

		Fc.clear();
		Fp.clear();
	}


	~State()
	{

        delete graph;
	}

	void SetState(SeedSteinerTreeSet* g, vector<int> terminals_,string ordering_)
	{
        terminals = terminals_;

		if(ordering_=="BFS")graph = CopyGraphBFSbased(g, terminals);

		else cerr;

		Fc.clear();
		Fp.clear();
	}

	void ComputeFrontier();
	void ComputeFrontier(int i);
	BDDinputGraph* CopyGraphBFSbased(SeedSteinerTreeSet* original_graph, vector<int> terminals_);
    BDDinputGraph* CopyGraphShortestbased(SeedSteinerTreeSet *original_graph, vector<int> terminals_);
    BDDinputGraph* CopyGraphKEdgebased(SeedSteinerTreeSet* original_graph, vector<int> terminals_,Index* index);
    BDDinputGraph* CopyGraph(SeedSteinerTreeSet* original_graph, vector<int> terminals_);

	bool FindElement(int edge_number, int value);



};
#endif //FRONTIER_NETWORK_STATE_H
