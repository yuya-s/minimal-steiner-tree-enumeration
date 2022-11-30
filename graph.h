#ifndef FRONTIER_NETWORK_GRAPH_H
#define FRONTIER_NETWORK_GRAPH_H

#include "utility.h"


class Index;

struct Edge {
public:
	int src;  // edge source
	int dest; // edge destination
    double weight;

	Edge(int src_a, int dest_a, double weight_a)
	{
		src = src_a;
		dest = dest_a;
        weight = weight_a;
	}
};


//******************************************************************************
// Original structure that are inputted by users
// undirected graphs。
class OriginalGraph {

public:
	int number_of_vertices; // the number of vertices
    int number_of_edges;
	int number_of_terminals;
	double correct_answer;
	vector<Edge> edge_list; // edge list
    vector<int> degrees;
    vector<int> terminals;
    vector<vector<int>> vertex_neighbor_edge;

	int max_degree;

	OriginalGraph(){}

	OriginalGraph(int _number_of_vertices){
		number_of_vertices = _number_of_vertices;
        //number_of_edges = 0;
        degrees.resize(number_of_vertices);
        vertex_neighbor_edge.resize(number_of_vertices);
		for(int i=0;i<number_of_vertices;i++){
            degrees[i]=0;
            vertex_neighbor_edge[i].clear();
        }
		max_degree=0;
		edge_list.clear();
        terminals.clear();
	}

	int GetNumberOfVertices() // return the number of vertices
	{
		return number_of_vertices;
	}

	vector<Edge>* GetEdgeList() // return edge list
	{
		return &edge_list;
	}

    void InputFile(std::string);
	std::vector<OriginalGraph> GraphPrune(Index *_index);
	void GraphReduction();
};

class SeedSteinerTreeSet : public OriginalGraph
{
public:
    int number_of_seedtrees;
    vector<double> seedtrees_answer;

    void InputFile(OriginalGraph* , vector<string>);

};

//******************************************************************************
// BDDinputGraph structure
// undirected graphs for inputing BDD
class BDDinputGraph {

public:
	int number_of_vertices; // The number of vertices
	vector<Edge> edge_list; // edge list
	vector<int> degrees;
	vector<bool> terminals;
	int max_degree;

	BDDinputGraph(){
	}

	vector<Edge>* GetEdgeList() // return edge list
	{
		return &edge_list;
	}

};

#endif //FRONTIER_NETWORK_GRAPH_H
