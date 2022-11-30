//
// Created by sasaki on 17/09/27.
//

#ifndef FRONTIER_NETWORK_FRONTIER_H
#define FRONTIER_NETWORK_FRONTIER_H

#include "bdd.h"
#include "utility.h"
#include "state.h"
#include "graph.h"
#include <bits/stdc++.h>
using namespace std;


class FrontierAlgorithm {

public:

    //Compute Probability and Delete upper part of BDD for memory efficiency
    BDD* ComputeSteinertree(State *state, unsigned long int max_elements, double _error, double weight);
    void MemoryTest(State* state, int max_elements);

private:

	FrontierNode* CheckTerminal(FrontierNode* n_hat, int i, int x, State* state, BDDNode& bddnode, double weight_);


	//void UpdateInfo(FrontierNode* n_hat, int i, int x, State* state);


	FrontierNode* UpdateInfo(FrontierNode *n_hat, int i, int x, State *state, BDDNode& bddnode, double weight_);


	FrontierNode* Find(FrontierNode* n_prime, const vector<FrontierNode*>& N_i, int i, State* state);

    void HashFrontier(FrontierNode* n_prime);

	bool Sampling(FrontierNode* n_prime, int i, State* state, mt19937 mt);

	//bool CheckConnectivity(vector< vector <int> > *neighborList, State* state);

	bool IsEquivalent(FrontierNode* node1, FrontierNode* node2, int i, State* state);

};

#endif //FRONTIER_NETWORK_FRONTIER_H
