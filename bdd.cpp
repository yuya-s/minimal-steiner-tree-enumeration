#include "bdd.h"


bool operator>(const BDDNodeonPriorityQueue& a, const BDDNodeonPriorityQueue& b)
{
    return a.cost > b.cost;
}


void BDD::Reduce(){

    vector<vector<bool>> reducedFlags;

    reducedFlags.resize(node_list_array.size());

    for(unsigned int i = 0; i < node_list_array.size();i++){
        reducedFlags[i].resize(node_list_array[i].size());
        for(unsigned int j = 0; j < node_list_array[i].size();j++){
            reducedFlags[i][j] = false;
        }
    }


    for(int for_edge=node_list_array.size()-1; for_edge>0;for_edge--) {

        for(int for_node=0;for_node<node_list_array[for_edge].size();for_node++) {
            bool existencePointerToFalse=false;
            bool nonexistencePointerToFalse=false;

            int pointerToNextNode = node_list_array[for_edge][for_node].existencePointer;
            if(pointerToNextNode==-2){
                existencePointerToFalse=true;
            }
            else if(pointerToNextNode>-1){
                if(reducedFlags[for_edge+1][pointerToNextNode])existencePointerToFalse=true;
            }

            int non_pointerToNextNode = node_list_array[for_edge][for_node].nonexistencePointer;
            if(non_pointerToNextNode==-2){
                nonexistencePointerToFalse=true;
            }
            else if(non_pointerToNextNode>-1){
                if(reducedFlags[for_edge+1][non_pointerToNextNode])nonexistencePointerToFalse=true;
            }

            if(existencePointerToFalse&&nonexistencePointerToFalse)reducedFlags[for_edge][for_node]=true;
            //cout<<for_edge<<","<<for_node<<":"<<reducedFlags[for_edge][for_node]<<":";
            //cout<<pointerToNextNode<<","<<non_pointerToNextNode<<endl;

            int a;
        }
    }

    for(int for_edge=0; for_edge<node_list_array.size();for_edge++) {
        for(int for_node=0;for_node<node_list_array[for_edge].size();for_node++) {
            //cout<<for_edge<<","<<for_node<<":"<<reducedFlags[for_edge][for_node]<<endl;
            if(reducedFlags[for_edge][for_node]) {
                node_list_array[for_edge][for_node].existencePointer=-2;
                node_list_array[for_edge][for_node].nonexistencePointer=-2;
            }
        }
    }
}

OriginalGraph BDD::Traversal(State* state, double threshold, int  steinernum){


    stack<BDDNodeTraversal> bddstack_next;
    stack<BDDNodeTraversal> bddstack_current;

    bddstack_current.push(BDDNodeTraversal(0));

    vector<vector<double>> node_visited_cost_current;
    vector<vector<double>> node_visited_cost_next;

    node_visited_cost_current.resize(1);
    node_visited_cost_current[0].resize(1);
    node_visited_cost_current[0][0]=0;

    vector<bool> visited_node_current;
    //vector<vector<int>> previsited_node_next;

    //node_visited_cost[0][0]=0;
    BDDNodeTraversal currentTop;
    double minimumCost =  std::numeric_limits<double>::max();
    double averageCost =0;
    double edgeCost;

    //for(unsigned int i = 0; i < node_list_array.size();i++){//node_list_array.size = # of edges.
    int treeNum=0;
    int sumEdgeCost=0;

    for(int for_edge=0; for_edge<node_list_array.size();for_edge++) {

        Edge edge = (*state->graph->GetEdgeList())[for_edge];
        edgeCost = edge.weight;
        sumEdgeCost+=edgeCost;
        if(for_edge<node_list_array.size()-1) {
            node_visited_cost_next.clear();
            node_visited_cost_next.resize(node_list_array[for_edge + 1].size());
            for (auto node_cost : node_visited_cost_next)node_cost.clear();

            visited_node_current.clear();
            visited_node_current.assign(node_list_array[for_edge + 1].size(), false);
        }

        //cout<<for_edge<<"("<<edgeCost<<","<<sumEdgeCost<<"):"<<bddstack_current.size()<<endl;

        while (1) {
            currentTop = bddstack_current.top();
            bddstack_current.pop();

            //cout<<currentTop.cost<<","<<pqBDDNode.size()<<endl;

            //if(currentTop.cost >= minimumCost)break;

            // BDD uses edgeID i
            int pointerToNextNode = node_list_array[for_edge][currentTop.nodeId].existencePointer;
            if (pointerToNextNode == -1) {
                for(auto cost : node_visited_cost_current[currentTop.nodeId]) {
                    if (cost + edgeCost > threshold)continue;
                    treeNum++;
                    averageCost+=cost+edgeCost;
                    cout << "T" << treeNum << "," << cost + edgeCost << endl;

                    if (minimumCost > cost + edgeCost) {
                        minimumCost = cost + edgeCost;
                        //goalNode = make_pair(currentTop.edgeIdinBDD, currentTop.nodeId);
                    }
                }
            } else if (pointerToNextNode != -2) {
                //double visited_cost = currentTop.cost[0] + edgeCost;
                //if (node_visited_cost[currentTop.edgeIdinBDD + 1][pointerToNextNode] > visited_cost) {


                for(auto cost : node_visited_cost_current[currentTop.nodeId]){
                    //double cost = node_visited_cost_current[currentTop.nodeId][0];
                    if(cost+edgeCost<threshold){
                        //if(for_edge>=84) node_visited_cost_next[pointerToNextNode].push_back(cost + edgeCost);
                        node_visited_cost_next[pointerToNextNode].push_back(cost + edgeCost);
                        //else if(for_edge>=node_list_array.size()/2&&cost+edgeCost<threshold-(node_list_array.size()-for_edge)*8) node_visited_cost_next[pointerToNextNode].push_back(cost + edgeCost);
                    }
                }
                if(!visited_node_current[pointerToNextNode]&&!node_visited_cost_next[pointerToNextNode].empty()){
                    bddstack_next.push(BDDNodeTraversal(pointerToNextNode));
                    visited_node_current[pointerToNextNode]=true;
                }
                //node_visited_cost[currentTop.edgeIdinBDD + 1][pointerToNextNode] = visited_cost;
                    //previsited_node[currentTop.edgeIdinBDD + 1][pointerToNextNode] = currentTop.nodeId;
                //}
            }

            // BDD does not use edgeID i
            pointerToNextNode = node_list_array[for_edge][currentTop.nodeId].nonexistencePointer;
            if (pointerToNextNode == -1) {
                for(auto cost : node_visited_cost_current[currentTop.nodeId]) {
                    if(cost>threshold)continue;
                    treeNum++;
                    //cout << treeNum << "," << cost << endl;
                }
//                if (minimumCost > currentTop.cost) {
//                    minimumCost = currentTop.cost;
//                    goalNode = make_pair(currentTop.edgeIdinBDD, currentTop.nodeId);
//                }
            } else if (pointerToNextNode != -2) {


                for(auto cost : node_visited_cost_current[currentTop.nodeId]){
                    if(cost<threshold){
                         node_visited_cost_next[pointerToNextNode].push_back(cost);
//                        if(for_edge>=84) node_visited_cost_next[pointerToNextNode].push_back(cost);
//                        else if(for_edge<node_list_array.size()/2&&cost<threshold-(node_list_array.size()-for_edge)*7) {
//                            node_visited_cost_next[pointerToNextNode].push_back(cost);
//                        }
//                        //else if(for_edge>=node_list_array.size()/2) node_visited_cost_next[pointerToNextNode].push_back(cost);
//                        else if(for_edge>=node_list_array.size()/2&&cost<threshold-(node_list_array.size()-for_edge)*8) node_visited_cost_next[pointerToNextNode].push_back(cost);

                    }
                }
                if(!visited_node_current[pointerToNextNode]&&!node_visited_cost_next[pointerToNextNode].empty()){
                    bddstack_next.push(BDDNodeTraversal(pointerToNextNode));
                    visited_node_current[pointerToNextNode]=true;
                }
                //double visited_cost = currentTop.cost;
                //if (node_visited_cost[currentTop.edgeIdinBDD + 1][pointerToNextNode] > visited_cost) {
                //    bddstack_next.push(BDDNodeonPriorityQueue(pointerToNextNode, visited_cost));
                //    node_visited_cost[currentTop.edgeIdinBDD + 1][pointerToNextNode] = visited_cost;
                //    previsited_node[currentTop.edgeIdinBDD + 1][pointerToNextNode] = currentTop.nodeId;
                //}
            }
            if (bddstack_current.empty())break;
        }

        if(bddstack_next.empty())break;
        else {
            bddstack_current=bddstack_next;
            bddstack_next=stack<BDDNodeTraversal>();

            for(int i=0; i < node_visited_cost_next.size();i++){
                //cout<<node_visited_cost_next[i].size()<<",";
                if(node_visited_cost_next[i].size()>=steinernum){
                    sort(node_visited_cost_next[i].begin(),node_visited_cost_next[i].end());
                    node_visited_cost_next[i].resize(steinernum);
                }
            }
//            cout<<endl;
//            for(auto costlist : node_visited_cost_next){
//                cout<<costlist.size()<<",";
//            }
//            cout<<endl;
            node_visited_cost_current=node_visited_cost_next;
            node_visited_cost_next.clear();
        }
    }
    cout<<"bddcost="<<minimumCost<<endl;

    //////////// //////////// //////////// ////////////
    //////////// certain Graph construction////////////////
//    double sumcost=0;
//    vector<Edge> edge_list;
//    edge_list.clear();
//
//    Edge ucedge = (*state->graph->GetEdgeList())[goalNode.first];
//    Edge cedge(ucedge.src,ucedge.dest,ucedge.weight);
//    edge_list.push_back(cedge);
//
//    sumcost+=cedge.weight;
//    int prenodeid = previsited_node[goalNode.first][goalNode.second];
//    int currentnodeid = goalNode.second;
//    for(int i=goalNode.first-1;i>=0;--i){
//        ucedge = (*state->graph->GetEdgeList())[i];
//        int e_pointer = node_list_array[i][prenodeid].existencePointer;
//        int ne_pointer = node_list_array[i][prenodeid].nonexistencePointer;
//        if(e_pointer==currentnodeid && ne_pointer != currentnodeid){
//            Edge cedge(ucedge.src,ucedge.dest,ucedge.weight);
//            edge_list.push_back(cedge);
//            sumcost+= cedge.weight;
//            //cout<<i<<":"<<endl;
//        }
//        if(i==0)break;
//        currentnodeid = prenodeid;
//        prenodeid = previsited_node[i][prenodeid];
//    }
//    cout<<"edge weight sum = "<<sumcost<<endl;
//
    OriginalGraph steiner_trees = OriginalGraph();
//    steiner_trees.number_of_vertices=state->graph->number_of_vertices;
//    steiner_trees.number_of_edges=edge_list.size();
//    steiner_trees.number_of_terminals=state->terminals.size();
//    steiner_trees.terminals = state->terminals;
//    steiner_trees.edge_list = edge_list;
//
//    steiner_trees.vertex_neighbor_edge.resize(state->graph->number_of_vertices);
//    for(int i=0;i<state->graph->number_of_vertices;i++)steiner_trees.vertex_neighbor_edge[i].clear();
//
//    sumcost=0;
//    for(int i=0;i<steiner_trees.number_of_edges;i++){
//        int v1 = edge_list[i].src;
//
//        int v2 = edge_list[i].dest;
//        steiner_trees.vertex_neighbor_edge[v1].push_back(i);
//        steiner_trees.vertex_neighbor_edge[v2].push_back(i);
//        cout<<edge_list[i].src<<"->"<<edge_list[i].dest<<endl;
//    }

    //////////// certain Graph construction////////////////
    //////////// //////////// //////////// //////////// //

    this->minimumCost=minimumCost;
    this->averageCost=averageCost/treeNum;
    this->treenum=treeNum;
    return steiner_trees;
}


OriginalGraph BDD::MinimumTraversal(State* state){

    vector<vector<double>> node_visited_cost;
    vector<vector<int>> previsited_node;
    pair<int, int> goalNode;//<edgeID, nodeId>

    node_visited_cost.resize(node_list_array.size());
    previsited_node.resize(node_list_array.size());

    for(unsigned int i = 0; i < node_list_array.size();i++){
        node_visited_cost[i].resize(node_list_array[i].size());
        previsited_node[i].resize(node_list_array[i].size());
        for(unsigned int j = 0; j < node_list_array[i].size();j++){
            node_visited_cost[i][j] = std::numeric_limits<double>::max();
            previsited_node[i][j] = -1;
        }
    }

    priority_queue<BDDNodeonPriorityQueue,std::vector<BDDNodeonPriorityQueue>,std::greater<BDDNodeonPriorityQueue>> pqBDDNode;
    //stack<BDDNodeonPriorityQueue> pqBDDNode;
    pqBDDNode.push(BDDNodeonPriorityQueue(0,0,0));
    node_visited_cost[0][0]=0;
    BDDNodeonPriorityQueue currentTop;
    double minimumCost =  std::numeric_limits<double>::max();
    double edgeCost;

    //for(unsigned int i = 0; i < node_list_array.size();i++){//node_list_array.size = # of edges.
    int treeNum=0;
    while(1){
        currentTop = pqBDDNode.top();
        pqBDDNode.pop();
        Edge edge = (*state->graph->GetEdgeList())[currentTop.edgeIdinBDD];
        edgeCost = edge.weight;

        if(currentTop.cost >= minimumCost)break;

        // BDD uses edgeID i
        int pointerToNextNode = node_list_array[currentTop.edgeIdinBDD][currentTop.nodeId].existencePointer;
        if(pointerToNextNode == -1){
            treeNum++;
            cout<<treeNum<<","<<currentTop.cost+edgeCost<<endl;
            if(minimumCost > currentTop.cost+edgeCost){
                minimumCost = currentTop.cost+edgeCost;
                goalNode=make_pair(currentTop.edgeIdinBDD,currentTop.nodeId);
            }
        }
        else if(pointerToNextNode != -2) {
            double visited_cost= currentTop.cost + edgeCost;
            if (node_visited_cost[currentTop.edgeIdinBDD + 1][pointerToNextNode] > visited_cost) {
                pqBDDNode.push(BDDNodeonPriorityQueue(currentTop.edgeIdinBDD + 1, pointerToNextNode, visited_cost));
                node_visited_cost[currentTop.edgeIdinBDD + 1][pointerToNextNode] =  visited_cost;
                previsited_node[currentTop.edgeIdinBDD + 1][pointerToNextNode] =  currentTop.nodeId;
            }
        }

        // BDD does not use edgeID i
        pointerToNextNode = node_list_array[currentTop.edgeIdinBDD][currentTop.nodeId].nonexistencePointer;
        if(pointerToNextNode == -1){
            cout<<treeNum<<","<<currentTop.cost<<endl;
            if(minimumCost > currentTop.cost) {
                minimumCost = currentTop.cost;
                goalNode = make_pair(currentTop.edgeIdinBDD, currentTop.nodeId);
            }
        }
        else if(pointerToNextNode != -2){
            double visited_cost= currentTop.cost;
            if (node_visited_cost[currentTop.edgeIdinBDD + 1][pointerToNextNode] > visited_cost) {
                pqBDDNode.push(BDDNodeonPriorityQueue(currentTop.edgeIdinBDD + 1, pointerToNextNode, visited_cost));
                node_visited_cost[currentTop.edgeIdinBDD + 1][pointerToNextNode] = visited_cost;
                previsited_node[currentTop.edgeIdinBDD + 1][pointerToNextNode] =  currentTop.nodeId;
            }
        }
        if(pqBDDNode.empty())break;
    }
    cout<<"bddcost="<<minimumCost<<endl;

    //////////// //////////// //////////// ////////////
    //////////// certain Graph construction////////////////
    double sumcost=0;
    vector<Edge> edge_list;
    edge_list.clear();

    Edge ucedge = (*state->graph->GetEdgeList())[goalNode.first];
    Edge cedge(ucedge.src,ucedge.dest,ucedge.weight);
    edge_list.push_back(cedge);

    sumcost+=cedge.weight;
    int prenodeid = previsited_node[goalNode.first][goalNode.second];
    int currentnodeid = goalNode.second;
    for(int i=goalNode.first-1;i>=0;--i){
        ucedge = (*state->graph->GetEdgeList())[i];
        int e_pointer = node_list_array[i][prenodeid].existencePointer;
        int ne_pointer = node_list_array[i][prenodeid].nonexistencePointer;
        if(e_pointer==currentnodeid && ne_pointer != currentnodeid){
            Edge cedge(ucedge.src,ucedge.dest,ucedge.weight);
            edge_list.push_back(cedge);
            sumcost+= cedge.weight;
            cout<<i+1<<":"<<cedge.weight<<endl;
        }
        if(i==0)break;
        currentnodeid = prenodeid;
        prenodeid = previsited_node[i][prenodeid];
    }
    cout<<"edge weight sum = "<<sumcost<<endl;

    OriginalGraph certaingraph = OriginalGraph();
    certaingraph.number_of_vertices=state->graph->number_of_vertices;
    certaingraph.number_of_edges=edge_list.size();
    certaingraph.number_of_terminals=state->terminals.size();
    certaingraph.terminals = state->terminals;
    certaingraph.edge_list = edge_list;

    certaingraph.vertex_neighbor_edge.resize(state->graph->number_of_vertices);
    for(int i=0;i<state->graph->number_of_vertices;i++)certaingraph.vertex_neighbor_edge[i].clear();

    sumcost=0;
    for(int i=0;i<certaingraph.number_of_edges;i++){
        int v1 = edge_list[i].src;
        int v2 = edge_list[i].dest;
        certaingraph.vertex_neighbor_edge[v1].push_back(i);
        certaingraph.vertex_neighbor_edge[v2].push_back(i);
        cout<<edge_list[i].src<<"->"<<edge_list[i].dest<<endl;
    }

    //////////// certain Graph construction////////////////
    //////////// //////////// //////////// //////////// //
    return certaingraph;
}


int BDD::GetSize(){

    vector<vector<bool>> reducedFlags;

    int BDDSize=0;
    for(int for_edge=node_list_array.size()-1; for_edge>0;for_edge--) {

        for(int for_node=0;for_node<node_list_array[for_edge].size();for_node++) {
            bool existencePointerToFalse=false;
            bool nonexistencePointerToFalse=false;

            int pointerToNextNode = node_list_array[for_edge][for_node].existencePointer;
            int non_pointerToNextNode = node_list_array[for_edge][for_node].nonexistencePointer;

            if(pointerToNextNode!=-2||non_pointerToNextNode!=-2){
                BDDSize++;
            }
        }
    }

    return BDDSize;
}
