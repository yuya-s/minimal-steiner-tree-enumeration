
#include "state.h"


struct Comp {
    constexpr bool operator()(pair<int, long double> const & a,
                              pair<int, long double> const & b) const noexcept
    { return a.second < b.second; }
};

struct Comp_weight {
    constexpr bool operator()(pair<int, double> const & a,
                              pair<int, double> const & b) const noexcept
    { return a.second > b.second; }
};

void State::ComputeFrontier(int i){

    vector<Edge>* edge_list = graph->GetEdgeList();
    Fp.clear();
        // copy frontier info from previous to current frontier info
    for (unsigned int j = 0; j < Fc.size(); ++j) {
        Fp.push_back(Fc[j]);
    }
    Edge edge = (*edge_list)[i-1];
    int src = edge.src;
    int dest = edge.dest;

    // Add frontiers
    if (!Contains(Fc, src)) {
        // if not contains, add src as new frontiers
        Fc.push_back(src);
    }
    if (!Contains(Fc, dest))
    {
         // if not contains, add dest as new frontiers
        Fc.push_back(dest);
    }


    //Delete frontiers
    if (!FindElement(i-1, src))
    {
        //Delete src if src does not appear in unprocessed edges
        Remove(&Fc, src);
    }
    if (!FindElement(i-1, dest))
    {
         //Delete dest if dest does not appear in unprocessed edges
        Remove(&Fc, dest);
    }
    sort(Fc.begin(),Fc.end());
    //for(int j=0;j<Fc.size();j++)cout<<"Frontier "<<i<<":"<<Fc[j]<<endl;
}



void State::ComputeFrontier() //Computing frontier information
{
    vector<Edge>* edge_list = graph->GetEdgeList();

    F = new vector<int>*[edge_list->size() + 1];
    F[0] = new vector<int>;

    for (unsigned int i = 0; i < edge_list->size(); ++i)
    {
        F[i + 1] = new vector<int>;

        for (unsigned int j = 0; j < F[i]->size(); ++j) {
            F[i + 1]->push_back((*F[i])[j]);
        }

        Edge edge = (*edge_list)[i];
        int src = edge.src;
        int dest = edge.dest;

        if (!Contains(*F[i + 1], src))
        {
            F[i + 1]->push_back(src);
        }
        if (!Contains(*F[i + 1], dest))
        {
            F[i + 1]->push_back(dest);
        }

        if (!FindElement(i, src))
        {
            Remove(F[i + 1], src);
        }
        if (!FindElement(i, dest))
        {
            Remove(F[i + 1], dest);
        }
        sort(F[i+1]->begin(),F[i+1]->end());
        //for(int j=0;j<F[i+1]->size();j++)cout<<"Frontier i="<<i<<":"<<(*F[i+1])[j]<<endl;
    }

}

BDDinputGraph* State::CopyGraph(SeedSteinerTreeSet* original_graph, vector<int> terminals_) {

    BDDinputGraph* graph_ = new BDDinputGraph();
    graph_->number_of_vertices = original_graph->number_of_vertices;
    graph_->max_degree = original_graph->max_degree;
    graph_->edge_list = original_graph->edge_list;
    graph_->degrees = original_graph->degrees;
    graph_->terminals.resize(original_graph->number_of_vertices);

    for(int i=0;i<graph_->number_of_vertices;i++){
        graph_->terminals[i]=false;
    }

    for(int i=0;i<terminals_.size();i++){
        graph_->terminals[terminals_[i]]=true;
    }

    return graph_;
}

BDDinputGraph* State::CopyGraphBFSbased(SeedSteinerTreeSet* original_graph, vector<int> terminals_){

    BDDinputGraph* graph_ = new BDDinputGraph();
    graph_->number_of_vertices = original_graph->number_of_vertices;
    graph_->max_degree = original_graph->max_degree;
    graph_->edge_list.clear();
    graph_->degrees.resize(original_graph->number_of_vertices);
    graph_->terminals.resize(original_graph->number_of_vertices);


    vector<bool> vertexInsertedList;
    vertexInsertedList.resize(original_graph->number_of_vertices);

    for(int i=0;i<graph_->number_of_vertices;i++){
        vertexInsertedList[i]=false;
        graph_->degrees[i]=0;
        graph_->terminals[i]=false;
    }
    int initial_vertex_id = terminals_[0];
    int minDegree=graph_->degrees[terminals_[0]];
    cout<<minDegree<<endl;

    vertexInsertedList[initial_vertex_id]=true;

    for(int i=0;i<terminals_.size();i++){
        graph_->terminals[terminals_[i]]=true;
    }

    vector<bool> edgeInsertedList;
    edgeInsertedList.resize(original_graph->edge_list.size());

    for(int i=0;i<original_graph->edge_list.size();i++){
        edgeInsertedList[i]=false;
    }


    std::queue<int> vertexQueue;
    vertexQueue.push(initial_vertex_id);
    int currentVertexId;
    vector< vector <int> >vertex_neighbor_edge;
    vertex_neighbor_edge.resize(original_graph->number_of_vertices);

    for(int i=0;i<original_graph->edge_list.size();i++) {
        int v1 = original_graph->edge_list[i].src;
        int v2 = original_graph->edge_list[i].dest;
        vertex_neighbor_edge[v1].push_back(i);
        vertex_neighbor_edge[v2].push_back(i);
    }


    while(1){

        currentVertexId = vertexQueue.front();
        vertexQueue.pop();

        for(int i =0; i <vertex_neighbor_edge[currentVertexId].size();i++){

            Edge edge_ = original_graph->edge_list[vertex_neighbor_edge[currentVertexId][i]];
            if(!edgeInsertedList[vertex_neighbor_edge[currentVertexId][i]]){
                graph_->edge_list.push_back(edge_);
                edgeInsertedList[vertex_neighbor_edge[currentVertexId][i]]=true;
                graph_->degrees[edge_.src]++;
                graph_->degrees[edge_.dest]++;
            }

            if(currentVertexId == edge_.src){
                if(!vertexInsertedList[edge_.dest]) {
                    vertexQueue.push(edge_.dest);
                    vertexInsertedList[edge_.dest] = true;
                }
            }
            else{
                if(!vertexInsertedList[edge_.src]) {
                    vertexQueue.push(edge_.src);
                    vertexInsertedList[edge_.src] = true;

                }
            }
        }
        if(vertexQueue.empty())break;

    }


    edgeInsertedList.resize(original_graph->edge_list.size());

    for(int i=0;i<original_graph->edge_list.size();i++){
        edgeInsertedList[i]=false;
    }

    vector<int> previousVertex;
    vector<int> previousEdge;
    vector<double> weight;
    previousVertex.resize(original_graph->number_of_vertices);
    previousEdge.resize(original_graph->number_of_vertices);
    weight.resize(original_graph->number_of_vertices);
    double currentVertexWeight;


    std::priority_queue<std::pair<int,long double>,vector<pair<int,double>>,Comp_weight> vertexQueue_p;

    for(int t=1;t<terminals_.size();t++) {
        int target_vertex = terminals_[t];

        for (int i = 0; i < original_graph->number_of_vertices; i++) {
            previousVertex[i] = -1;
            previousEdge[i] = -1;
            weight[i] = std::numeric_limits<double>::max();
        }

        previousVertex[initial_vertex_id] = -2;
        weight[initial_vertex_id] = 0;

        vertexQueue_p = priority_queue<std::pair<int,long double>,vector<pair<int,double>>,Comp_weight>();
        vertexQueue_p.push(make_pair(initial_vertex_id,1));

        while (1) {

            currentVertexId = vertexQueue_p.top().first;
            currentVertexWeight = vertexQueue_p.top().second;

            vertexQueue_p.pop();

            if (weight[target_vertex] <= currentVertexWeight && previousVertex[target_vertex] > -1){break;}

            for (int i = 0; i < vertex_neighbor_edge[currentVertexId].size(); i++) {

                Edge edge_ = original_graph->edge_list[vertex_neighbor_edge[currentVertexId][i]];

                if (currentVertexId == edge_.src) {

                    if (currentVertexWeight + edge_.weight < weight[edge_.dest]) {
                        weight[edge_.dest] = currentVertexWeight + edge_.weight;
                        previousVertex[edge_.dest] = currentVertexId;
                        previousEdge[edge_.dest] = vertex_neighbor_edge[currentVertexId][i];
                        vertexQueue_p.push(make_pair(edge_.dest, weight[edge_.dest]));
                        //cout<<currentVertexId<<"x"<<edge_.dest<<endl;
                    }
                } else {

                    if (currentVertexWeight + edge_.weight < weight[edge_.src]) {
                        weight[edge_.src] = currentVertexWeight + edge_.weight;
                        previousVertex[edge_.src] = currentVertexId;
                        previousEdge[edge_.src] =  vertex_neighbor_edge[currentVertexId][i];
                        vertexQueue_p.push(make_pair(edge_.src, weight[edge_.src]));
                        //cout<<currentVertexId<<"x"<<edge_.src<<endl;
                    }
                }
            }
            if (vertexQueue_p.empty()){break;}
        }
        int backVertex=target_vertex;
        stack<int> backEdge;
        while(1){
            int edgeId = previousEdge[backVertex];
            backVertex=previousVertex[backVertex];
            //cout<<edgeId<<endl;
            if(!edgeInsertedList[edgeId]) {
                 edgeInsertedList[edgeId] = true;
            }
            if(backVertex==initial_vertex_id)break;
        }
    }

    return graph_;

}

bool State::FindElement(int edge_number, int value)
{
    vector<Edge>* edge_list = graph->GetEdgeList();
    for (unsigned int i = edge_number + 1; i < edge_list->size(); ++i)
    {
        if (value == (*edge_list)[i].src || value == (*edge_list)[i].dest)
        {
            return true;
        }
    }
    return false;
}

