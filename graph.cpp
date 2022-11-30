
#include "graph.h"



void OriginalGraph::InputFile(std::string input_file_name)
{

    const char *graphfile=input_file_name.c_str();
    std::ifstream input(graphfile);

    if(!input){
        std::cout<<"error: cannot open graph file"<<std::endl;
        exit(1);
    }

    input>>number_of_vertices>>number_of_edges>>number_of_terminals>>correct_answer;
    edge_list.clear();


    degrees.resize(number_of_vertices);
    for(int i=0;i<number_of_vertices;i++)degrees[i]=0;


    int v1,v2;
    double w,p;
    for(int i =0; i < number_of_edges;i++){

        input >> v1 >> v2 >> w;
        if(input.fail())break;

        if(v1>v2){
            int tempv;
            tempv=v1;
            v1=v2;
            v2=tempv;
        }
        Edge edge(v1, v2, w);

        edge_list.push_back(edge);
        degrees[v1]++;
        degrees[v2]++;
    }

    int t;
    for(int i =0; i < number_of_terminals;i++){
        input >> t;
        terminals.push_back(t);
    }

    vertex_neighbor_edge.resize(number_of_vertices);
    for(int i=0;i<number_of_vertices;i++)vertex_neighbor_edge[i].clear();

    for(int i=0;i<number_of_edges;i++){
        v1 = edge_list[i].src;
        v2 = edge_list[i].dest;
        vertex_neighbor_edge[v1].push_back(i);
        vertex_neighbor_edge[v2].push_back(i);
    }

    for(int i=0;i<number_of_vertices;i++){
        if(max_degree<degrees[i])max_degree=degrees[i];
    }

    for(int tc=0;tc<terminals.size()-1;tc++){
        //cout<<terminals[tc]<<"'s degrees = "<<degrees[terminals[tc]]<<endl;
        for(int tcc=tc+1;tcc<terminals.size();tcc++){
            if(degrees[terminals[tc]] > degrees[terminals[tcc]]){
                int tempter = terminals[tc];
                terminals[tc] =terminals[tcc];
                terminals[tcc] = tempter;
            }
        }
    }
}



void SeedSteinerTreeSet::InputFile(OriginalGraph* graph_, vector<string> filenames_)
{

    number_of_vertices=graph_->number_of_vertices;
	number_of_terminals=graph_->number_of_terminals;
	correct_answer=graph_->correct_answer;

    terminals=graph_->terminals;

    degrees.clear();
    degrees.assign(number_of_vertices,0);

    number_of_seedtrees=filenames_.size();
    seedtrees_answer.clear();

    vector<bool> existEdges;
    existEdges.assign(graph_->number_of_edges,false);

    for(int i =0;i < filenames_.size();i++){
        std::ifstream input(filenames_[i]);

        if(!input){
            std::cout<<"error: cannot open seed tree file"<< filenames_[i] <<std::endl;
            exit(1);
        }

        string tmp;
        bool firstline=true;
        while(getline(input, tmp)){

            if(firstline){
                seedtrees_answer.push_back(stod(tmp));
                firstline=false;
            }
            else{
                existEdges[stoi(tmp)]=true;
            }
        }
    }

    edge_list.clear();
    for(int i=0;i < existEdges.size();i++){

        if(existEdges[i]){
            edge_list.push_back(graph_->edge_list[i]);
        }
    }

    int v1,v2;
    vertex_neighbor_edge.resize(number_of_vertices);
    for(int i=0;i<number_of_vertices;i++)vertex_neighbor_edge[i].clear();

    for(int i=0;i<number_of_edges;i++){
        v1 = edge_list[i].src;
        v2 = edge_list[i].dest;
        vertex_neighbor_edge[v1].push_back(i);
        vertex_neighbor_edge[v2].push_back(i);
        degrees[v1]++;
        degrees[v2]++;
    }

    for(int i=0;i<number_of_vertices;i++){
        if(max_degree<degrees[i])max_degree=degrees[i];
    }

    for(int tc=0;tc<terminals.size()-1;tc++){
        //cout<<terminals[tc]<<"'s degrees = "<<degrees[terminals[tc]]<<endl;
        for(int tcc=tc+1;tcc<terminals.size();tcc++){
            if(degrees[terminals[tc]] > degrees[terminals[tcc]]){
                int tempter = terminals[tc];
                terminals[tc] =terminals[tcc];
                terminals[tcc] = tempter;
            }
        }
    }

}


///////////////////////////////////////////
/////Prune/////////////////////////////
///////////////////////////////////////////

void OriginalGraph::GraphReduction() {


    std::cout<<"GraphReduction Start: "<<edge_list.size()<<std::endl;

    bool continueFlag=false;
    vector<bool> edgeInsertFlag;
    int count=edge_list.size();
    edgeInsertFlag.resize(edge_list.size());
    vector< vector < pair <int,int> >> adjacencyList;
    vector<int> checkedAdjacencyCount(number_of_vertices,0);

    adjacencyList.resize(number_of_vertices);

    for(int i=0;i<edge_list.size();i++){
        if(edge_list[i].src!=edge_list[i].dest){
            adjacencyList[edge_list[i].src].push_back(make_pair(edge_list[i].dest,i));

            adjacencyList[edge_list[i].dest].push_back(make_pair(edge_list[i].src,i));
            edgeInsertFlag[i]=true;
        }
        else edgeInsertFlag[i]=false;
    }

    //cout<<number_of_vertices<<","<<degrees.size()<<endl;

    while(1) {
        continueFlag = false;

        //cout << "1 : " << count << endl;

        for (int i = 0; i < number_of_vertices; i++) {

            if (adjacencyList[i].size() == 2 && !FindVector(terminals, i)) {

                if(adjacencyList[i][0].first==adjacencyList[i][1].first)continue;
                //if(FindVector(terminals, adjacencyList[i][0].first)||FindVector(terminals, adjacencyList[i][1].first))continue;

                edgeInsertFlag[adjacencyList[i][0].second] = false;
                edgeInsertFlag[adjacencyList[i][1].second] = false;
                int vertex1 = adjacencyList[i][0].first;
                int vertex2 = adjacencyList[i][1].first;


                if (vertex1 != vertex2) {
                    Edge i_edge(vertex1, vertex2,0);
                    i_edge.weight=edge_list[adjacencyList[i][0].second].weight +
                                                  edge_list[adjacencyList[i][1].second].weight;//check
                    edge_list.push_back(i_edge);
                    edgeInsertFlag.push_back(true);
                    adjacencyList[vertex1].push_back(make_pair(vertex2, edge_list.size() - 1));

                    adjacencyList[vertex2].push_back(make_pair(vertex1, edge_list.size() - 1));
                    for(int j=0;j<adjacencyList[vertex1].size();j++){
                        if(adjacencyList[vertex1][j].first==i){
                             adjacencyList[vertex1].erase(adjacencyList[vertex1].begin()+j);
                        }
                    }
                    for(int j=0;j<adjacencyList[vertex2].size();j++){
                        if(adjacencyList[vertex2][j].first==i){
                             adjacencyList[vertex2].erase(adjacencyList[vertex2].begin()+j);
                        }
                    }
                }
                else{
                    cout<<"!!!!"<<endl;
                }
                adjacencyList[i].clear();
                continueFlag=true;
                //cout<<edge_list.size()<<","<<edgeInsertFlag.size()<<endl;
                count--;
            }
        }


        //cout << "2 : " << count << endl;

        for (int i = 0; i < number_of_vertices; i++) {

            for (int j = 0; j < adjacencyList[i].size(); j++) {

                if(i == adjacencyList[i][j].first){
                    edgeInsertFlag[adjacencyList[i][j].second]=false;
                    adjacencyList[i].erase(adjacencyList[i].begin()+j);
                    j--;
                    count--;
                    continue;
                }


            }
            checkedAdjacencyCount[i]=adjacencyList[i].size();

        }

        //cout<<"3 : "<<count<<endl;
        //break;
        if(!continueFlag)break;
    }

    vector<Edge> tempEdgeList;

    //std::cout<<"GraphReduction Done: "<<edge_list.size()<<std::endl;

    for(int i=0;i < edge_list.size();i++){

        if(edgeInsertFlag[i]){
            tempEdgeList.push_back(edge_list[i]);
        }
    }
    edge_list.clear();
    for(int i=0;i < tempEdgeList.size();i++){
        edge_list.push_back(tempEdgeList[i]);
    }


    number_of_edges=edge_list.size();

    vertex_neighbor_edge.resize(number_of_vertices);
    for(int i=0;i<number_of_vertices;i++)vertex_neighbor_edge[i].clear();

    for(int i=0;i<number_of_edges;i++){
        int v1 = edge_list[i].src;
        int v2 = edge_list[i].dest;
        vertex_neighbor_edge[v1].push_back(i);
        vertex_neighbor_edge[v2].push_back(i);
    }

 //   edge_list=tempEdgeList;

    std::cout<<"GraphReduction Done: "<<edge_list.size()<<std::endl;

}

