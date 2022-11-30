
#include "frontier.h"
#include "bdd.h"




int globalCount=0;

long double findTime=0;
long double checkTerminalTime=0;
long double updateTime=0;
long double makecopyTime=0;

//******************************************************************************
using namespace std;
typedef vector<int> frontierInfo;

namespace std {
  template<>
  class hash<frontierInfo> {
  public:
    size_t operator()(const frontierInfo &p) const {
        size_t seed = 0;

        auto t_hash =0;


        int size = p.size();
        for(int i=0;i < size;i++) {
            t_hash += hash<int>()(p[i]);

        }

        seed ^= t_hash + 0x9e3779b9 + (seed << 6) + (seed >> 2);

        return seed;
    }
  };
}

bool cmp(const FrontierNode* a, const FrontierNode* b)
{
    return a->number_of_existence < b->number_of_existence;
}



BDD* FrontierAlgorithm::ComputeSteinertree(State *state, unsigned long int max_elements, double _error, double weight_){

    //initialize global variants for time//
    struct timespec startTime, endTime;
    long double totalTime=0;
    long double connectCount=0;


    bool layerSampleFlag=false;

    list<double> probabilityList;

    bool pushFlag=false;



    clock_gettime(CLOCK_REALTIME, &startTime);
    findTime=0;
    checkTerminalTime=0;
    updateTime=0;
    makecopyTime=0;

    std::random_device rnd;
    std::uniform_real_distribution<> rand(0,1);
    mt19937 mt(rnd());

    vector<Edge>* edge_list = state->graph->GetEdgeList();
    unsigned int edgeSize = edge_list->size();

    priority_queue < FrontierNode*, vector<FrontierNode*>, greater<FrontierNode*> > pq_current;
    priority_queue < FrontierNode*, vector<FrontierNode*>, greater<FrontierNode*> > pq_next;

    vector<FrontierNode*> N_current;
    vector<FrontierNode*> N_next;



    FrontierNode* n_root = FrontierNode::CreateRootNode(state)->MakeCopy();

    N_current.push_back(n_root);



    google::dense_hash_map<long int, int> hashSpace;
    hashSpace.set_empty_key(-1);




    BDD* bdd = new BDD;
    bdd->Initialize(edgeSize);

    BDDNode bddnode=BDDNode();

    bdd->node_list_array[0].push_back(bddnode);

    for (unsigned int i = 1; i <= edgeSize; ++i) { // for edge i
        state->ComputeFrontier(i);


        findTime=0;
        checkTerminalTime=0;
        updateTime=0;
        makecopyTime=0;
        int hashCol=0;

        unsigned int size=N_current.size();
        //int size=pq_current.size();
        Edge edge = (*state->graph->GetEdgeList())[i - 1];

        std::cout<<"["<<i<<" edge weight="<< edge.weight<<"] size : "<<size<<", #frontier"<<state->Fc.size()<<endl;

        //cout<<edge.src<<"->"<<edge.dest<<":"<<edge.weight<<endl;

        if(max_elements<size)size=max_elements;

        bool breakFlag=false;
        BDDNode bddnode=BDDNode();

        for (unsigned int j = 0; j < size; ++j) { // process for layer i edge
            FrontierNode* n_hat = N_current.back();




            for (int x = 1; x >= 0; --x) { // process for edge fixing/delete
                struct timespec checkterminal_startTime, checkterminal_endTime;


                clock_gettime(CLOCK_REALTIME, &checkterminal_startTime);
                //cout<<"j="<<j<<","<<x<<endl;
                FrontierNode* n_prime = CheckTerminal(n_hat, i, x, state, bdd->node_list_array[i-1][bdd->node_list_array[i-1].size()-1-j], weight_);
                clock_gettime(CLOCK_REALTIME, &checkterminal_endTime);

                checkTerminalTime += (checkterminal_endTime.tv_sec - checkterminal_startTime.tv_sec) + (checkterminal_endTime.tv_nsec - checkterminal_startTime.tv_nsec) * pow(10, -9);

                if (n_prime != NULL) { // not to sinks

                    struct timespec find_startTime, find_endTime;
                    clock_gettime(CLOCK_REALTIME, &find_startTime);

                    vector<int>& frontier = state->Fc;
                    int fSize = frontier.size();

                    long int hashValue=17;
                    int tnum=0;

                    bool testflag=false;
                    long double heuristicFunction;
                    for(int f_size =0;f_size < fSize;f_size++){
                        hashValue*=31;
                        hashValue+=n_prime->comp[f_size];
                        if(testflag)cout<<hashValue<<endl;
                        hashValue*=31;
                        if(testflag)cout<<hashValue<<endl;
                        if(n_prime->terminal[f_size]>0) {
                            hashValue+=17;
                            hashValue*=31;
                            if(testflag)cout<<hashValue<<endl;
                        }
                        hashValue+=n_prime->addedge[f_size];
                        hashValue*=31;

                    }

                    int n_primeprime = hashSpace[hashValue];

                    clock_gettime(CLOCK_REALTIME, &find_endTime);

                    if (n_primeprime != 0)
                    {


                        if(x==0)bdd->node_list_array[i-1][bdd->node_list_array[i-1].size()-1-j].nonexistencePointer=n_primeprime-1;
                        else if(x==1)bdd->node_list_array[i-1][bdd->node_list_array[i-1].size()-1-j].existencePointer=n_primeprime-1;

                        N_next[n_primeprime-1]->number_of_existence=max(n_prime->number_of_existence,N_next[n_primeprime-1]->number_of_existence);
                        N_next[n_primeprime-1]->weight=min(n_prime->weight,N_next[n_primeprime-1]->weight);

                        delete n_prime;
                    }
                    else
                    {

                       int size_Next = pq_next.size();
                       if(size_Next >= max_elements && max_elements!=0){
                           cout<<"over maximum width: width = "<< size_Next<<endl;
                       }
                       n_prime->HueristicValue_num_of_existence();
                       n_prime->nodeId = bdd->node_list_array[i].size();

                       N_next.push_back(n_prime);

                       BDDNode newbddnode=BDDNode();
                       if(x==0)bdd->node_list_array[i-1][bdd->node_list_array[i-1].size()-1-j].nonexistencePointer= bdd->node_list_array[i].size();
                       else if(x==1)bdd->node_list_array[i-1][bdd->node_list_array[i-1].size()-1-j].existencePointer= bdd->node_list_array[i].size();
                       if(globalCount==2374)bdd->node_list_array[i-1][bdd->node_list_array[i-1].size()-1-j].nonexistencePointer= -2;

                        bdd->node_list_array[i].push_back(newbddnode);

                       hashSpace[hashValue] = bdd->node_list_array[i].size();

                    }
                }

            }

            if(n_hat !=NULL){delete n_hat;n_hat=NULL;}
            N_current.pop_back();
            if(breakFlag)break;

        }

        N_current.clear();
        N_current = N_next;
        N_next.clear();
        hashSpace.clear();


        clock_gettime(CLOCK_REALTIME, &endTime);
        totalTime = (endTime.tv_sec - startTime.tv_sec) + (endTime.tv_nsec - startTime.tv_nsec) * pow(10, -9);

        //cout << "Total time : "<< totalTime<<"(FindVector, Checkterminal (makecopy, update) : "<<findTime<<","<< checkTerminalTime<<" ("<<makecopyTime<<" , "<<updateTime<<" )"<<endl;

        if(N_current.empty())break;
        if(breakFlag){

            break;
        }
    }

    cout << "BDD construction Done"<<endl;


    return bdd;
}


FrontierNode* FrontierAlgorithm::CheckTerminal(FrontierNode* n_hat, int i, int x, State* state, BDDNode& bddnode, double weight_)
{
    Edge edge = (*state->graph->GetEdgeList())[i - 1];

    struct timespec makecopy_startTime, makecopy_endTime, update_startTime, update_endTime;



    clock_gettime(CLOCK_REALTIME, &update_startTime);


    FrontierNode* n_prime = FrontierAlgorithm::UpdateInfo(n_hat, i, x, state, bddnode, weight_);
    clock_gettime(CLOCK_REALTIME, &update_endTime);
    updateTime += (update_endTime.tv_sec - update_startTime.tv_sec) + (update_endTime.tv_nsec - update_startTime.tv_nsec) * pow(10, -9);


    if(x==1) {

        if(n_prime == NULL){

            bddnode.nonexistencePointer=-2;
            return NULL;
        }
        if(i == static_cast<int>(state->graph->GetEdgeList()->size())) {
             cout<<"output1: "<<n_prime->weight<<endl;
             delete n_prime;
             bddnode.existencePointer=-1;
             return NULL;
        }

        for(int j=0;j<state->Fc.size();++j){
            if (n_prime->terminal[j]==state->terminals.size()) {

                bddnode.existencePointer=-1;
                cout<<"output2: "<<n_prime->weight<<endl;
                delete n_prime;
                return NULL;
            }
        }
    }
    else {
        if(n_prime == NULL) {

            bddnode.nonexistencePointer=-2;
            return NULL;

        }
        else if(i == static_cast<int>(state->graph->GetEdgeList()->size())){
            bddnode.nonexistencePointer=-2;
            if(n_prime!=NULL)delete n_prime;
            return NULL;
        }
    }


    return n_prime;

}


FrontierNode* FrontierAlgorithm::UpdateInfo(FrontierNode *n_hat, int i, int x, State *state, BDDNode& bddnode, double weight_)
{
    FrontierNode* node = new FrontierNode();

    int fSize = state->Fc.size();
	node->deg.resize(fSize);// = new int[fSize];
	node->comp.resize(fSize);// = new int[fSize];
    node->terminal.resize(fSize);// =  new int[fSize];
    node->addedge.resize(fSize);

    for(int j=0;j<fSize;j++){
        node->deg[j]=0;// = new int[fSize];
	    node->comp[j]=0;// = new int[fSize];
        node->terminal[j]=0;// =  new int[fSize];
        node->addedge[j]=0;
    }
    node->number_of_existence = n_hat->number_of_existence+x;
    node->weight = n_hat->weight;

    globalCount++;


    int i1=0;
    int i2=0;

    int frontierSize = state->Fc.size();
    int previousFrontierSize = state->Fp.size();

    int u,v;

    while(1){

        if(i1>=previousFrontierSize){//already copy all previous frontiers to current frontiers
            while(1){// all new frontiers
                if(i2>=frontierSize)break;
                v = state->Fc[i2];
                node->comp[i2]=v;
                node->deg[i2]=state->graph->degrees[v];
                if(state->graph->terminals[v])node->terminal[i2]=1;
                else node->terminal[i2]=0;
                i2++;
            }
        }
        if(i2>=frontierSize)break;

        u = state->Fp[i1];//i1-th previous frontier id
        v = state->Fc[i2];//i2-th current frontier id

        if(u==v){// if same
            node->comp[i2]=n_hat->comp[i1];
            node->deg[i2]=n_hat->deg[i1];
            node->terminal[i2]=n_hat->terminal[i1];
            node->addedge[i2]=n_hat->addedge[i1];
            i1++;
            i2++;
        }
        else if(u>v){//if v is not included in the set of previous frontiers
            node->comp[i2]=v;
            node->deg[i2]=state->graph->degrees[v];
            if(state->graph->terminals[v])node->terminal[i2]=1;
            else node->terminal[i2]=0;
            i2++;
            if(i2>=frontierSize)break;
        }
        else if(u<v){//if u is not included in the set of current frontiers
            if(!state->graph->terminals[u]&&n_hat->addedge[i1]==1&&x==0){
                //cout<<globalCount<<" added unnecessary edges for steriner trees"<<endl;
                delete node;
                bddnode.nonexistencePointer=-2;
                return NULL;
            }
            if(!state->graph->terminals[u]&&n_hat->addedge[i1]==0&&x==1){
                //cout<<globalCount<<" added unnecessary edges for steriner trees"<<endl;
                delete node;
                bddnode.existencePointer=-2;
                return NULL;
            }
            i1++;
        }
        if(i2>=frontierSize)break;

    }


    Edge edge = (*state->graph->GetEdgeList())[i - 1];

    if(x==0) { // edge is non-exist

        vector<int> id;
        vector<int> preid;
        vector<int> precompid;
        vector<int> predegree;
        vector<bool> preterminal;
        id.resize(2);
        preid.resize(2);

        precompid.resize(2);
        preterminal.resize(2);
        predegree.resize(2);

        id[0]=-1;
        id[1]=-1;
        preid[0]=-1;
        preid[1]=-1;

        precompid[0]=0;
        precompid[1]=0;
        preterminal[0]=false;
        preterminal[1]=false;
        predegree[0]=0;
        predegree[1]=0;

        for (int y = 0; y <= 1; ++y) {
            int u = (y == 0 ? edge.src : edge.dest);

            for (int j = 0; j < frontierSize; j++) {
                if (state->Fc[j] == u) {
                    id[y] = j;
                    break;
                }
            }

            if (id[y] == -1) {//id[y]==-1 -> u is not current Frontier

                for (int j = 0; j < previousFrontierSize; j++) {
                    if (state->Fp[j] == u) { // preid[y]==-1 -> u is not previous Frontier as well as current
                        preid[y]=j;
                        precompid[y]=n_hat->comp[j];
                        predegree[y]=n_hat->deg[j];

                        if(n_hat->terminal[j]>0)preterminal[y]=true;

                        for (int jj = 0; jj < frontierSize; jj++) {
                            if (node->comp[jj] == n_hat->comp[preid[y]]) {
                                node->deg[jj]--;
                            }
                        }
                        break;
                    }
                }
                if(preid[y]==-1){
                    if(state->graph->terminals[u]){
                        delete node;
                        bddnode.nonexistencePointer=-2;
                        return NULL;
                    }
                }
            }
            else {
                for (int j = 0; j < frontierSize; j++) {
                    if (node->comp[j] == node->comp[id[y]]) {
                        node->deg[j]--;
                    }
                }
            }
        }

        if(id[0]==-1&&preid[0]==-1&&preterminal[0]){//fridge
            delete node;
            bddnode.nonexistencePointer=-2;
            return NULL;
        }
        if(id[1]==-1&&preid[1]==-1&&preterminal[1]){//fridge
            delete node;
            bddnode.nonexistencePointer=-2;
            return NULL;
        }
        if(id[0]==-1&&id[1]==-1&&precompid[0]==precompid[1]&&preterminal[0]&&predegree[0]==2){
            delete node;
            bddnode.nonexistencePointer=-2;
            return NULL;
        }
        if(id[0]==-1&&predegree[0]==1&preterminal[0]){//separate component
            delete node;
            bddnode.nonexistencePointer=-2;
            return NULL;
        }
        if(id[1]==-1&&predegree[1]==1&&preterminal[1]){//separate component
            delete node;
            bddnode.nonexistencePointer=-2;
            return NULL;
        }

        for(int j=0;j<frontierSize;j++){
            if(node->terminal[j]>0&&node->deg[j]==0){
                delete node;
                bddnode.nonexistencePointer=-2;
                return NULL;
            }
        }

        bool fFlag=false;
        for(int j=0;j<frontierSize;j++){
            if(node->terminal[j]>0){fFlag=true;break;}
            if(node->terminal[j]>0&&node->deg[j]==0){
                assert(1);
            }
        }
        if(!fFlag){
            assert(1);
        }

    }
    else if (x == 1) // edge is exist
    {
        node->weight += edge.weight;
        if(node->weight >= weight_){
            delete node;
            bddnode.existencePointer=-2;
            return NULL;
        }

        int preTerminalTotalCount=0;
        for (int j = 0; j < previousFrontierSize; j++) {
            preTerminalTotalCount+=n_hat->terminal[j];
        }


        vector<int> id;
        vector<int> compid;
        vector<int> terminal;
        vector<int> degree;
        vector<int> preid;
        vector<int> precompid;
        vector<int> preterminal;
        vector<int> predegree;

        id.resize(2);
        compid.resize(2);
        terminal.resize(2);
        degree.resize(2);

        preid.resize(2);
        precompid.resize(2);
        preterminal.resize(2);
        predegree.resize(2);

        id[0]=-1;
        id[1]=-1;
        compid[0]=-1;
        compid[1]=-1;
        terminal[0]=0;
        terminal[1]=0;
        degree[0]=0;
        degree[1]=0;

        preid[0]=-1;
        preid[1]=-1;
        precompid[0]=-1;
        precompid[1]=-1;
        preterminal[0]=0;
        preterminal[1]=0;
        predegree[0]=0;
        predegree[1]=0;

        for (int y = 0; y <= 1; ++y) {
            int u = (y == 0 ? edge.src : edge.dest);

            for (int j = 0; j < frontierSize; j++) {
                if (state->Fc[j] == u) {
                    id[y] = j;
                    break;
                }
            }

            if (id[y] == -1) {

                for (int j = 0; j < previousFrontierSize; j++) {
                    if (state->Fp[j] == u) {
                        preid[y]=j;
                        precompid[y] = n_hat->comp[j];
                        predegree[y]=n_hat->deg[j];

                        if(n_hat->terminal[j]>0)preterminal[y]=n_hat->terminal[j];

                        for (int jj = 0; jj < frontierSize; jj++) {
                            if (node->comp[jj] == precompid[y]) {
                                node->deg[jj]--;
                            }
                        }
                        break;
                    }
                }
                if(preid[y]==-1){
                    predegree[y]=state->graph->degrees[u];
                    if(state->graph->terminals[u]){
                        preterminal[y]=1;
                    }
                }
            }
            else {
                for (int j = 0; j < frontierSize; j++) {
                    if (node->comp[j] == node->comp[id[y]]) {
                        node->deg[j]--;
                    }
                }
                compid[y]=node->comp[id[y]];
                terminal[y]=node->terminal[id[y]];
                degree[y]=node->deg[id[y]];
            }
        }

        if(id[0]!=-1)node->addedge[id[0]]++;
        if(id[1]!=-1)node->addedge[id[1]]++;


        int c_min = std::min(compid[0], compid[1]);
        int c_max = std::max(compid[0], compid[1]);

        if(precompid[0]==precompid[1]&&precompid[0]!=-1&&precompid[1]!=-1){//if they are previously in the same component, the component become having cycles
            delete node;
            bddnode.existencePointer=-2;
            return NULL;
        }
        if(compid[0]==precompid[1]&&compid[0]!=-1&&precompid[1]!=-1){//if they are previously in the same component, the component become having cycles
            delete node;
            bddnode.existencePointer=-2;
            return NULL;
        }
        if(precompid[0]==compid[1]&&precompid[0]!=-1&&compid[1]!=-1){//if they are previously in the same component, the component become having cycles
            delete node;
            bddnode.existencePointer=-2;
            return NULL;
        }
        if(compid[0]==compid[1]&&compid[0]!=-1&&compid[1]!=-1){//if they are previously in the same component, the component become having cycles
            delete node;
            bddnode.existencePointer=-2;
            return NULL;
        }

        if(id[0]==-1&&id[1]==-1){//both ends of edge are not frontiers
            if(precompid[0]!=precompid[1]&&preterminal[0]+preterminal[1]==state->terminals.size()){
                for(int j=0;j<frontierSize;j++) {
                    node->comp[j]=precompid[0];
                    node->terminal[j] = preterminal[0]+preterminal[1];
                }
                return node;//connected all terminals
            }

            if(precompid[0]!=precompid[1]&&(preterminal[0]>0||preterminal[1]>0)&&predegree[0]+predegree[1]==2){//if either of them are terminals and they are connected only by a bride
                delete node;
                bddnode.existencePointer=-2;
                return NULL;
            }

            if(precompid[0]==precompid[1]&&preterminal[0]>0&&predegree[0]==2){//if they are terminals, are included same component, and become isolated
                delete node;
                bddnode.existencePointer=-2;
                return NULL;
            }


        }



        if(id[0]==-1 || id[1]==-1){// handling not frontier vertices

            if(id[0]==-1 && preid[0]==-1&&id[1]!=-1){ // 0 is fridge
                for(int j=0;j<frontierSize;j++) {
                    if(node->comp[j]==node->comp[id[1]]) {
                        node->terminal[j] += preterminal[0];
                    }
                }
            }
            else if(id[1]==-1 && preid[1]==-1&&id[0]!=-1){ //1 is fridge
                for(int j=0;j<frontierSize;j++) {
                    if(node->comp[j]==node->comp[id[0]]) {
                        node->terminal[j] += preterminal[1];
                    }
                }
            }
            else if(id[0]==-1&&id[1]==-1&&preid[0]==-1&&preid[1]!=-1&&predegree[0]==1){//both 0  are fridges
                if(preterminal[0]>0) {
                    for (int j = 0; j < frontierSize; j++) {//compute degree
                        if (node->comp[j] == precompid[1]) {
                            node->terminal[j] += preterminal[0];
                            //cout<<node->terminal[j]<<","<<node->terminal[min_id]<<endl;
                        }
                    }
                }
            }
            else if(id[0]==-1&&id[1]==-1&&preid[0]!=-1&&preid[1]==-1&&predegree[1]==1){//for fridge
                if(preterminal[1]>0) {
                    for (int j = 0; j < frontierSize; j++) {//compute degree
                        if (node->comp[j] == precompid[0]) {
                            node->terminal[j] += preterminal[1];
                            //cout<<node->terminal[j]<<","<<node->terminal[min_id]<<endl;
                        }
                    }
                }
            }
            else if(id[0]==-1 && preid[0]!=-1){

                int c_maxDeg=0;
                int c_maxTer=0;

                int otherComp;
                int otherDeg;
                int otherTer;
                if(id[1]!=-1){
                    otherComp=compid[1];
                    otherDeg=degree[1];
                    otherTer=terminal[1];
                }
                else {
                    otherComp=precompid[1];
                    otherDeg=predegree[1]-1;
                    otherTer=preterminal[1];
                }

                c_min=min(otherComp,precompid[0]);
                c_max=max(otherComp,precompid[0]);

                if(c_min!=c_max) {
                    for (int j = 0; j < frontierSize; j++) {//compute degree
                        if (node->comp[j] == precompid[0]) {
                            //c_maxDeg = node->deg[j];
                            //c_maxTer = node->terminal[j];
                            node->comp[j] = otherComp;
                            node->deg[j] = otherDeg;
                            node->terminal[j] = otherTer;
                            //cout<<node->terminal[j]<<","<<node->terminal[min_id]<<endl;
                        }
                    }
                    for (int j = 0; j < frontierSize; j++) {
                        if (node->comp[j] == otherComp) {
                            node->deg[j] += predegree[0] - 1;
                            node->terminal[j] += preterminal[0];
                        }
                    }


                    for (int j = 0; j < frontierSize; j++) {
                        if (node->comp[j] == c_max) {
                            node->comp[j] = c_min;
                        }
                    }
                }

            }
            else if(id[1]==-1 && preid[1]!=-1){


                int c_maxDeg=0;
                int c_maxTer=0;

                int otherComp;
                int otherDeg;
                int otherTer;
                if(id[0]!=-1){
                    otherComp=compid[0];
                    otherDeg=degree[0];
                    otherTer=terminal[0];
                }
                else {
                    otherComp=precompid[0];
                    otherDeg=predegree[0]-1;
                    otherTer=preterminal[0];
                }

                c_min=min(otherComp,precompid[1]);
                c_max=max(otherComp,precompid[1]);

                if(c_min!=c_max) {
                    for (int j = 0; j < frontierSize; j++) {//compute degree
                        if (node->comp[j] == precompid[1]) {
                            c_maxDeg = node->deg[j];
                            c_maxTer = node->terminal[j];
                            node->comp[j] = otherComp;
                            node->deg[j] = otherDeg;
                            node->terminal[j] = otherTer;
                            //cout<<node->terminal[j]<<","<<node->terminal[min_id]<<endl;
                        }
                    }
                    for (int j = 0; j < frontierSize; j++) {
                        if (node->comp[j] == otherComp) {
                            node->deg[j] += predegree[1] - 1;
                            node->terminal[j] += preterminal[1];
                        }
                    }


                    for (int j = 0; j < frontierSize; j++) {
                        if (node->comp[j] == c_max) {
                            node->comp[j] = c_min;
                        }
                    }
                }

            }

        }
        else if(c_min!=c_max) { // both are frontiers

            int min_id;
            if(node->comp[id[0]]==c_min)min_id=id[0];
            else min_id = id[1];


            int c_maxDeg=0;
            int c_maxTer=0;


            for (int j = 0; j < frontierSize; j++) {//compute degree
                if (node->comp[j] == c_max) {
                    c_maxDeg = node->deg[j];
                    c_maxTer = node->terminal[j];
                    node->comp[j] = c_min;
                    node->deg[j] = node->deg[min_id];
                    node->terminal[j]=node->terminal[min_id];
                    //cout<<node->terminal[j]<<","<<node->terminal[min_id]<<endl;
                }
            }
            for (int j = 0; j < frontierSize; j++) {
                if (node->comp[j] == c_min) {
                    node->deg[j] += c_maxDeg;
                    node->terminal[j] += c_maxTer;
                }
            }
        }

    }


    if(i < state->graph->GetEdgeList()->size()) {
        //int fSize = state->F[i]->size();
        for (unsigned int j = 0; j < frontierSize; ++j) {
            int u = state->Fc[j];
            int comp_u = node->comp[j];
            //cout << u << "," << n_hat->comp[u] << std::endl;
            if (!Contains(state->Fc, comp_u)) {
                for (int jj = 0; jj < frontierSize; jj++) {
                    if (node->comp[jj] == comp_u)node->comp[jj] = u;//Frontier以外のcomp IDをFrontier NodeのIDにする。
                }
            }
        }
    }

    return node;


}

FrontierNode* FrontierAlgorithm::Find(FrontierNode* n_prime, const vector<FrontierNode*>& N_i, int i, State* state)
{
    for (unsigned int j = 0; j < N_i.size(); ++j) {
        FrontierNode* n_primeprime = N_i[j];

        if (IsEquivalent(n_prime, n_primeprime, i, state))
        {
            return n_primeprime;
        }
    }
    return NULL;
}


void FrontierAlgorithm::HashFrontier(FrontierNode* n_prime)
{
    return;
}


bool FrontierAlgorithm::IsEquivalent(FrontierNode* node1, FrontierNode* node2, int i, State* state)
{
    vector<int>& frontier = state->Fc;

    int fSize = frontier.size();
    for (unsigned int j = 0; j < fSize; ++j) {
        int v = frontier[j];
        if (node1->deg[v] != node2->deg[v]) {
            return false;
        }
        if (node1->comp[v] != node2->comp[v]) {
            return false;
        }
        if (node1->terminal[v] != node2->terminal[v]) {
            return false;
        }
    }
    return true;
}

