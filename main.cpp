
#include "frontier.h"
#include "bdd.h"
#include "utility.h"
#include "state.h"
#include "graph.h"

using namespace std;


int main(int argc, char* argv[])
{
    string resultFile;
    string inputFile;

    bool graphReduce = false;
    int maxBDDNode = 0;//stoi(argv[2]);
    double errorTolerance = 0;// stod(argv[3]);

    int seedtreeNumber=0;
    int specifiedtreenum=1;
    double threshold = 0;
    vector<string> seedtreefiles;

    string samplingmethod="NONE";
    string ordering="BFS";

    try {
        for (int i = 1; i < argc; ++i) {
            std::string s = argv[i];
            if (s[0] == '-') {
                s = s.substr(1);

                if (s == "r") {
                    graphReduce = true;
                }
                else if (s == "i") {
                    inputFile = argv[++i];
                }
                else if (s == "o") {
                    resultFile = argv[++i];
                }
                else if (s == "m"){
                    maxBDDNode = stoi(argv[++i]);
                }
                else if(s=="k"){
                    specifiedtreenum=stoi(argv[++i]);
                }
                else if(s=="th"){
                    threshold=stod(argv[++i]);
                }

                else if(s=="s"){
                    seedtreeNumber=stoi(argv[++i]);
                    seedtreefiles.clear();
                    for(int j=0;j<seedtreeNumber;j++){
                        seedtreefiles.push_back(argv[++i]);
                    }
                }
                else {
                    throw std::exception();
                }
            }
            else {
                throw std::exception();
            }
        }
    }
    catch (std::exception& e) {
        cout<<"error input"<<endl;
        return 1;
    }

    Result result;
    result.outputFile = "./result/"+resultFile;
    result.inputFile = inputFile;
    result.maxwidth = maxBDDNode;
    result.seednum=seedtreeNumber;
    result.specifiedtreenum=specifiedtreenum;
    result.threshold=threshold;

    struct rlimit rl;
	int stacklimit = getrlimit(RLIMIT_STACK, &rl);
	rl.rlim_cur = RLIM_INFINITY;
	stacklimit = setrlimit(RLIMIT_STACK, &rl);

	struct rusage r1;
	if(getrusage(RUSAGE_SELF, &r1) != 0) {
		/*Failure*/
	}
	printf("maxrss=%ld\n", r1.ru_maxrss);


    struct timespec startTime, endTime, graphReducestartTime, graphReduceendTime;
	FrontierNode::Initialize();

	//BDDinputGraph graph;
    OriginalGraph originalgraph = OriginalGraph();
    OriginalGraph querygraph;


    originalgraph.InputFile(inputFile);
    result.optimal_answer = originalgraph.correct_answer;


    //////////////////////////////
    //Input Set of Steiner Trees//
    //////////////////////////////
    OriginalGraph output;
    SeedSteinerTreeSet seedtrees = SeedSteinerTreeSet();
    seedtrees.InputFile(&originalgraph,seedtreefiles);



    std::mt19937 terminalRand(1);

    double approximate_answer;
    double responseTime;


    querygraph=originalgraph;
    double graphPruneTime=0;
    double graphReduceTime=0;
    double graphSize=0;
    int graphSetSize=1;
    int maxwidth=0;
    approximate_answer=0;

    clock_gettime(CLOCK_REALTIME, &startTime);

    //////////////////////////////
    //Preprocess//
    //////////////////////////////

    BDD* bdd;
    State state;

    result.seedtreesize=seedtrees.edge_list.size();
    clock_gettime(CLOCK_REALTIME, &graphReducestartTime);
    if (graphReduce)seedtrees.GraphReduction();
    clock_gettime(CLOCK_REALTIME, &graphReduceendTime);
    graphReduceTime +=(graphReduceendTime.tv_sec - graphReducestartTime.tv_sec) + (graphReduceendTime.tv_nsec - graphReducestartTime.tv_nsec) * pow(10, -9);
    graphSize=seedtrees.edge_list.size();
    result.reduceseedtreesize=seedtrees.edge_list.size();

    //////////////////////////////
    //Construct BDD//
    //////////////////////////////

    clock_gettime(CLOCK_REALTIME, &startTime);
    state.SetState(&seedtrees, originalgraph.terminals, ordering);
    FrontierAlgorithm steinertreeenumeration;
    //pre_approximate_answer = 2*pre_approximate_answer;

    bdd = steinertreeenumeration.ComputeSteinertree(&state, maxBDDNode, errorTolerance, threshold);
    clock_gettime(CLOCK_REALTIME, &endTime);
    responseTime =
        (endTime.tv_sec - startTime.tv_sec) + (endTime.tv_nsec - startTime.tv_nsec) * pow(10, -9);
    result.bdd_construct_time = responseTime;


    //////////////////////////////
    //Traverse BDD//
    //////////////////////////////

    bdd->Reduce();
    clock_gettime(CLOCK_REALTIME, &startTime);
    //output = bdd->MinimumTraversal(&state);
    output = bdd->Traversal(&state, threshold, specifiedtreenum);

    clock_gettime(CLOCK_REALTIME, &endTime);
    responseTime = (endTime.tv_sec - startTime.tv_sec) + (endTime.tv_nsec - startTime.tv_nsec) * pow(10, -9);
    result.bdd_enumeration_time = responseTime;
    result.minimum_output_answers=bdd->minimumCost;
    result.average_output_answers=bdd->averageCost;
    result.enumeratetreenum =bdd->treenum;
    result.bddsize=bdd->GetSize();

    cout<<"BDDsize = "<<result.bddsize;



    std::cout << "response time:" << responseTime << std::endl;



    result.ResultOutput();
    if(getrusage(RUSAGE_SELF, &r1) != 0) {
		/*Failure*/
	}
	printf("maxrss=%ld GB \n", r1.ru_maxrss/(1024*1024));

//    result.maxRSS = r1.ru_maxrss;
	return 0;
}
