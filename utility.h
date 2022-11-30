
#ifndef FRONTIER_NETWORK_UTLITY_H
#define FRONTIER_NETWORK_UTLITY_H



#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <queue>
#include <stack>
#include <list>
#include <limits>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <omp.h>
#include <sys/resource.h>
#include <sparsehash/sparse_hash_map>
#include <google/dense_hash_map>


using namespace std;

typedef long long int int64; // 64 integer




//namespace mp = boost::multiprecision;
//typedef mp::number<mp::gmp_float<10000> > cpp_dec_float_100000;
//******************************************************************************

static bool Contains(const vector<int>& vec, int element)
{
	return std::find(vec.begin(), vec.end(), element) != vec.end();
}

static void Remove(vector<int>* vec, int element)
{
	vec->erase(std::remove(vec->begin(), vec->end(), element), vec->end());
}


static bool AllSameValues(int* comp, vector<int>& vec){

    for(unsigned int i=0; i < vec.size()-1;i++){
        if( comp[ vec[i] ] != comp[ vec[i+1] ])return false;
    }
    return true;
}

static bool ContainOneGivenValues(int* deg, vector<int>& vec, int _value){

    for(unsigned int i=0; i < vec.size();i++){
        if( deg[ vec[i] ] == _value )return true;
    }
    return false;
}

inline bool FindVector(std::vector<int> array, int object){

    for(unsigned int i=0; i < array.size();i++){
        if(array[i]==object)return true;
    }
    return false;
}

class Result {
public:


    string inputFile;
    string outputFile;

    int maxwidth;
    int seednum;
    int specifiedtreenum;
    double threshold;
    long double optimal_answer;
    long double average_output_answers;
    long double minimum_output_answers;
    int enumeratetreenum;
    int seedtreesize;
    int reduceseedtreesize;

    double bdd_construct_time;
    double bdd_enumeration_time;
    double total_time;

    int bddsize;


    void ResultOutput() {

        //const char *outputFile = str.c_str();
        std::ofstream fout(outputFile, ios::app);

        //double ratio_shortestbased_correct = shortestbased_answer/optimal_answer;
        //double ratio_bddshortestbased_corret= bdd_shortestbased_answer/optimal_answer;
        total_time = bdd_construct_time+bdd_enumeration_time;
        fout  << std::fixed<< inputFile<<"\t"<<maxwidth<<"\t"<<seednum<<"\t"<<specifiedtreenum<<"\t"<<threshold<<"\t" <<seedtreesize<<"\t" <<reduceseedtreesize<<"\t" << optimal_answer << "\t" << enumeratetreenum<<"\t" << average_output_answers << "\t" << minimum_output_answers <<"\t" << bdd_construct_time <<"\t" <<bdd_enumeration_time << "\t" <<total_time<<"\t"<<bddsize<<std::endl;

    }


};
//******************************************************************************

#endif //FRONTIER_NETWORK_UTLITY_H
