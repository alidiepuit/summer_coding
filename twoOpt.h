//==================================================================
// File			: twoOpt.h
// Author		: rsagalyn
// Date			: Aug 25, 2013
// Description		: Perform optimizations on graph
//==================================================================
#ifndef MYGRAPH_H
#define MYGRAPH_H

#include <vector>
#include <stack>
#include <iostream>
//#include <cstdlib>
//#include <cmath>

using namespace std;

// Non-looping version of two-opt optimization heuristic
long long twoOpt(long long **graph, vector<long long> &path, long long &pathLength, long long n);

// 2-Opt helper function: swap two nodes
long long is_path_shorter(long long **graph, long long v1, long long v2, long long v3, long long v4, long long &total_distance);

void reverse(vector<long long> &path, long long start, long long end, long long n);




// move this to tsp class
long long get_path_length(long long **graph, vector<long long> &path, long long size);


#endif
