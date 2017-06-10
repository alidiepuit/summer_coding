//==================================================================
// File			: tsp.h
// Author		: rsagalyn
// Date			: Aug 18, 2013
// Description	:
//==================================================================
#ifndef MWM_H_
#define MWM_H_

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <pthread.h>
#include <queue>
#include <stack>
#include <string>
#include <stdio.h>
#include <vector>
#include <map>

#include "twoOpt.h"

using namespace std;

// Toggle printing debugging info to console
#define DEBUG 1

// Number of threads to use to fill N x N cost matrix
#define THREADS 1

// Calcualte lowest index controlled by thread id
#define START_AT(id,p,n) ((id)*(n)/(p))

// Calculate highest index controlled by thread id
#define END_AT(id,p,n) (START_AT((id)+1,p,n)-1)

//const direction
static const int DIRX[4] = { 0,-1, 0, 1};
static const int DIRY[4] = {-1, 0, 1, 0};

#define MAXVAL 10000
#define ll long long

class TSP
{
private:

	// x and y coords of a node
	struct City
	{
		int x;
		int y;
	};

	// Filename supplied on command line to read from
	string inFname;

	// Program-generated filename to output to
	string outFname;

	// List of odd nodes
	vector<long long>odds;

	// Smaller cost matrix used to store distances between odd nodes
	// Used to find minimum matching on odd nodes
	long long **cost;

	// Initialization function
	void getNodeCount();

	// Find odd vertices in graph
	void findOdds();

	// Prim helper function
	long long minKey(long long key[], bool mstSet[]);


protected:


public:
	// Number of cities
	long long n;

	// Number of gas stations
	int _numGasStation;

	//Gas tank size
	int _tankSize;

	//Map size
	int _row; //row
	int _col; //column

	//Map
	int **_originMap;


	long long **_graphId;

	vector<pair<int,int> > *_nearestGasStation;
	long long **_costGasStationToCity;
	long long **_graphIdGasStation;

	// euler circuit
	vector<long long>circuit;

	// Store cities and coords read in from file
	vector<City>cities;

	// Full n x n cost matrix of distances between each city
	long long **graph;

	// Current shortest path length
	long long pathLength;


	// Adjacency list
	// Array of n dynamic arrays, each holding a list of nodes it's index is attached to
	vector<long long> *adjlist;


	long long start_idx[THREADS];

	long long end_idx[THREADS];

	// n x 2 array to store length of TSP path starting at each node
	// col 0 => starting index   col 1 => path length from that node
	long long **path_vals;


	// Constructor
	TSP(string in, string out);

	// Destructor
	~TSP();

	// Calculate distance
	long long get_distance(struct City c1, struct City c2);


	// Initialization functions
	void readInput();
	void fillMatrix_threads();

	// Find MST using Prim's algorithm
	void findMST_old();


	// Find perfect matching
	void perfect_matching();

	// Find best node to start euler at
	// Doesn't create tour, just checks
	long long find_best_path(ll);

	// Create tour starting at specified node
	void create_tour(ll);

	// Private functions implemented by create_tour() and find_best_path()
	void euler (long long pos, vector<long long> &);
	//void euler(int);
	void make_hamilton(vector<long long> &, ll&);

	// Calls twoOpt function
	void make_shorter();


	// Debugging functions
	void printCities();
	void printAdjList();
	void printResult();
	void printEuler();
	void printPath();

	// Get node count
	long long get_size() {return n;};

	void initGraph();
	void floatMatrix(long long x, long long y, long long tankSize);
	bool isValidPosition(int, int);
	long long getMaxValue();
};

#endif /* MWM_H_ */
