//==================================================================
// File			: main.cpp
// Author		: Rebecca Sagalyn
// Date			: Aug 25, 2013
// Description	: Driver for tsp.h
//==================================================================
#include <iostream>
#include <climits>
#include "tsp.h"
#include "usage.h"
#include "twoOpt.h"
#include "MyThread.h"		// thread wrapper class
// The length was annoying me.
#define CPS CLOCKS_PER_SEC

#define NUM_THREADS 1

int main(int argc, char** argv) {
	// Check that user entered filename on command line
	if (argc < 2)
	{
		usage();
		exit(-1);
	}

	// Read file names from input
	string f, o;
	f = o = argv[1];
	o.append(".tour");

	// Create new tsp object
	TSP tsp(f, o);

	// Start timing
	clock_t t2;

	// Read cities from file
	// if (DEBUG)
	// 	cout << "Reading cities" << endl;
	// tsp.readInput();
	// if (DEBUG)
	// 	cout << "Time to read cities: "
	// 			<< ((float) (clock() - t)) / CLOCKS_PER_SEC << " s\n";

	cout << "number of cities: " << tsp.n << endl;
	cout << "number of gas station: " << tsp._numGasStation << endl;
	// tsp.printCities();

	// Fill N x N matrix with distances between nodes
	if (DEBUG)
		cout << "\nFilling matrix" << endl;
	t2 = clock();
	tsp.fillMatrix_threads();
	tsp.initGraph();
	// tsp.initGraph();
	if (DEBUG)
		cout << "Time to fill matrix: " << ((float) (clock() - t2)) / CPS
				<< " s\n";

	

	// Find a MST T in graph G
	if (DEBUG)
		cout << "\nFinding mst" << endl;
	t2 = clock();
	tsp.findMST_old();
	if (DEBUG)
		cout << "Time to find mst: " << ((float) (clock() - t2)) / CPS
				<< " s\n";



	//print result to file
	tsp.printResult();

	//print detail path
	

	return 0;
}
