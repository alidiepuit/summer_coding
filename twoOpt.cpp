//==================================================================
// File			: twoOpt.cpp
// Author		: rsagalyn
// Date			: Aug 25, 2013
// Description	:
//==================================================================

#include "twoOpt.h"

// Input: edge 1's v, edge 2's u
// Remove edge 1 and edge 2, reconnect using new path
void reverse(vector<long long> &path, long long start, long long end, long long n)
{
	while (end - start > 0)
	{
		// Start, end is relative value to the start,
		// the index is start|slut % size
		long long temp = path[start % n];
		path[start % n] = path[end % n];
		path[end % n] = temp;
		start++;
		end--;
	}
}


long long is_path_shorter(long long **graph, long long v1, long long v2, long long v3, long long v4, long long &total_dist)
{
	if ((graph[v1][v3] + graph[v2][v4]) < (graph[v1][v2] + graph[v3][v4]))
	{
		total_dist -= (graph[v1][v2] + graph[v3][v4] - graph[v1][v3]
				- graph[v2][v4]);
		return 1;
	}
	return 0;
}


// Non-looping version of two-opt heuristic
long long twoOpt(long long **graph, vector<long long> &path, long long &pathLength, long long n)
{
	long long counter = 0;
	long long term_cond = 5;
	long long old_distance = pathLength;
	//long long size = path.size();
	long long v1, v2, u1, u2;

	// Iterate over each edge
	for (long long i = 0; i < n; i++)
	{
		// first edge
		u1 = i;
		v1 = (i+1) % n; // wrap around to first node if u1 is last node

		// Skip adjacent edges, start with node one past v1
		for (long long j = i + 2; (j + 1) % n != i; j++)
		{
			// mod by length to go back to beginning
			u2 = j % n;
			v2 = (j+1) % n;

			// Check if new edges would shorten the path length
			// If so, decreases pathLength
			if (is_path_shorter(graph, path[u1], path[v1], path[u2],
					path[v2], pathLength))
			{

				// Swap u1--v1 and u2--v2
				reverse(path, i + 1, j, n); // v1, u2

				if (pathLength - old_distance < 5 && counter == term_cond)
					break;

				// reset i
				if (pathLength - old_distance > term_cond )
					i = 0;

				old_distance = pathLength;
				counter++;
			}
		}
	}
	return pathLength;
}



long long get_path_length(long long **graph, vector<long long> &path, long long size){
	long long length = 0;
	for (long long i = 0; i < size-1; i++)
	{
		length += graph[path[i]][path[i+1]];
	}
	length += graph[path[size-1]] [path[0]]; // back home
	return length;
}


