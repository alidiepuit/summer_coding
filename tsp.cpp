#include "tsp.h"

struct thread_data {
	long long tid;
	TSP *tsp;
};
struct thread_data *data;

TSP::TSP(string in, string out){
	/////////////////////////////////////////////////////
	// Constructor
	/////////////////////////////////////////////////////
	inFname = in; outFname = out;
	// set n to number of lines read from input file
	readInput();

	// Allocate memory
	graph = new ll*[n];
	for (long long i = 0; i < n; i++) {
		graph[i] = new ll[n];
		for (long long j = 0; j < n; j++) {
			graph[i][j] = 0;
		}
	}

	cost = new ll*[n];
	for (long long i = 0; i < n; i++) {
		cost[i] = new ll[n];
	}

	path_vals = new ll*[n];
	for (long long i = 0; i < n; i++) {
		path_vals[i] = new ll[n];
	}

	// Adjacency lsit
	adjlist = new vector<long long> [n];
};

TSP::~TSP(){
	/////////////////////////////////////////////////////
	// Destructor
	/////////////////////////////////////////////////////

	for (long long i = 0; i < n; i++) {
		delete [] graph[i];
		delete [] _graphId[i];
		delete [] cost[i];
		delete [] path_vals[i];
	}
	delete [] path_vals;
	delete [] graph;
	delete [] _graphId;
	delete [] cost;
	delete [] adjlist;

	// for (long long i = 0; i < _row; i++) {
	// 	delete [] _originMap;
	// }
	// delete [] _originMap;
}

void TSP::readInput(){
	ifstream inStream;
	inStream.open(inFname.c_str(), ios::in);

	if (!inStream) {
	  cerr << "Can't open input file " << inFname << endl;
	  exit(1);
	}
	std::string unused;
	
	//read init position
	long long beg_x, beg_y;
	inStream >> beg_x >> beg_y;

	//read gas tank size
	inStream >> _tankSize;

	//read map size
	inStream >> _row >> _col;

	// Allocate memory
	_originMap = new ll*[_row];
	_graphId = new ll*[_row];
	for (long long i = 0; i < _row; i++) {
		_originMap[i] = new ll[_col];
		_graphId[i] = new ll[_col];
		for (long long j = 0; j < _col; j++) {
			_originMap[i][j] = 0;
			_graphId[i][j] = -1;
		}
	}

	long long count = 1;
	//first position
	struct City c = {beg_x-1, beg_y-1};
	cities.push_back(c);
	_graphId[beg_x-1][beg_y-1] = 0;
	for(long long i = 0; i < _row; i++)
		for(long long j = 0; j < _col; j++) {
			inStream >> _originMap[i][j];
			if (_originMap[i][j] == 3) {
				_graphId[i][j] = count;
				// Push back new city to vector
				struct City c = {i, j};
				cities.push_back(c);
				count++;
			}
			
		}

	n = count;
	inStream.close();
};

long long TSP::get_distance(struct TSP::City c1, struct TSP::City c2) {
	/////////////////////////////////////////////////////
	// Calculate distance between c1 and c2
	/////////////////////////////////////////////////////
	long long dx = pow((float)(c1.x - c2.x), 2);
	long long dy = pow((float)(c1.y - c2.y), 2);
	return (floor((float) (sqrt(dx + dy)) + 0.5));
};

void *F(void* args){
	struct thread_data *my_data = (struct thread_data *) args;
	long long tid = my_data->tid;
	TSP *tsp = my_data->tsp;
	long long **graph = tsp->graph;
	long long start, end;
	//start = START_AT(tid, THREADS, tsp->n);
	//end = END_AT(tid, THREADS, tsp->n);

	start = tsp->start_idx[tid];
	end = tsp->end_idx[tid];
	//cout << "thread " << setw(4) << left << tid << setw(8) << left << " start: " << setw(5) << left << start;
	//cout << setw(6) << left << " end: " << setw(5) << left << end << " load: " << end- start + 1 << endl;

	//clock_t t = clock();
	// fill matrix with distances from every city to every other city
	for (long long i = start; i <= end; i++) {
		for (long long j = i; j < tsp->n; j++) {
			// Don't delete this line  it's supposed to be there.
			// graph[i][j] = graph[j][i] =  tsp->get_distance(tsp->cities[i], tsp->cities[j]);
			graph[i][j] = graph[j][i] = tsp->getMaxValue();
		}
	}

	// cout << __LINE__ << "fuck" << endl;
	tsp->initGraph();

	//t = clock() - t;
	//t = clock();
	//cout << "thread " << tid << " time: " << 1000*(((float)clock())/CLOCKS_PER_SEC) << " s"<< endl;
	pthread_exit(NULL);
}

long long TSP::getMaxValue() {
	return 10000 < _tankSize ? _tankSize : 10000;
}

void TSP::initGraph() {
	for (vector<City>::iterator it = cities.begin(); it != cities.end(); ++it) {
		// cout << it->x << " " << it->y << endl;
		floatMatrix(it->x, it->y, _tankSize);
	}
}

void TSP::floatMatrix(long long x, long long y, long long tankSize) {
    queue<pair<int, City> > vec;
    struct City c = {x, y};
    vec.push(make_pair(0, c));
	map<int,bool> visitedGraph;
    long long idCity = _graphId[x][y];

    while (vec.size() > 0) {
        pair<int, City> p = vec.front();
        vec.pop();
        City city = p.second;
        long long maxV = p.first;
        // cout << maxV << " " << city.x << " " << city.y << " " << _graphId[city.x][city.y] << endl;
        //can't go
        if (maxV == tankSize) continue;

        visitedGraph[city.x*_row+city.y] = true;
        
        for(long long i = 0; i < 4; i++) {
            long long newx = city.x + DIRX[i];
            long long newy = city.y + DIRY[i];
            struct City c = {newx, newy};
            // cout << "check isValidPosition " << newx << " " << newy << " " << isValidPosition(newx, newy)  << " " << visitedGraph[newx*_row+newy] << endl;
            if (isValidPosition(newx, newy) && !visitedGraph[newx*_row+newy]) {
                visitedGraph[newx*_row+newy] = true;
                vec.push(make_pair(maxV+1,c));

                //check whether visit city
                if (_originMap[newx][newy] == 3) {
                	graph[idCity][_graphId[newx][newy]] = maxV+1;
                	graph[_graphId[newx][newy]][idCity] = maxV+1;
                }
            }
        }
    }

 //    cout << "done" << endl;
	// for(long long i = 0; i < n; i++) {
 //    	for(long long j = 0; j < n; j++) {
 //    		cout << graph[i][j] << " ";
 //    	}
 //    	cout << endl;
 //    }
    
}

bool TSP::isValidPosition(long long x, long long y) {
	return 0 <= x && x < _row && 0 <= y && y < _col && _originMap[x][y] != 0;
}

void TSP::fillMatrix_threads(){
	/////////////////////////////////////////////////////
	/////////////////////////////////////////////////////
	long long amount = (n / THREADS) * 0.2;
	long long x = (n / THREADS) - amount;		// min amount given to threads
	long long rem = n - (x * THREADS);
	long long half = THREADS/2 + 1;

	long long pos = 0;
	for (long long i = 0; i < half; i++) {
		start_idx[i] = pos;
		pos += (x - 1);
		end_idx[i] = pos;
		pos++;
	}
	long long remainingThreads = THREADS - half + 1;
	long long extraToEach = rem / remainingThreads;
	// Divide remainer among second half of threads
	for (long long i = half; i < THREADS; i++) {
		start_idx[i] = pos;
		pos += (x + extraToEach - 1);
		end_idx[i] = pos;
		pos++;
	}
	end_idx[THREADS-1] = n - 1;

	long long rc; void *status;
	data = new struct thread_data[n];

	// allocate space for n thread ids
	pthread_t *thread = new pthread_t[n];
	pthread_attr_t attr;

	// Initialize and set thread detached attribute
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	for (long t = 0; t < THREADS; t++) {
		//printf("Creating thread %ld\n", t);
		data[t].tid = t;
		data[t].tsp = this;
		rc = pthread_create(&thread[t], &attr, F, (void*)&data[t]);
		if (rc) {
			printf("ERROR; return code from pthread_create() is %llu\n", rc);
			exit(-1);
		}
	}

	// Free attribute and wait for the other threads
	pthread_attr_destroy(&attr);
	for (long t = 0; t < THREADS; t++) {
		rc = pthread_join(thread[t], &status);
		if (rc) {
			printf("ERROR; return code from pthread_join() is %llu\n", rc);
			exit(-1);
		}
		 //printf("Completed join with thread %ld having a status of %ld\n",t,(long)status);
	}
	delete [] data;
}

void TSP::findMST_old() {
	/////////////////////////////////////////////////////
	// In each iteration, we choose a minimum-weight
	// edge (u, v), connecting a vertex v in the set A to
	// the vertex u outside of set A
	/////////////////////////////////////////////////////
	long long key[n];   // Key values used to pick minimum weight edge in cut
	bool in_mst[n];  // To represent set of vertices not yet included in MST
	long long parent[n];

	// For each vertex v in V
	for (long long v = 0; v < n; v++) {
		// Initialize all keys to infinity
		key[v] = std::numeric_limits<long long>::max();

		// Mark as not being in mst yet
		in_mst[v] = false;
	}

	// Node 0 is the root node so give it the lowest distance (key)
	key[0] = 0;
	parent[0] = -1; // First node is always root of MST

	for (long long i = 0; i < n - 1; i++) {
		// Find closest remaining (not in tree) vertex
		// TO DO : This would be better represented by heap/pqueue
		long long v = minKey(key, in_mst);

		// Add vertex v to the MST
		in_mst[v] = true;

		// Look at each vertex u adjacent to v that's not yet in mst
		for (long long u = 0; u < n; u++) {
			if (graph[v][u] && in_mst[u] == false && graph[v][u] < key[u]) {
				// Update parent index of u
				parent[u] = v;

				// Update the key only if dist is smaller than key[u]
				key[u] = graph[v][u];
			}
		}
	}

	// map relations from parent array onto matrix
	for (long long v1 = 0; v1 < n; v1++) {
		// there is an edge between v1 and parent[v1]
		long long v2 = parent[v1];
		if (v2 != -1) {
			adjlist[v1].push_back(v2);
			adjlist[v2].push_back(v1);
		}
	}
};

// findMST helper function
long long TSP::minKey(long long key[], bool mstSet[]) {
	// Initialize min value
	long long min = std::numeric_limits<long long>::max();
	long long min_index;
	for (long long v = 0; v < n; v++)
		if (mstSet[v] == false && key[v] < min) {
			min = key[v];
			min_index = v;
		}
	return min_index;
};

void TSP::findOdds() {
	/////////////////////////////////////////////////////
	// Find nodes with odd degrees in T to get subgraph O
	/////////////////////////////////////////////////////

	// store odds in new vector for now
	for (long long r = 0; r < n; r++) {
		//cities[r].isOdd = ((adjlist[r].size() % 2) == 0) ? 0 : 1;
		if ((adjlist[r].size() % 2) != 0 ) {
			odds.push_back(r);
		}
	}
}

void TSP::perfect_matching() {
	/////////////////////////////////////////////////////
	// find a perfect matching M in the subgraph O using greedy algorithm
	// but not minimum
	/////////////////////////////////////////////////////
	long long closest, length; //long long d;
	std::vector<long long>::iterator tmp, first;

	// Find nodes with odd degrees in T to get subgraph O
	findOdds();

	// for each odd node
	while (!odds.empty()) {
		first = odds.begin();
		vector<long long>::iterator it = odds.begin() + 1;
		vector<long long>::iterator end = odds.end();
		length = std::numeric_limits<long long>::max();
		for (; it != end; ++it) {
			// if this node is closer than the current closest, update closest and length
			if (graph[*first][*it] < length) {
				length = graph[*first][*it];
				closest = *it;
				tmp = it;
			}
		}	// two nodes are matched, end of list reached
		adjlist[*first].push_back(closest);
		adjlist[closest].push_back(*first);
		odds.erase(tmp);
		odds.erase(first);
	}
}


// Take reference to a path vector
// so can either modify actual euler path or a copy of it
void TSP::euler (long long pos, vector<long long> &path) {
	/////////////////////////////////////////////////////////
	// Based on this algorithm:
	//	http://www.graph-magics.com/articles/euler.php
	// we know graph has 0 odd vertices, so start at any vertex
	// O(V+E) complexity
	/////////////////////////////////////////////////////////

	// make copy of original adjlist to use/modify
	vector<long long> *temp = new vector<long long> [n];
	for (long long i = 0; i < n; i++) {
		temp[i].resize(adjlist[i].size());
		temp[i] = adjlist[i];
	}

	path.clear();

	// Repeat until the current vertex has no more neighbors and the stack is empty.
	stack<long long> stk;
	while (!stk.empty() || temp[pos].size() > 0 ) {
		// If current vertex has no neighbors -
		if (temp[pos].size() == 0) {
			// add it to circuit,
			path.push_back(pos);
			// remove the last vertex from the stack and set it as the current one.
			long long last = stk.top();
			stk.pop();
			pos = last;
		}
		// Otherwise (in case it has neighbors)
		else {
			// add the vertex to the stack,
			stk.push(pos);
			// take any of its neighbors,
			long long neighbor = temp[pos].back();
			// remove the edge between selected neighbor and that vertex,
			temp[pos].pop_back();
	        for (unsigned long long i = 0; i < temp[neighbor].size(); i++)
	            if (temp[neighbor][i] == pos) { // find position of neighbor in list
	        	    temp[neighbor].erase (temp[neighbor].begin() + i); // remove it
	                break;
	            }
			// and set that neighbor as the current vertex.
	        pos = neighbor;
		}
	}
	path.push_back(pos);
}


void TSP::make_hamilton(vector<long long> &path, long long &path_dist) {
	// remove visited nodes from Euler tour
	bool visited[n]; // boolean value for each node if it has been visited yet
	memset(visited, 0, n * sizeof(bool));

	path_dist = 0;

	long long root = path.front();
	vector<long long>::iterator curr = path.begin();
	vector<long long>::iterator next = path.begin()+1;
	visited[root] = true;

	// loop until the end of the circuit list is reached
	while ( next != path.end() ) {
		// if we haven't been to the next city yet, go there
		if (!visited[*next]) {
			path_dist += graph[*curr][*next];
			curr = next;
			visited[*curr] = true;
			next = curr + 1;
		}else {
			next = path.erase(next); // remove it
		}
	}

	// add the distance back to the root
	path_dist += graph[*curr][*next];
}

void TSP::create_tour(long long pos){

	// call euler with actual circuit vector
	euler(pos, circuit);

	// make it hamiltonian
	// pass actual vars
	make_hamilton(circuit, pathLength);
}


// Does euler and hamilton but doesn't modify any variables
// Just finds path length from the node specified and returns it
long long TSP::find_best_path (long long pos) {

	// create new vector to pass to euler function
	vector<long long>path;
	euler(pos, path);

	// make it hamiltonian, pass copy of vars
	long long length;
	make_hamilton(path, length);

	// Optimize
	twoOpt(graph, path, length, n);
	twoOpt(graph, path, length, n);
	twoOpt(graph, path, length, n);
	twoOpt(graph, path, length, n);
	twoOpt(graph, path, length, n);

	return length;
}


void TSP::make_shorter(){
	// Modify circuit & pathLength
	twoOpt(graph, circuit, pathLength, n);
}




//================================ PRlong long FUNCTIONS ================================//

void TSP::printResult(){
	ofstream outputStream;
	outputStream.open(outFname.c_str(), ios::out);
	outputStream << pathLength-getMaxValue() << endl;
	for (vector<long long>::iterator it = circuit.begin(); it != circuit.end(); ++it) {
	//for (vector<long long>::iterator it = circuit.begin(); it != circuit.end()-1; ++it) {
		outputStream << *it << endl;
	}
	//outputStream << *(circuit.end()-1);
	outputStream.close();
};

void TSP::printPath(){
	cout << endl;
	for (vector<long long>::iterator it = circuit.begin(); it != circuit.end()-1; ++it) {
		cout << *it << " to " << *(it+1) << " ";
		cout << graph[*it][*(it+1)] << endl;
	}
	cout << *(circuit.end()-1) << " to " << circuit.front();

	cout << "\nLength: " << pathLength-getMaxValue() << endl << endl;
};

void TSP::printEuler() {
	for (vector<long long>::iterator it = circuit.begin(); it != circuit.end(); ++it)
		cout << *it << endl;
}

void TSP::printAdjList() {
	for (long long i = 0; i < n; i++) {
		cout << i << ": "; //prlong long which vertex's edge list follows
		for (long long j = 0; j < (int)adjlist[i].size(); j++) {
			cout << adjlist[i][j] << " "; //prlong long each item in edge list
		}
		cout << endl;
	}
};

void TSP::printCities(){
	cout << endl;
	long long i = 0;
	for (vector<City>::iterator it = cities.begin(); it != cities.end(); ++it)
		cout << i++ << ":  " << it->x << " " << it->y << endl;
}


