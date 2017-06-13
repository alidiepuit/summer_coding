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
		// _nearestGasStation[i] = new vector<pair<int,int> >();
		for (long long j = 0; j < n; j++) {
			graph[i][j] = 0;
		}
	}

	_nearestGasStation = new vector<pair<int,int> > [n];
	_nearestCity = new vector<pair<int,int> > [_numGasStation];

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


	_costGasStationToCity = new ll*[_numGasStation];
	for (long long i = 0; i < _numGasStation; i++) {
		_costGasStationToCity[i] = new ll[n];
		for (long long j = 0; j < n; j++) {
			_costGasStationToCity[i][j] = 0;
		}
	}

	_graphGasStation = new ll*[_numGasStation];
	for (ll i = 0; i < _numGasStation; i++) {
		_graphGasStation[i] = new ll[_numGasStation];
		for (ll j = 0; j < _numGasStation; j++) {
			_graphGasStation[i][j] = 0;
		}
	}
};

TSP::~TSP(){
	/////////////////////////////////////////////////////
	// Destructor
	/////////////////////////////////////////////////////

	for (long long i = 0; i < n; i++) {
		delete [] graph[i];
		// delete [] _graphId[i];
		delete [] cost[i];
		delete [] path_vals[i];
	}
	delete [] path_vals;
	delete [] graph;
	// delete [] _graphId;
	delete [] cost;

	delete [] adjlist;

	for (long long i = 0; i < n; i++) {
		// delete [] _originMap[i];
	}
	// delete [] _originMap;

	for (long long i = 0; i < _numGasStation; i++) {
		// delete [] _costGasStationToCity[i];
	}
	// delete [] _costGasStationToCity;
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
	_positionStart = make_pair(beg_x-1, beg_y-1);

	//read gas tank size
	inStream >> _tankSize;

	//read map size
	inStream >> _row >> _col;

	// Allocate memory
	_originMap = new int*[_row];
	_graphId = new ll*[_row];
	_graphIdGasStation = new ll*[_row];
	for (long long i = 0; i < _row; i++) {
		_originMap[i] = new int[_col];
		_graphId[i] = new ll[_col];
		_graphIdGasStation[i] = new ll[_col];
		for (long long j = 0; j < _col; j++) {
			_originMap[i][j] = 0;
			_graphId[i][j] = -1;
			_graphIdGasStation[i][j] = -1;
		}
	}

	long long count = 1;
	long long countGasStation = 0;
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

			if (_originMap[i][j] == 2) {
				struct City c = {i, j};
				_listGasStation.push_back(c);
				_graphIdGasStation[i][j] = countGasStation;
				countGasStation++;
			}
			
		}


	n = count;
	_numGasStation = countGasStation;


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
	// long long **graph = tsp->graph;
	long long start, end;
	//start = START_AT(tid, THREADS, tsp->n);
	//end = END_AT(tid, THREADS, tsp->n);

	start = tsp->start_idx[tid];
	end = tsp->end_idx[tid];
	//cout << "thread " << setw(4) << left << tid << setw(8) << left << " start: " << setw(5) << left << start;
	//cout << setw(6) << left << " end: " << setw(5) << left << end << " load: " << end- start + 1 << endl;

	tsp->initGraphGasStation(start, end);

	//clock_t t = clock();
	// fill matrix with distances from every city to every other city
	// for (long long i = start; i <= end; i++) {
	// 	for (long long j = i; j < tsp->n; j++) {
			// Don't delete this line  it's supposed to be there.
			// graph[i][j] = graph[j][i] =  tsp->get_distance(tsp->cities[i], tsp->cities[j]);
	// 	}
	// }

 //    cout << "print graph" << endl;
	// for(long long i = 0; i < tsp->_numGasStation; i++) {
 //    	for(long long j = 0; j < tsp->_numGasStation; j++) {
 //    		cout << tsp->_graphGasStation[i][j] << " ";
 //    	}
 //    	cout << endl;
 //    }

	//t = clock() - t;
	//t = clock();
	//cout << "thread " << tid << " time: " << 1000*(((float)clock())/CLOCKS_PER_SEC) << " s"<< endl;
	pthread_exit(NULL);
}

void TSP::initGraph(ll start, ll end) {
	for (ll i = start; i <= end; i++) {
		City city = cities[i];
		// cout << "begin " << city.x << " " << city.y << endl;
		floatMatrix(city.x, city.y, _tankSize);
	}

}

void TSP::floatMatrix(int x, int y, long long tankSize) {
    queue<pair<int, City> > vec;
    struct City c = {x, y};
    vec.push(make_pair(0, c));
	int visitedGraph[_row+1][_col+1];
	memset(visitedGraph, -1, (_row+1) * (_col+1) * sizeof(int));
    long long idCity = _graphId[x][y];

    long long maxV = 0;
    // cout << "fuck " << c.x << " " << c.y << endl;
    visitedGraph[c.x][c.y] = maxV;
    while (vec.size() > 0) {
        pair<ll, City> p = vec.front();
        vec.pop();
        City city = p.second;
        maxV = p.first;
        //can't go
        // cout << "visit " << city.x << " " << city.y << endl;
        // if (maxV == tankSize) {
        	// cout << "read quota " << maxV << " " << tankSize << endl;
        	// continue;
		// }
        
        for(int i = 0; i < 4; i++) {
            int newx = city.x + DIRX[i];
            int newy = city.y + DIRY[i];
            struct City c = {newx, newy};
            // cout << "check isValidPosition " << newx << " " << newy << " " << isValidPosition(newx, newy)  << " " << visitedGraph[newx][newy] << endl;
            if (isValidPosition(newx, newy) && visitedGraph[newx][newy]==-1) {
                visitedGraph[newx][newy] = maxV+1;
                vec.push(make_pair(maxV+1,c));

                //check whether visit city
                if (_originMap[newx][newy] == 3) {
                	graph[idCity][_graphId[newx][newy]] = maxV+1;
                	graph[_graphId[newx][newy]][idCity] = maxV+1;
                }

                //check whether visit gasstation
                if (_originMap[newx][newy] == 2) {
                	long long idGasStation = _graphIdGasStation[newx][newy];
                	_costGasStationToCity[idGasStation][idCity] = maxV+1;
                	_nearestGasStation[idCity].push_back(make_pair(newx,newy));
                }
            }
        }
    }
  //   if (idCity == 1) {
	 //    cout << "done" << endl;
		// for(long long i = 0; i < _row; i++) {
	 //    	for(long long j = 0; j < _col; j++) {
	 //    		cout << visitedGraph[i][j] << " ";
	 //    	}
	 //    	cout << endl;
	 //    }

	 //    cout << "_nearestGasStation to city " << idCity << endl;
	 //    for(vector<pair<int,int> >::iterator it = _nearestGasStation[idCity].begin(); it != _nearestGasStation[idCity].end(); it++) {
	 //    	long long idGasStation = _graphIdGasStation[(*it).first][(*it).second];
	 //    	cout << (*it).first << " " << (*it).second << " = " << _costGasStationToCity[idGasStation][idCity] << endl;
	 //    }
  //   }
}

void TSP::initGraphGasStation(ll start, ll end) {
	// Allocate memory
	
	for (ll i = start; i <= end; i++) {
		City city = _listGasStation[i];
		// cout << "begin " << city.x << " " << city.y << endl;
		floatMatrixGasStation(city.x, city.y, _tankSize);
	}
}

void TSP::floatMatrixGasStation(int x, int y, ll tankSize) {
    queue<pair<int, City> > vec;
    struct City c = {x, y};
    vec.push(make_pair(0, c));
	int visitedGraph[_row+1][_col+1];
	memset(visitedGraph, -1, (_row+1) * (_col+1) * sizeof(int));
    long long idGasStation = _graphIdGasStation[x][y];

    long long maxV = 0;
    // cout << "fuck " << c.x << " " << c.y << endl;
    visitedGraph[c.x][c.y] = maxV;
    while (vec.size() > 0) {
        pair<ll, City> p = vec.front();
        vec.pop();
        City city = p.second;
        maxV = p.first;

        if (maxV == _tankSize) {
        	return;
        }
        
        for(int i = 0; i < 4; i++) {
            int newx = city.x + DIRX[i];
            int newy = city.y + DIRY[i];
            struct City c = {newx, newy};
            // cout << "check isValidPosition " << newx << " " << newy << " " << isValidPosition(newx, newy)  << " " << visitedGraph[newx][newy] << endl;
            if (isValidPosition(newx, newy) && visitedGraph[newx][newy]==-1 && maxV+1 <= _tankSize) {
                visitedGraph[newx][newy] = maxV+1;
                vec.push(make_pair(maxV+1,c));

                //check whether visit gas station
                if (_originMap[newx][newy] == 2) {
                	_graphGasStation[idGasStation][_graphIdGasStation[newx][newy]] = maxV+1;
                	_graphGasStation[_graphIdGasStation[newx][newy]][idGasStation] = maxV+1;
                }

                //check whether visit city
                if (_originMap[newx][newy] == 3 || (newx == _positionStart.first && newy == _positionStart.second)) {
                	long long idCity = _graphId[newx][newy];
                	_costGasStationToCity[idGasStation][idCity] = maxV+1;
                	_nearestCity[idGasStation].push_back(make_pair(newx,newy));
                	// cout << idGasStation << " " << idCity << " " << maxV+1 << endl;
                	if (maxV+1 <= tankSize / 2 || (newx == _positionStart.first && newy == _positionStart.second)) {
                		_nearestGasStation[idCity].push_back(make_pair(x,y));
            		}
                }
            }
        }
    }
}

bool TSP::isValidPosition(int x, int y) {
	return 0 <= x && x < _row && 0 <= y && y < _col && _originMap[x][y] != 0;
}

void TSP::fillMatrix_threads(){
	/////////////////////////////////////////////////////
	/////////////////////////////////////////////////////
	long long amount = (_numGasStation / THREADS) * 0.2;
	long long x = (_numGasStation / THREADS) - amount;		// min amount given to threads
	long long rem = _numGasStation - (x * THREADS);
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
	end_idx[THREADS-1] = _numGasStation - 1;

	long long rc; void *status;
	data = new struct thread_data[_numGasStation];

	// allocate space for _numGasStation thread ids
	pthread_t *thread = new pthread_t[_numGasStation];
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
	long long key[_numGasStation];   // Key values used to pick minimum weight edge in cut
	bool in_mst[_numGasStation];  // To represent set of vertices not yet included in MST
	long long parent[_numGasStation];
	memset(parent, 0, _numGasStation * sizeof(ll));

	// For each vertex v in V
	for (long long v = 0; v < _numGasStation; v++) {
		// Initialize all keys to infinity
		key[v] = std::numeric_limits<long long>::max();

		// Mark as not being in mst yet
		in_mst[v] = false;
	}

	// cout << __LINE__ << " fuck " << endl;
	// for(vector<pair<int,int> >::iterator it = _nearestGasStation[0].begin(); it != _nearestGasStation[0].end(); it++) {
	// 	cout << (*it).first << " " << (*it).second << endl;
	// }
	pair<int,int> firstGasStation = *(_nearestGasStation[0].begin());
	ll idFirstGasStation = _graphIdGasStation[firstGasStation.first][firstGasStation.second];
	// cout << "fuck " << idFirstGasStation << endl;

	// Node 0 is the root node so give it the lowest distance (key)
	key[idFirstGasStation] = 0;
	parent[idFirstGasStation] = -1; // First node is always root of MST

	for (long long i = 0; i < _numGasStation - 1; i++) {
		// Find closest remaining (not in tree) vertex
		// TO DO : This would be better represented by heap/pqueue
		long long v = minKey(key, in_mst);

		// Add vertex v to the MST
		in_mst[v] = true;

		// cout << "find v " << v << endl;

		// Look at each vertex u adjacent to v that's not yet in mst
		for (long long u = 0; u < _numGasStation; u++) {
			if (_graphGasStation[v][u] > 0 && in_mst[u] == false && _graphGasStation[v][u] < key[u]) {
				// Update parent index of u
				parent[u] = v;
				// cout << "find u " << u << " " << v << " " << _graphGasStation[v][u] << endl;
				// Update the key only if dist is smaller than key[u]
				key[u] = _graphGasStation[v][u];
			}
		}
	}

	// cout << n << endl;
	// for (long long v1 = 0; v1 < _numGasStation; v1++) cout << "fuck " << v1 << " " << parent[v1] << endl;

	// map relations from parent array onto matrix
	for (long long v1 = 0; v1 < _numGasStation; v1++) {
		// there is an edge between v1 and parent[v1]
		long long v2 = parent[v1];
		if (v2 != -1) {
			// cout << __LINE__ << " fuck " << v1 << " " << v2 << endl;
			adjlist[v1].push_back(v2);
			adjlist[v2].push_back(v1);
		}
	}

	bool *mark = new bool[_numGasStation];
	memset(mark, 0, _numGasStation * sizeof(bool));
	ll t = _numGasStation-1;
	for(int i = _numGasStation-1; i >= 0; i--)
		if (!mark[i]) {
			vector<ll> tmp;
			t = i;
			while (t != -1) {
				tmp.push_back(t);
				// cout << t << " ";
				mark[t] = 1;
				t = parent[t];
			}
			// cout << "======" << endl;
			std::reverse(tmp.begin(),tmp.end());
			_connectedGasStation.push_back(tmp);
		}

	memset(mark, 0, _numGasStation * sizeof(bool));

	bool *markCity = new bool[n];
	memset(markCity, 0, n * sizeof(bool));
	_finalPath.push_back(_positionStart);
	markCity[0] = 1;
	ll idGasStation = 0;
	bool hasCity = false;
	for(vector< vector<ll> >::iterator it = _connectedGasStation.begin(); it != _connectedGasStation.end(); it++) {
		vector<ll> connected = *it;
		hasCity = false;
		vector<pair<int,int> > tmpPath;
		vector<ll> tempGasStation;
		cout << "temp gas station: ";
		for(vector<ll>::iterator gasStation = connected.begin(); gasStation != connected.end(); gasStation++) {
			idGasStation =  *gasStation;
			cout << idGasStation << " ";
			tempGasStation.push_back(idGasStation);
				// cout << "visited " << idGasStation << endl;
				City coorGasStation = _listGasStation[idGasStation];
				tmpPath.push_back(make_pair(coorGasStation.x, coorGasStation.y));
				for(vector<pair<int,int> >::iterator city = _nearestCity[idGasStation].begin(); city != _nearestCity[idGasStation].end(); city++) {
					ll idCity = _graphId[(*city).first][(*city).second];
					if (_costGasStationToCity[idGasStation][idCity] <= _tankSize/2 && !markCity[idCity]) {
						tmpPath.push_back(*city);
						tmpPath.push_back(make_pair(coorGasStation.x, coorGasStation.y));
						markCity[idCity] = true;
						hasCity = true;
					}
				}
				mark[idGasStation] = true;
		}
		cout << endl;
		if (hasCity) {
			cout << "has city" << endl;
			for(vector<pair<int,int> >::iterator itTmp = tmpPath.begin(); itTmp != tmpPath.end(); itTmp++) {
				_finalPath.push_back(*itTmp);
			}
			reverse(tempGasStation.begin(),tempGasStation.end());
			cout << "back to root ";
			for(vector<ll>::iterator i = tempGasStation.begin(); i != tempGasStation.end(); i++) {
				cout << *i << " ";
				ll idGasStation =  *i;
				City coorGasStation = _listGasStation[idGasStation];
				_finalPath.push_back(make_pair(coorGasStation.x, coorGasStation.y));
			}
			cout << endl;
		}
	}
	for(vector<pair<int,int> >::iterator rit = _finalPath.begin(); rit != _finalPath.end(); rit++) {
		pair<int,int> coor = *rit;
		if (_graphIdGasStation[coor.first][coor.second] > 0) {
			ll idGasStation = _graphIdGasStation[coor.first][coor.second];
			for(vector<pair<int,int> >::iterator city = _nearestCity[idGasStation].begin(); city != _nearestCity[idGasStation].end(); city++) {
				ll idCity = _graphId[(*city).first][(*city).second];
				if (!markCity[idCity]) {
					_finalPath.push_back(*city);
					stop = true;
					break;
				}
			}
			if (stop) break;
			_finalPath.push_back(make_pair(coor.first,coor.second));
		}
	}
};

// findMST helper function
long long TSP::minKey(long long key[], bool mstSet[]) {
	// Initialize min value
	long long min = std::numeric_limits<int>::max();
	long long min_index;
	for (long long v = 0; v < _numGasStation; v++)
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
	for (long long r = 0; r < _numGasStation; r++) {
		//cities[r].isOdd = ((adjlist[r].size() % 2) == 0) ? 0 : 1;
		// cout << r << " " << adjlist[r].size() << endl;
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

	// for (vector<long long>::iterator it = odds.begin(); it != odds.end(); ++it) {
	// 	cout << *it << endl;
	// }

	// for each odd node
	while (!odds.empty()) {
		first = odds.begin();
		vector<long long>::iterator it = odds.begin() + 1;
		vector<long long>::iterator end = odds.end();
		length = std::numeric_limits<int>::max();
		closest = -1;
		tmp = first;
		for (; it != end; ++it) {
			// if this node is closer than the current closest, update closest and length
			if (_graphGasStation[*first][*it] > 0 && _graphGasStation[*first][*it] < length) {
				length = _graphGasStation[*first][*it];
				closest = *it;
				tmp = it;
			}
		}	// two nodes are matched, end of list reached

		if (closest >= 0) {
			adjlist[*first].push_back(closest);
			adjlist[closest].push_back(*first);
		}
		if (tmp != first) {
			odds.erase(tmp);
		}
		odds.erase(first);
	}
}


// Take reference to a path vector
// so can either modify actual euler path or a copy of it
void TSP::euler(long long pos, vector<long long> &path) {
	/////////////////////////////////////////////////////////
	// Based on this algorithm:
	//	http://www.graph-magics.com/articles/euler.php
	// we know graph has 0 odd vertices, so start at any vertex
	// O(V+E) complexity
	/////////////////////////////////////////////////////////

	// make copy of original adjlist to use/modify
	vector<long long> *temp = new vector<long long> [_numGasStation];
	for (long long i = 0; i < _numGasStation; i++) {
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
			// cout << "vertext has no neighbors " << pos << endl;
			// remove the last vertex from the stack and set it as the current one.
			long long last = stk.top();
			stk.pop();
			pos = last;
		}
		// Otherwise (in case it has neighbors)
		else {
			// add the vertex to the stack,
			// cout << "Otherwise " << pos << endl;
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
	bool visited[_numGasStation]; // boolean value for each node if it has been visited yet
	memset(visited, 0, _numGasStation * sizeof(bool));

	path_dist = 0;

	long long root = path.front();
	vector<long long>::iterator curr = path.begin();
	vector<long long>::iterator next = path.begin()+1;
	visited[root] = true;

	// loop until the end of the circuit list is reached
	while ( next != path.end() ) {
		// if we haven't been to the next city yet, go there
		if (!visited[*next]) {
			// cout << *curr << " " << *next << " " << graph[*curr][*next] << endl;
			path_dist += _graphGasStation[*curr][*next];
			curr = next;
			visited[*curr] = true;
			next = curr + 1;
		}else {
			next = path.erase(next); // remove it
		}
	}

	// for summer coding, we don't need go back to root
	// add the distance back to the root
	// path_dist += graph[*curr][*next];
}

void TSP::create_tour(long long pos){

	// call euler with actual circuit vector
	// euler(pos, circuit);

	// printEuler();

	// make it hamiltonian
	// pass actual vars
	// make_hamilton(circuit, pathLength);
}


// Does euler and hamilton but doesn't modify any variables
// Just finds path length from the node specified and returns it
long long TSP::find_best_path(long long pos) {
	vector<pair<int,int> > tempPath;
	long long length = 0;
	// ll gasRemain = _tankSize;
	_fixPath = new ll[circuit.size()];
	ll sizePath = circuit.size();
	ll i = 0;
	for (vector<ll>::iterator it = circuit.begin(); it != circuit.end(); ++it) {
		_fixPath[i] = *it;
		i++;
	}

	stop = false;
	maxNumCitiesFound = 0;
	try_to_find_optimal_solution(tempPath, _tankSize, 0, sizePath, 0);

	// cout << __LINE__ << " done" << endl;
	// if (_finalPath.size() > 0) {
	// 	for (vector<pair<int,int> >::iterator it = _finalPath.begin(); it != _finalPath.end(); ++it) {
	// 		cout << (*it).first << " " << (*it).second << endl;
	// 	}
	// } else {
		// cout << __LINE__ << " fuck" << endl;
	// }

	return length;
}

//trace circuit path
//find gas station fit maximum number of heart
//length from gas station to heart <= gas tank size / 2
void TSP::find_best_path_2() {
	_fixPath = new ll[circuit.size()];
	ll sizePath = circuit.size();
	ll i = 0;
	for (vector<ll>::iterator it = circuit.begin(); it != circuit.end(); ++it) {
		_fixPath[i] = *it;
		i++;
	}

	ll pos = 0;
	ll gasRemain = _tankSize;
	pair<int,int> prevGasStation;
	while (pos < sizePath) {
		//DEBUG
		cout << "============" << endl;
		cout << "first index path " << pos << endl;

		//pair<path vector<ll>, gas station pair<int,int> >
		pair<vector<ll>, pair<int,int> > solution = try_to_find_maximum_heart_has_same_gas_station(pos, sizePath);
		vector<ll> listCity = solution.first;
		pair<int,int> gasStation = solution.second;

		//DEBUG
		cout << "solution: ";
		for(vector<ll>::iterator it = listCity.begin(); it != listCity.end(); it++) {
			cout << *it << " ";
		}
		cout << endl;

		//put indexCurrentCity to backlog
		if (listCity.empty()) {
			//TODO
			// ll indexCurrentCity = _fixPath[pos];
			pos++;
		} else {
			//first group
			if (!_finalPath.empty()) {
				pair<int,int> lastCity = *(_finalPath.end()-1);
				ll idLastCity = _graphId[lastCity.first][lastCity.second];
				ll idFirstCityInNextGroup = *(listCity.begin());
				ll costLastCityToNextGroup = graph[idLastCity][idFirstCityInNextGroup];
				ll idGasStation = _graphIdGasStation[gasStation.first][gasStation.second];
				ll costGasStationToNextCity = _costGasStationToCity[idGasStation][idFirstCityInNextGroup];

				//--found list cities
				cout << "gas remain " << gasRemain << endl;
				cout << "costLastCityToNextGroup " << idLastCity << " " << idFirstCityInNextGroup << " " << costLastCityToNextGroup << endl;
				cout << "costGasStationToNextCity " << costGasStationToNextCity << endl;

				//not direct path
				if (gasRemain - costLastCityToNextGroup - costGasStationToNextCity < 0) {
					_finalPath.push_back(prevGasStation);
					gasRemain = _tankSize - costGasStationToNextCity;
				} else {
					gasRemain -= costLastCityToNextGroup;
				}
			}

			//find group single gas station
			pair<vector<pair<int,int> >, ll >  solution= greedy_single_gas_station(listCity, gasStation, gasRemain);
			vector<pair<int,int> > pathVisitGroup = solution.first;
			//DEBUG
			cout << "path visit: ";
			for(vector<pair<int,int> >::iterator it = pathVisitGroup.begin(); it != pathVisitGroup.end(); it++) {
				cout << "(" << (*it).first << "," << (*it).second << ")" << " ";
			}
			cout << endl;

			for(vector<pair<int,int> >::iterator it = pathVisitGroup.begin(); it != pathVisitGroup.end(); it++) {
				_finalPath.push_back(*it);
			}

			pos += listCity.size();
			gasRemain = solution.second;
			prevGasStation = gasStation;
		}
	}

	cout << "done" << endl;
	if (_finalPath.size() > 0) {
		for (vector<pair<int,int> >::iterator it = _finalPath.begin(); it != _finalPath.end(); ++it) {
			cout << (*it).first << " " << (*it).second << endl;
		}
	} else {
		cout << __LINE__ << " fuck" << endl;
	}
}//end: find_best_path_2

//pair<path, indexGasStation>
pair<vector<ll>, pair<int,int> > TSP::try_to_find_maximum_heart_has_same_gas_station(ll pos, ll sizePath) {
	set<ll> setGasStation;
	vector<ll> res;
	cout << "try_to_find_maximum_heart_has_same_gas_station" << endl;
	while (pos < sizePath) {
		int idCurrentCity = _fixPath[pos];
		cout << "current city " << idCurrentCity << endl;
		set<ll> tempGasStation;
		//find all gas station nearest current city with cost <= _tankSize / 2
		for (vector<pair<int,int> >::iterator it = _nearestGasStation[idCurrentCity].begin(); it != _nearestGasStation[idCurrentCity].end(); ++it) {
			ll idGasStation = _graphIdGasStation[(*it).first][(*it).second];
			if (_costGasStationToCity[idGasStation][idCurrentCity] <= _tankSize / 2)
				tempGasStation.insert(idGasStation);
		}

		cout << "nearest gas station: ";
		for (set<ll>::iterator it = tempGasStation.begin(); it != tempGasStation.end(); ++it) {
			cout << *it << " ";
		}
		cout << endl;

		if (setGasStation.empty()) {
			setGasStation = tempGasStation;
		}

		std::vector<ll> v_intersection;
		std::set_intersection(tempGasStation.begin(), tempGasStation.end(),
                          setGasStation.begin(), setGasStation.end(),
                          std::back_inserter(v_intersection));
		if (v_intersection.empty() || tempGasStation.empty()) {
			pair<int,int> gasStation = make_pair(-1,-1);
			if (!setGasStation.empty()) {
				City city = _listGasStation[*(setGasStation.begin())];
				gasStation = make_pair(city.x,city.y);
			}
			return make_pair(res, gasStation);
		}

		res.push_back(idCurrentCity);

		setGasStation.clear();
		for (vector<ll>::iterator it = v_intersection.begin(); it != v_intersection.end(); ++it) {
			setGasStation.insert(*it);
		}

		//next city
		pos++;
	}

	//visited gasStation 
	City city = _listGasStation[*(setGasStation.begin())];
	pair<int,int> gasStation = make_pair(city.x,city.y);
	return make_pair(res, gasStation);
} //end: try_to_find_maximum_heart_has_same_gas_station

//pair<path, gas remain>
pair<vector<pair<int,int> >, ll> TSP::greedy_single_gas_station(vector<ll> pathCities, pair<int,int> gasStation, ll gasRemain) {
	ll i = 0;
	ll length = pathCities.size();
	vector<pair<int,int> > res;
	while (i < length) {
		ll indexCurrentCity = pathCities[i];
		City currentCity = cities[indexCurrentCity];

		res.push_back(make_pair(currentCity.x,currentCity.y));

		//visited last city
		if (i == length-1) {
			break;
		}

		ll indexNextCity = pathCities[i+1];
		ll indexGasStation = _graphIdGasStation[gasStation.first][gasStation.second];
		ll costCurrentCityToNextCity = graph[indexCurrentCity][indexNextCity];
		ll costGasStationToNextCity = _costGasStationToCity[indexGasStation][indexNextCity];

		cout << "cost current city to next city " << indexCurrentCity << " " << indexNextCity << " " << costCurrentCityToNextCity << endl;
		cout << "cost gas station to next city " << costGasStationToNextCity << endl;

		if (gasRemain - costCurrentCityToNextCity - costGasStationToNextCity >= 0) {
			gasRemain -= costCurrentCityToNextCity;
		} else {
			res.push_back(make_pair(gasStation.first,gasStation.second));
			gasRemain = _tankSize - costGasStationToNextCity;
		}

		//visit next city
		i++;
	}

	//insert gas station to path
	// res.push_back(make_pair(gasStation.first,gasStation.second));
	return make_pair(res, gasRemain);
}

void TSP::try_to_find_optimal_solution(vector<pair<int,int> > tempPath, ll gasRemain, ll pos, ll sizePath, ll numCities) {
	if (stop) return;
	City city = cities[_fixPath[pos]];
	// cout << "======================" << endl;
	// cout << "city " << city.x << " " << city.y << endl;
	// cout << "gas remain " << gasRemain << endl;
	tempPath.push_back(make_pair(city.x,city.y));
	if (numCities > maxNumCitiesFound) {
		// cout << "haha" << endl;
		// for (vector<pair<int,int> >::iterator it = tempPath.begin(); it != tempPath.end(); ++it) {
		// 	cout << (*it).first << " " << (*it).second << endl;
		// }
		_finalPath = tempPath;
		maxNumCitiesFound = numCities;
	}
	if (pos == sizePath-1) {
		// cout << "done" << endl;
		// for (vector<pair<int,int> >::iterator it = tempPath.begin(); it != tempPath.end(); ++it) {
		// 	cout << (*it).first << " " << (*it).second << endl;
		// }
		// stop = true;
		return;
	} else {
		for(ll i = pos+1; i < sizePath; i++) {
			ll idCurrentCity = _fixPath[pos];
			ll idNextCity = _fixPath[i];

			if (graph[idCurrentCity][idNextCity] > 0 && gasRemain - graph[idCurrentCity][idNextCity] >= 0) {
				// cout << idCurrentCity << " go directly to city " << idNextCity << " = " << graph[idCurrentCity][idNextCity] << endl;
				try_to_find_optimal_solution(tempPath, gasRemain - graph[idCurrentCity][idNextCity], i, sizePath, numCities+1);
			}

			// cout << idCurrentCity << " can not go directly to city " << idNextCity << " = " << graph[idCurrentCity][idNextCity] << endl;

			pair<int,int> gasStation = find_optimal_gas_station(gasRemain, idCurrentCity, idNextCity);
			//find optimal gas station
			// cout << idCurrentCity << " visit optimal gas station " << gasStation.first << " " << gasStation.second << endl;
			
			//find solution to city i
			if (gasStation.first != -1) {
				tempPath.push_back(gasStation);
				ll idGasStation = _graphIdGasStation[gasStation.first][gasStation.second];
				// cout << "cost gas station to city " << idNextCity << " = " << _costGasStationToCity[idGasStation][idNextCity] << endl;
				ll remain = _tankSize - _costGasStationToCity[idGasStation][idNextCity];
				try_to_find_optimal_solution(tempPath, remain, i, sizePath, numCities+1);
				tempPath.pop_back();
			}
			
		}
	}
} //end: try_to_find_optimal_solution

pair<int,int> TSP::find_optimal_gas_station(ll gasRemain, ll idCurrentCity, ll idNextCity) {
	ll maxGasRemainAtNextCity = gasRemain - graph[idCurrentCity][idNextCity];
	pair<int,int> optimalGasStation = make_pair(-1,-1);
	for (vector<pair<int,int> >::iterator it = _nearestGasStation[idCurrentCity].begin(); it != _nearestGasStation[idCurrentCity].end(); ++it) {
		pair<int,int> p = *it;
		ll idGasStation = _graphIdGasStation[p.first][p.second];
		// cout << "try gas station " << p.first << " " << p.second << " cost to next city " << idNextCity << " = " << _costGasStationToCity[idGasStation][idNextCity] << endl;
		if (_costGasStationToCity[idGasStation][idNextCity] > 0 && _costGasStationToCity[idGasStation][idCurrentCity] > 0
				&& gasRemain >= _costGasStationToCity[idGasStation][idCurrentCity]) {

			ll remain = _tankSize - _costGasStationToCity[idGasStation][idNextCity];
			if (maxGasRemainAtNextCity < remain) {
				maxGasRemainAtNextCity = remain;
				optimalGasStation = make_pair(p.first,p.second);
			} else {

			}

		}
	}
	return optimalGasStation;
}//end: find_optimal_gas_station


void TSP::make_shorter(){
	// Modify circuit & pathLength
	twoOpt(graph, circuit, pathLength, n);
}




//================================ PRlong long FUNCTIONS ================================//

void TSP::printResult(){
	ofstream outputStream;
	outputStream.open(outFname.c_str(), ios::out);
	// outputStream << pathLength << endl;
	for (vector<pair<int,int> >::iterator it = _finalPath.begin(); it != _finalPath.end(); ++it) {
		outputStream << (*it).first << " " << (*it).second << endl;
	// for (vector<long long>::iterator it = circuit.begin(); it != circuit.end(); ++it) {
	// 	outputStream << _listGasStation[(*it)].x << " " << _listGasStation[(*it)].y << endl;
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

	cout << "\nLength: " << pathLength << endl << endl;
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
