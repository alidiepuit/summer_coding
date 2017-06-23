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
	bool **visitedGraph;
    visitedGraph = new bool*[_row+1];
	for (long long i = 0; i < _row+1; i++) {
		visitedGraph[i] = new bool[_col+1];
		for (long long j = 0; j < _col+1; j++) {
			visitedGraph[i][j] = 0;
		}
	}
    long long idGasStation = _graphIdGasStation[x][y];

    long long maxV = 0;
    // cout << "fuck " << c.x << " " << c.y << endl;
    visitedGraph[c.x][c.y] = 1;
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
            if (isValidPosition(newx, newy) && visitedGraph[newx][newy]==0 && maxV+1 <= _tankSize) {
                visitedGraph[newx][newy] = 1;
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
                	// cout << idGasStation << " " << newx << " " << newy << " " << maxV+1 << endl;
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
	memset(parent, -1, _numGasStation * sizeof(ll));

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

	_connectedGraph = new vector<long long> [_numGasStation];
	_numGasStationConnected = 0;
	for(ll i = 0; i < _numGasStation; i++) {
		for(ll j = 0; j < _numGasStation; j++) {
			if (parent[j] == i) {
				_connectedGraph[i].push_back(j);
			}
		}
	}

	// for(int i = 0; i < _row; i++) {
	// 	if (_graphIdGasStation[i][130] >= 0) {
	// 		cout << __LINE__ << " " << i << endl;
	// 	}
	// }
	// ll idGasStation = _graphIdGasStation[235][130];
	// ll idCity1 = _graphId[89][28];
	// ll idCity2 = _graphId[81][36];
	// cout << __LINE__ << " " << idGasStation << endl;

	_markCity = new bool[n];
	memset(_markCity, 0, n * sizeof(bool));
	_markGasStation = new bool[_numGasStation];
	memset(_markGasStation, 0, _numGasStation * sizeof(bool));

	_finalPath.push_back(_positionStart);
	_markCity[_graphId[_positionStart.first][_positionStart.second]] = 1;

	// cout << __LINE__ << " fuck" << endl;
	recursive_connected_graph(idFirstGasStation);

	/*
	bool hasLastHeart = false;
	vector<pair<int,int> > tmp;
	for(vector<pair<int,int> >::iterator it = _finalPath.begin(); it != _finalPath.end(); it++) {
		pair<int,int> coor = *it;
		if (_graphIdGasStation[coor.first][coor.second] > 0) {
			ll idGasStation = _graphIdGasStation[coor.first][coor.second];
			for(vector<pair<int,int> >::iterator city = _nearestCity[idGasStation].begin(); 
												city != _nearestCity[idGasStation].end(); city++) {
				ll idCity = _graphId[(*city).first][(*city).second];
				if (!_markCity[idCity]) {
					tmp.push_back(*city);
					hasLastHeart = true;
					break;
				}
			}
			if (hasLastHeart) break;
			tmp.push_back(make_pair(coor.first,coor.second));
		}
	}

	if (hasLastHeart) {
		for(vector<pair<int,int> >::iterator it = tmp.begin(); it != tmp.end(); it++)
			_finalPath.push_back(*it);
	}
	*/
};

vector<pair<int,int> > TSP::greedy_single_gas_station(vector<pair<int,int> > pathCities, 
	pair<int,int> gasStation) {
	ll i = 0;
	ll length = pathCities.size();
	vector<pair<int,int> > res;
	ll gasRemain = _tankSize;
	ll indexGasStation = _graphIdGasStation[gasStation.first][gasStation.second];
	// res.push_back(gasStation);
	vector<ll> cityInPath;
	while (i < length) {
		pair<int,int> coor = pathCities[i];
		// cout << __LINE__ << " " << i  << " " << coor.first << " " << coor.second << endl;
		ll indexCurrentCity = _graphId[coor.first][coor.second];
		City currentCity = cities[indexCurrentCity];
		// cout << "mark city " << indexCurrentCity << " " << _markCity[indexCurrentCity] << endl;
		if (!_markCity[indexCurrentCity] && _costGasStationToCity[indexGasStation][indexCurrentCity] <= _tankSize/2) {
			// cout << __LINE__ << endl;
			if (res.size() == 0) {
				// cout << __LINE__ << endl;

				gasRemain -= _costGasStationToCity[indexGasStation][indexCurrentCity];
				res.push_back(make_pair(currentCity.x,currentCity.y));
				cityInPath.push_back(indexCurrentCity);

				if (indexCurrentCity == 215) {
					cout << __LINE__ << " " << indexCurrentCity << " " << gasRemain << endl;
				}
				// cout << __LINE__ << endl;
			} else {
				// cout << i << " fuck " << endl;
				City coor = cities[*(cityInPath.end()-1)];
				ll indexPrevCity = _graphId[coor.x][coor.y];
				ll costPrevCityToCurrentCity = graph[indexCurrentCity][indexPrevCity];
				// cout << costPrevCityToCurrentCity << " fuck " << endl;
				ll costGasStationToCurrentCity = _costGasStationToCity[indexGasStation][indexCurrentCity];

				if (costPrevCityToCurrentCity <= 0 || costGasStationToCurrentCity <= 0) {
					cout << __LINE__ << " " << indexPrevCity << " " << graph[indexCurrentCity][indexPrevCity] << " " << costGasStationToCurrentCity << endl;
				}

				if (costPrevCityToCurrentCity > 0 && gasRemain - costPrevCityToCurrentCity - costGasStationToCurrentCity >= 0) {
					gasRemain -= costPrevCityToCurrentCity;
					res.push_back(make_pair(currentCity.x,currentCity.y));

					cityInPath.push_back(indexCurrentCity);
				} else {
					res.push_back(gasStation);
					res.push_back(make_pair(currentCity.x,currentCity.y));

					cityInPath.push_back(indexCurrentCity);
					gasRemain = _tankSize - costGasStationToCurrentCity;
				}
			}
	// cout << __LINE__ << endl;
			_markCity[indexCurrentCity] = 1;
		}

	// cout << __LINE__ << endl;
		//visit next city
		i++;
	}
	// cout << __LINE__ << endl;
	//insert gas station to path
	// res.push_back(gasStation);
	return res;
}

void TSP::recursive_connected_graph(ll idGasStation) {
	City coorGasStation = _listGasStation[idGasStation];
	_finalPath.push_back(make_pair(coorGasStation.x,coorGasStation.y));

	

	// vector<pair<int,int> > path = greedy_single_gas_station(_nearestCity[idGasStation], 
												// make_pair(coorGasStation.x,coorGasStation.y));



	// if (coorGasStation.x == 234 && coorGasStation.y == 130) {
	// 	ll idCity1 = _graphId[249][139];
	// 	ll idCity2 = _graphId[253][128];
	// 	cout << __LINE__ << " " << idCity1 << " " << idCity2 << " " << graph[idCity1][idCity2] << endl;
	// 	for(vector<pair<int,int> >::iterator i = path.begin(); i != path.end(); i++) {
	// 		ll idCity = _graphId[(*i).first][(*i).second];
	// 		cout << (*i).first << "," << (*i).second << " " << _costGasStationToCity[idGasStation][idCity] << endl;	
	// 	}
	// }

	for(vector<pair<int,int> >::iterator i = _nearestCity[idGasStation].begin(); i != _nearestCity[idGasStation].end(); i++) {
		ll idCity = _graphId[(*i).first][(*i).second];
		if (!_markCity[idCity] && _costGasStationToCity[idGasStation][idCity] <= _tankSize/2) {
			_finalPath.push_back(make_pair((*i).first,(*i).second));
			_finalPath.push_back(make_pair(coorGasStation.x,coorGasStation.y));
			_markCity[idCity] = 1;
		}
	}

	// cout << "visited path size: " << path.size() << endl;
	// if (path.size() > 0) {
	// 	for(vector<pair<int,int> >::iterator i = path.begin(); i != path.end(); i++)
	// 		_finalPath.push_back(*i);
	// 	_finalPath.push_back(make_pair(coorGasStation.x,coorGasStation.y));
	// }

	// cout << "has child gas station: " << _connectedGraph[idGasStation].size() << endl;
	for(vector<ll>::iterator i = _connectedGraph[idGasStation].begin();
				i != _connectedGraph[idGasStation].end(); i++) {
		recursive_connected_graph(*i);
		_finalPath.push_back(make_pair(coorGasStation.x,coorGasStation.y));
	}
}

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
	euler(pos, circuit);

	// printEuler();

	// make it hamiltonian
	// pass actual vars
	make_hamilton(circuit, pathLength);
}



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
		if (it != _finalPath.begin()) {
			pair<int,int> pPrev = *(it-1);
			pair<int,int> p = *it;
			if (p.first == pPrev.first && p.second == pPrev.second) {
				continue;
			}
		}
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

void TSP::initGraph() {
    for (ll i = 0; i < n; i++) {
        City city = cities[i];
        // cout << "begin " << city.x << " " << city.y << endl;
        floatMatrix(city.x, city.y, _tankSize);
    }
}

void TSP::floatMatrix(int x, int y, long long tankSize) {
    queue<pair<int, City> > vec;
    struct City c = {x, y};
    vec.push(make_pair(0, c));
    bool **visitedGraph;
    visitedGraph = new bool*[_row+1];
	for (long long i = 0; i < _row+1; i++) {
		visitedGraph[i] = new bool[_col+1];
		for (long long j = 0; j < _col+1; j++) {
			visitedGraph[i][j] = 0;
		}
	}
    long long idCity = _graphId[x][y];

    long long maxV = 0;
    if (_graphId[x][y] == 196) {
    	// cout << "fuck " << c.x << " " << c.y << endl;
    } 
    visitedGraph[c.x][c.y] = 1;
    while (vec.size() > 0) {
        pair<ll, City> p = vec.front();
        vec.pop();
        City city = p.second;
        maxV = p.first;
        
        if (maxV >= _tankSize) {
        	continue;
        }

        for(int i = 0; i < 4; i++) {
            int newx = city.x + DIRX[i];
            int newy = city.y + DIRY[i];
            struct City c = {newx, newy};
            if (_graphId[x][y] == 196) {
            	// cout << "check isValidPosition " << newx << " " << newy << " " << isValidPosition(newx, newy)  << " " << visitedGraph[newx][newy] << endl;
            }
            if (isValidPosition(newx, newy) && visitedGraph[newx][newy]==0 && maxV+1 <= _tankSize) {
                visitedGraph[newx][newy] = 1;
                vec.push(make_pair(maxV+1,c));

                //check whether visit city
                if (_originMap[newx][newy] == 3) {
                    graph[idCity][_graphId[newx][newy]] = maxV+1;
                    graph[_graphId[newx][newy]][idCity] = maxV+1;
                    if (_graphId[x][y] == 215 && _graphId[newx][newy] == 196) {
				    	// cout << __LINE__ << " fuck " << maxV+1 << endl;
				    }  
                }
            }
        }
    }
}
