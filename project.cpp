#include <cstdio>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include <queue>
#include <stack>
#include <assert.h>
using namespace std;

typedef pair<int, double> id;
typedef pair<double, double> point;
typedef vector<int> vi;
typedef vector<id> vid;
typedef vector<vid> graph;

#define PI acos(-1.0)
#define EARTH_RAD (6371009) 

//-----------------------------------------
//	Helper Methods Signatures
//-----------------------------------------

double dist(point node1, point node2);
double greatCircleDistance(point node1, point node2);
double computeTourDistance(vi tour, vector<point> nodeList);
vi mutateTour(vi tour);
void print_vector(vi v);

//-----------------------------------------
//	Main Methods
//-----------------------------------------

// Returns a random graph in the form of an adjacency list
graph randomGraph(int n, double D, double p)
{
	vector<point> nodePosition;
	nodePosition.resize(n);

	for (int i = 0; i < n; ++i)
	{
		nodePosition[i].first = (double)rand() / RAND_MAX * D;
		nodePosition[i].second = (double)rand() / RAND_MAX * D;
	}

	graph randomResult; randomResult.resize(n);
	for (int i = 0; i < n; ++i)
	{
		for (int j = i + 1; j < n; ++j)
		{
			double d = greatCircleDistance(nodePosition[i],nodePosition[j]);
			//printf("%d %d %lf %lf %lf %lf %lf\n", i, j, nodePosition[i].first, nodePosition[i].second, nodePosition[j].first, nodePosition[j].second, d);
			randomResult[i].push_back(make_pair(j,d));
			randomResult[j].push_back(make_pair(i,d));
			
		}
	}
	return randomResult;
}

vector<point> randomGraphNodeList(int n, double D, double p)
{
	vector<point> nodePosition;
	nodePosition.resize(n);

	for (int i = 0; i < n; ++i)
	{
		nodePosition[i].first = (double)rand() / RAND_MAX * D;
		nodePosition[i].second = (double)rand() / RAND_MAX * D;
	}

	return nodePosition;
}

graph nodeListToAdjList(vector<point> nodeList)
{
	graph adjList;
	adjList.resize(nodeList.size());
	for (int i = 0; i < nodeList.size(); ++i)
	{
		for (int j = 0; j < nodeList.size(); ++j)
		{
			if (i != j)
				adjList[i].push_back(make_pair(j,greatCircleDistance(nodeList[i],nodeList[j])));
		}
	}
	return adjList;
}

vi twoOptAlgorithm(vector<point> nodeList, int runTime){
	int startTime = time(NULL);
	int V = (int)nodeList.size();
	
	// Start with a tour from 0 to V-1
	vi tour;
	for(int i = 0; i < V; i++){
		tour.push_back(i);
	}

	double bestDistance = computeTourDistance(tour, nodeList);
	while(time(NULL) - startTime < runTime){
		//printf("here\n");
		int beststart = -1;
		int bestend = -1;
		vi currentTour = mutateTour(tour);
		double currentDistance = computeTourDistance(currentTour, nodeList);

		for(int i = 1; i < V; i++){
  			for(int j = i+1; j < V; j++){
  			//	print_vector(currentTour);
  			//	print_vector(new_tour);
  				double newDist = currentDistance - greatCircleDistance(nodeList[currentTour[i-1]], nodeList[currentTour[i]])
  				 - greatCircleDistance(nodeList[currentTour[j]], nodeList[currentTour[j+1==V?0:j]]) 
  				 + greatCircleDistance(nodeList[currentTour[i-1]], nodeList[currentTour[j]]) 
  				 + greatCircleDistance(nodeList[currentTour[i]], nodeList[currentTour[j+1==V?0:j]]) ;
			//	printf("%d\n", newDist);
  			//	printf("%lf\n", newDist);
  				if(newDist < bestDistance){
  					beststart = i;
  					bestend = j;
  					bestDistance = newDist;
           			//print_vector(tour);
  				}
  			}
  		}
  		if(beststart != -1){
  			vi new_tour;
			for(int k = 0; k < beststart; k++){
				new_tour.push_back(currentTour[k]);
			}
			for(int k = beststart; k <= bestend; k++){
				new_tour.push_back(currentTour[bestend-(k-beststart)]);
			}
			for(int k = bestend+1; k < V; k++){
				new_tour.push_back(currentTour[k]);
			}
			tour = new_tour;
  		}
	}
	return tour;
}

vi nearestNeighbourHeuristic(graph adjList){
	int V = (int)adjList.size();
	vi tour;
	vi selected;
	selected.assign(V, 0);
	tour.push_back(0);
	selected[0] = 1;
	for(int i = 1; i < V; i++){
		int prev = tour[i-1];
		//printf("%d\n", prev);
		double best = 9999999999999999;
		int bestIndex = 0;
		//print_vector(selected);
		for(int j = 0; j < adjList[prev].size(); j++){
			int node = adjList[prev][j].first;
			double currentWeight = adjList[prev][j].second;
			//printf("%d %lf\n", node, currentWeight);
			if(selected[node] == 0 && currentWeight < best){
				best = currentWeight;
				bestIndex = node;
			}
		}
		selected[bestIndex] = 1;
		tour.push_back(bestIndex);
	}
	return tour;
}

graph findMST(graph adjList) {
	graph MST;
	priority_queue<pair<double,pair<int,int> > > pq;
	vector<bool> visited;

	MST.resize(adjList.size());
	visited.assign(adjList.size(),false);
	pq.push(make_pair(0,make_pair(0,-1)));

	while (!pq.empty()) {
		double dist = pq.top().first;
		int vertex_now = pq.top().second.first;
		int vertex_prev = pq.top().second.second;
		pq.pop();
		if (visited[vertex_now]) continue;
		visited[vertex_now] = true;
		if (vertex_prev != -1) {
			MST[vertex_now].push_back(make_pair(vertex_prev,dist));
			MST[vertex_prev].push_back(make_pair(vertex_now,dist));
		}
		for (int i = 0; i < adjList[vertex_now].size(); ++i) {
			double dist = adjList[vertex_now][i].second;
			int vertex_next = adjList[vertex_now][i].first;
			pq.push(make_pair(-dist,make_pair(vertex_next,vertex_now)));
		}
	}
	return MST;
}

vi twoApproxAlgorithm(graph adjList){
	vector<bool> visited;
	vector<int> tour;
	stack<int> DFS; //let's do DFS without recursion, shall we

	visited.assign(adjList.size(),false);
	graph MST = findMST(adjList);
	DFS.push(0);

	while (!DFS.empty()) {
		int vertex_now = DFS.top();
		DFS.pop();
		if (visited[vertex_now]) continue;
		visited[vertex_now] = true;
		tour.push_back(vertex_now);
		for (int i = 0; i < MST[vertex_now].size(); ++i) {
			int vertex_next = MST[vertex_now][i].first;
			DFS.push(vertex_next);
		}
	}

	return tour;
}

double memo[18][1 << 18];

double heldKarpTSP(int pos, int mask, vector<point> nodeList){
	//printf("%d %d\n", pos, mask);
	int N = nodeList.size();
	if(mask == (1 << N)-1) {
		//printf("here\n");
		return greatCircleDistance(nodeList[0], nodeList[pos]);
	}
	else if(memo[pos][mask] >= 0.0) return memo[pos][mask];

	double ans = 100000000000.0;
	for(int i = 0; i < N; i++){
		if((mask & (1 << i)) == 0 && i != pos) {
			//printf("%lf\n", greatCircleDistance(nodeList[pos], nodeList[i]) + heldKarpTSP(i, mask + (1 << i), nodeList));
			ans = min(ans, greatCircleDistance(nodeList[pos], nodeList[i]) + heldKarpTSP(i, mask | (1 << i), nodeList));
		}
	}
	return memo[pos][mask] = ans;
}

double optimalTour(vector<point> nodeList){
	int N = nodeList.size();
	for(int i = 0; i < N; i++){
		for(int j = 0; j < (1 << N); j++){
			memo[i][j] = -1.0;
		}
	}
	return heldKarpTSP(0, 1, nodeList);
}

void outputGraphToFile(string fileName, vector<point> nodeList) 
{
	freopen(fileName.c_str(), "w", stdout);	
	for (int i = 0; i < nodeList.size(); ++i)
	{
		printf("%.6lf %.6lf\n",nodeList[i].first,nodeList[i].second);
	}
}

vector<point> inputGraphFromFile(string fileName) 
{
	freopen(fileName.c_str(), "r", stdin);	
	vector<point> graphResult;
	double x,y;
	while (scanf("%lf %lf",&x,&y)!=EOF) {
		graphResult.push_back(make_pair(x,y));
	}
	return graphResult;
}

void generateRandomGraphs(void)
{
	outputGraphToFile("data_random/random1.txt",randomGraphNodeList(10,100,1));
	outputGraphToFile("data_random/random2.txt",randomGraphNodeList(20,100,1));
	outputGraphToFile("data_random/random3.txt",randomGraphNodeList(500,100,1));
	outputGraphToFile("data_random/random4.txt",randomGraphNodeList(1000,1,1));
	outputGraphToFile("data_random/random5.txt",randomGraphNodeList(1000,1000,1));
}

void generateSmallRandomGraphs(void)
{
	outputGraphToFile("data_small_random/random1.txt",randomGraphNodeList(15,10,1));
	outputGraphToFile("data_small_random/random2.txt",randomGraphNodeList(15,10,1));
	outputGraphToFile("data_small_random/random3.txt",randomGraphNodeList(15,10,1));
	outputGraphToFile("data_small_random/random4.txt",randomGraphNodeList(15,10,1));
	outputGraphToFile("data_small_random/random5.txt",randomGraphNodeList(15,10,1));
	outputGraphToFile("data_small_random/random6.txt",randomGraphNodeList(15,10,1));
	outputGraphToFile("data_small_random/random7.txt",randomGraphNodeList(15,10,1));
	outputGraphToFile("data_small_random/random8.txt",randomGraphNodeList(15,10,1));
	outputGraphToFile("data_small_random/random9.txt",randomGraphNodeList(15,10,1));
	outputGraphToFile("data_small_random/random10.txt",randomGraphNodeList(15,10,1));
}

void generateMidRandomGraphs(void)
{
	outputGraphToFile("data_mid_random/random1.txt",randomGraphNodeList(1000,10,1));
	outputGraphToFile("data_mid_random/random2.txt",randomGraphNodeList(1000,10,1));
	outputGraphToFile("data_mid_random/random3.txt",randomGraphNodeList(1000,10,1));
	outputGraphToFile("data_mid_random/random4.txt",randomGraphNodeList(1000,10,1));
	outputGraphToFile("data_mid_random/random5.txt",randomGraphNodeList(1000,10,1));
	outputGraphToFile("data_mid_random/random6.txt",randomGraphNodeList(1000,10,1));
	outputGraphToFile("data_mid_random/random7.txt",randomGraphNodeList(1000,10,1));
	outputGraphToFile("data_mid_random/random8.txt",randomGraphNodeList(1000,10,1));
	outputGraphToFile("data_mid_random/random9.txt",randomGraphNodeList(1000,10,1));
	outputGraphToFile("data_mid_random/random10.txt",randomGraphNodeList(1000,10,1));
}

void generateBigRandomGraphs(void)
{
	outputGraphToFile("data_big_random/random1.txt",randomGraphNodeList(10000,10,1));
	outputGraphToFile("data_big_random/random2.txt",randomGraphNodeList(10000,10,1));
	outputGraphToFile("data_big_random/random3.txt",randomGraphNodeList(10000,10,1));
	outputGraphToFile("data_big_random/random4.txt",randomGraphNodeList(10000,10,1));
	outputGraphToFile("data_big_random/random5.txt",randomGraphNodeList(10000,10,1));
	outputGraphToFile("data_big_random/random6.txt",randomGraphNodeList(10000,10,1));
	outputGraphToFile("data_big_random/random7.txt",randomGraphNodeList(10000,10,1));
	outputGraphToFile("data_big_random/random8.txt",randomGraphNodeList(10000,10,1));
	outputGraphToFile("data_big_random/random9.txt",randomGraphNodeList(10000,10,1));
	outputGraphToFile("data_big_random/random10.txt",randomGraphNodeList(10000,10,1));
}

void test(string inputFile)
{
	vector<point> nodeList = inputGraphFromFile(inputFile);
	graph adjList = nodeListToAdjList(nodeList);
	vector<int> tour = twoApproxAlgorithm(adjList);
	printf("%.3lf\n",computeTourDistance(tour,nodeList));
}

void experiment1(){
	string files[] = {"data_small_random/random1.txt",
					"data_small_random/random2.txt",
					"data_small_random/random3.txt",
					"data_small_random/random4.txt",
					"data_small_random/random5.txt",
					"data_small_random/random6.txt",
					"data_small_random/random7.txt",
					"data_small_random/random8.txt",
					"data_small_random/random9.txt",
					"data_small_random/random10.txt"};

	string testname[] = {"Random Graph 1",
						 "Random Graph 2",
						 "Random Graph 3",
						 "Random Graph 4",
						 "Random Graph 5",
						 "Random Graph 6",
						 "Random Graph 7",
						 "Random Graph 8",
						 "Random Graph 9",
						 "Random Graph 10"};
	srand (time(NULL)); 
	for (int i = 0; i < 10; ++i)
	{
		printf("ON DATASET %d: %s\n",i,testname[i].c_str());
		fflush(stdout);
		vector<point> nodeList = inputGraphFromFile(files[i]);
		graph adjList = nodeListToAdjList(nodeList);

		double startTime;

		printf("2-approximation algorithm :\n");
		startTime = time(NULL);
		vector<int> twoApproxTour = twoApproxAlgorithm(adjList);
		printf("Running time   : %.3lf s\n",time(NULL) - startTime);
		printf("Distance of tour produced : %.3lf\n",computeTourDistance(twoApproxTour,nodeList));
		puts("");

		printf("Held-Karp algorithm :\n");
		startTime = time(NULL);
		printf("Running time   : %.3lf ms\n",time(NULL) - startTime);
		printf("Distance of tour produced : %.3lf\n",optimalTour(nodeList));
		puts("");


		printf("\n");
	}
}

void experiment2(){
	string files[] = {"data_mid_random/random1.txt",
					"data_mid_random/random2.txt",
					"data_mid_random/random3.txt"};

	string testname[] = {"Random Graph 1",
						 "Random Graph 2",
						 "Random Graph 3"};
	srand (time(NULL)); 
	for (int i = 0; i < 3; ++i)
	{
		for(int j = 10; j <= 120; j += 10){

			printf("ON DATASET %d: %s\n",i,testname[i].c_str());
			fflush(stdout);
			vector<point> nodeList = inputGraphFromFile(files[i]);
			graph adjList = nodeListToAdjList(nodeList);

			double startTime;

			printf("2-OPT algorithm %d seconds:\n", j);
			startTime = time(NULL);
			vector<int> twoOptTour = twoOptAlgorithm(nodeList,j);
			printf("Running time   : %.3lf ms\n",time(NULL) - startTime);
			printf("Distance of tour produced : %.3lf\n",computeTourDistance(twoOptTour,nodeList));
			puts("");
		}

		printf("\n");
	}
}

void experiment3()
{
	string files[] = {"data_random/random1.txt",
					  "data_random/random2.txt",
					  "data_random/random3.txt",
					  "data_random/random4.txt",
					  "data_random/random5.txt",
					  "data_reduced/911.txt",
					  "data_reduced/nypd.txt"};

	string testname[] = {"Random Graph 1",
						 "Random Graph 2",
						 "Random Graph 3",
						 "Random Graph 4",
						 "Random Graph 5",
						 "Real Graph 1",
						 "Real Graph 2" };

	//freopen("output.txt", "w", stdout);
	srand (time(NULL)); // Randomize seed
	//test(files[5]);
	//return 0;

	for (int i = 0; i < 7; ++i)
	{
		printf("ON DATASET %d: %s\n",i,testname[i].c_str());
		fflush(stdout);
		vector<point> nodeList = inputGraphFromFile(files[i]);
		graph adjList = nodeListToAdjList(nodeList);

		double startTime;

		printf("2-approximation algorithm :\n");
		startTime = time(NULL);
		vector<int> twoApproxTour = twoApproxAlgorithm(adjList);
		printf("Running time   : %.3lf s\n",time(NULL) - startTime);
		printf("Distance of tour produced : %.3lf\n",computeTourDistance(twoApproxTour,nodeList));
		puts("");

		printf("2-OPT algorithm :\n");
		startTime = time(NULL);
		vector<int> twoOptTour = twoOptAlgorithm(nodeList,30);
		printf("Running time   : %.3lf ms\n",time(NULL) - startTime);
		printf("Distance of tour produced : %.3lf\n",computeTourDistance(twoOptTour,nodeList));
		puts("");

		printf("Nearest neighbour heuristic:\n");
		startTime = time(NULL);
		vector<int> nearestNeighbourTour = nearestNeighbourHeuristic(adjList);
		printf("Running time   : %.3lf s\n",time(NULL) - startTime);
		printf("Distance of tour produced : %.3lf\n",computeTourDistance(nearestNeighbourTour,nodeList));
		puts("");


		printf("\n");
	}

}

int main(){
	experiment2();
}

//-----------------------------------------
//	Helper Methods
//-----------------------------------------

double dist(point node1, point node2)
{
	return sqrt((node1.first - node2.first)*(node1.first - node2.first) + (node1.second - node2.second)*(node1.second - node2.second));
}

double greatCircleDistance(point node1, point node2){
	if(node1.first == node2.first && node1.second == node2.second) return 0.0;
	double pLat = node1.first;
	double pLong = node1.second;
	double qLat = node2.first;
	double qLong = node2.second;
	double radius = EARTH_RAD;
	pLat *= PI / 180; pLong *= PI / 180;
  	qLat *= PI / 180; qLong *= PI / 180;
  	return radius * acos(cos(pLat)*cos(pLong)*cos(qLat)*cos(qLong) +
                       cos(pLat)*sin(pLong)*cos(qLat)*sin(qLong) +
                       sin(pLat)*sin(qLat));
}

bool isValidTour(vector<int> tour)
{
	sort(tour.begin(),tour.end());
	for (int i = 0; i < tour.size(); ++i)
		if (tour[i] != i) return false;
	return true;
}

double computeTourDistance(vi tour, vector<point> nodeList){
	double result = 0.0;
	assert(tour.size() == nodeList.size());
	assert(isValidTour(tour));
	int V = (int)nodeList.size();
	for(int i = 0; i < V-1; i++){
		result += greatCircleDistance(nodeList[tour[i]], nodeList[tour[i+1]]);
		//printf("(%.3lf,%.3lf)(%.3lf,%.3lf)%.3lf\n",nodeList[tour[i]].first,nodeList[tour[i]].second,nodeList[tour[i+1]].first,nodeList[tour[i+1]].second,result);
	}
	result += greatCircleDistance(nodeList[V-1], nodeList[0]);
	//printf("%lf\n", result);
	return result;
}

vi mutateTour(vi tour){
	int V = (int)tour.size();
	vi finalTour = tour;
	for(int tt = 0; tt < 2; tt++){
		vi newTour;
		int x = rand() % V;
		int y = rand() % V;
		if(x > y) swap(x,y);
		for(int k = 0; k < x; k++){
			newTour.push_back(finalTour[k]);
		}
		for(int k = x; k <= y; k++){
			newTour.push_back(finalTour[y-(k-x)]);
		}
		for(int k = y+1; k < V; k++){
			newTour.push_back(finalTour[k]);
		}
		finalTour = newTour;
	}
	// print_vector(finalTour);
	return finalTour;
}

void print_vector(vi v){
  printf("[");
  for(int i = 0; i < (int)v.size(); i++){
    printf("%d ", v[i]);
  }
  printf("]");
  printf("\n");
}