#include <bits/stdc++.h>
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
			double d = dist(nodePosition[i],nodePosition[j]);
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
		vi currentTour = mutateTour(tour);
		for(int i = 0; i < V; i++){
  			for(int j = i+1; j < V; j++){
  				vi new_tour;
  				for(int k = 0; k < i; k++){
  					new_tour.push_back(currentTour[k]);
  				}
  				for(int k = i; k <= j; k++){
  					new_tour.push_back(currentTour[j-(k-i)]);
  				}
  				for(int k = j+1; k < V; k++){
  					new_tour.push_back(currentTour[k]);
  				}

  				double newDist = computeTourDistance(new_tour, nodeList);
			//	printf("%d\n", newDist);
  			//	printf("%lf\n", newDist);
  				if(newDist < bestDistance){
  					tour = new_tour;
  					bestDistance = newDist;
           			print_vector(tour);
  				}
  			}
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
		double best = 999999999;
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

int main()
{
	srand (time(NULL)); // Randomize seed
	//randomGraph(10,2,1);
	print_vector(twoOptAlgorithm(randomGraphNodeList(10, 100, 1),1));
	print_vector(nearestNeighbourHeuristic(randomGraph(10, 2, 1)));
}

//-----------------------------------------
//	Helper Methods
//-----------------------------------------

double dist(point node1, point node2)
{
	return sqrt((node1.first - node2.first)*(node1.first - node2.first) + (node1.second - node2.second)*(node1.second - node2.second));
}

double greatCircleDistance(point node1, point node2){
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

double computeTourDistance(vi tour, vector<point> nodeList){
	double result = 0.0;
	int V = (int)nodeList.size();
	for(int i = 0; i < V-1; i++){
		result += greatCircleDistance(nodeList[i], nodeList[i+1]);
	}
	result += greatCircleDistance(nodeList[V-1], nodeList[0]);
	return result;
}

vi mutateTour(vi tour){
	int V = (int)tour.size();
	vi finalTour = tour;
	for(int tt = 0; tt < 2; tt++){
		vi newTour;
		int x = rand() % V;
		int y = rand() % V;
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