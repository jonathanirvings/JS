#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
using namespace std;

typedef pair<int, double> id;
typedef vector<id> vid;
typedef vector<vid> graph;

double distance(pair<double,double> node1, pair<double,double> node2)
{
	return sqrt((node1.first - node2.first) * (node1.second - node2.second));
}

graph randomGraph(int n, double D, double p)
{
	vector<pair<double, double> > nodePosition;
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
			double dist = distance(nodePosition[i],nodePosition[j]);
			if ((double)rand() / RAND_MAX < p)
			{
				randomResult[i].push_back(make_pair(j,dist));
				randomResult[j].push_back(make_pair(i,dist));
			}
		}
	}
	return randomResult;
}

int main()
{
	randomGraph(10,2,1);
}