// Implementation of the KMeans Algorithm
// reference: http://mnemstudio.org/clustering-k-means-example-1.htm

#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>

using namespace std;

class Point
{
private:
	int id_point, id_cluster;
	vector<double> values;
	int dimension;
	string name;

public:
	static double l2_distance(Point& a, Point& b) {
		int dimension = a.dimension;
		double sum = 0.0;

		for(int i = 0; i < dimension; i++)
		{
			sum += pow(a.getValue(i) - b.getValue(i), 2.0);
		}

		return sum;
	}

	static double gaussian_kernel_distance(Point& a, Point& b) {
		double squared_sum = l2_distance(a, b);
		double sigma = 1.0;

		return exp(- squared_sum / (2 * sigma));
	}

	Point(int id_point, vector<double>& values, string name = "")
	{
		this->id_point = id_point;
		dimension = values.size();

		for(int i = 0; i < dimension; i++)
			this->values.push_back(values[i]);

		this->name = name;
		id_cluster = -1;
	}

	int getID()
	{
		return id_point;
	}

	void setCluster(int id_cluster)
	{
		this->id_cluster = id_cluster;
	}

	int getCluster()
	{
		return id_cluster;
	}

	double getValue(int index)
	{
		return values[index];
	}

	int getTotalValues()
	{
		return dimension;
	}

	void addValue(double value)
	{
		values.push_back(value);
	}

	string getName()
	{
		return name;
	}
};

class Cluster
{
private:
	int id_cluster;
	vector<Point> points;

public:
	static double distance_btw_cluster_and_point(Cluster& cluster, Point& point) {
		vector<Point>& cluster_points = cluster.getPoints();
		int cluster_size = cluster_points.size();

		double result = 0.0;

		for (Point xk : cluster_points) {
			result -= 2 * Point::gaussian_kernel_distance(point, xk) / cluster_size;
		}

		for (Point xk : cluster_points) {
			for (Point xl : cluster_points) {
				result += Point::gaussian_kernel_distance(xl, xk) / pow(cluster_size, 2.0);
			}
		}

		// Point cluster_center(-1, cluster.getCentralValues());
		return result;
	}
	
	Cluster(int id_cluster, Point point)
	{
		this->id_cluster = id_cluster;

		int dimension = point.getTotalValues();

		points.push_back(point);
	}

	void addPoint(Point point)
	{
		points.push_back(point);
	}

	bool removePoint(int id_point)
	{
		int total_points = points.size();

		for(int i = 0; i < total_points; i++)
		{
			if(points[i].getID() == id_point)
			{
				points.erase(points.begin() + i);
				return true;
			}
		}
		return false;
	}

	Point getPoint(int index)
	{
		return points[index];
	}

	vector<Point>& getPoints()
	{
		return points;
	}

	int getTotalPoints()
	{
		return points.size();
	}

	int getID()
	{
		return id_cluster;
	}
};

std::ostream& operator<<(std::ostream& os, const std::vector<int> &input)
	{
		for (auto const& i: input) {
			os << i << " ";
		}
		return os;
	}

class KMeans
{
private:
	int K; // number of clusters
	int dimension, total_points, max_iterations;
	vector<Cluster> clusters;

	// return ID of nearest center (uses euclidean distance)
	int getIDNearestCenter(Point point)
	{
		double min_dist;
		int id_cluster_center = 0;

		double dist = Cluster::distance_btw_cluster_and_point(clusters[0], point);

		min_dist = dist;

		for(int i = 1; i < K; i++)
		{
			dist = Cluster::distance_btw_cluster_and_point(clusters[i], point);

			if(dist < min_dist)
			{
				min_dist = dist;
				id_cluster_center = i;
			}
		}

		return id_cluster_center;
	}

public:
	KMeans(int K, int total_points, int dimension, int max_iterations)
	{
		this->K = K;
		this->total_points = total_points;
		this->dimension = dimension;
		this->max_iterations = max_iterations;
	}

	void run(vector<Point> & points)
	{
		cout << "start\n";

		if(K > total_points)
			return;

		vector<int> prohibited_indexes;

		// choose K distinct values for the centers of the clusters
		for(int i = 0; i < K; i++)
		{
			while(true)
			{
				int index_point = rand() % total_points;

				if(find(prohibited_indexes.begin(), prohibited_indexes.end(),
						index_point) == prohibited_indexes.end())
				{
					prohibited_indexes.push_back(index_point);
					points[index_point].setCluster(i);
					Cluster cluster(i, points[index_point]);
					clusters.push_back(cluster);
					break;
				}
			}
		}

		int iter = 1;

		vector<int> id_new_clusters;
		vector<int> id_old_clusters;

		for (int i = 0; i < total_points; i++) {
			id_new_clusters.push_back(-1);
			id_old_clusters.push_back(-1);
		} 

		while(true)
		{
			bool done = true;

			// associates each point to the nearest center
			for(int i = 0; i < total_points; i++)
			{
				id_old_clusters[i] = points[i].getCluster();
				id_new_clusters[i] = getIDNearestCenter(points[i]);
			}

			for (int i = 0; i < total_points; i++) {
				int id_old_cluster = id_old_clusters[i];
				int id_new_center = id_new_clusters[i];

				if(id_old_cluster != id_new_center)
				{
					if(id_old_cluster != -1) {
						clusters[id_old_cluster].removePoint(points[i].getID());
					}
						
					points[i].setCluster(id_new_center);
					clusters[id_new_center].addPoint(points[i]);
					done = false;
				}
			}
			

			if(done == true || iter >= max_iterations)
			{
				cout << "Break in iteration " << iter << "\n\n";
				break;
			}

			iter++;
		}

		// shows elements of clusters
		for(int i = 0; i < K; i++)
		{
			int total_points_cluster =  clusters[i].getTotalPoints();

			cout << "Cluster " << clusters[i].getID() + 1 << endl;
			for(int j = 0; j < total_points_cluster; j++)
			{
				cout << "Point " << clusters[i].getPoint(j).getID() + 1 << ": ";
				for(int p = 0; p < dimension; p++)
					cout << clusters[i].getPoint(j).getValue(p) << " ";

				string point_name = clusters[i].getPoint(j).getName();

				if(point_name != "")
					cout << "- " << point_name;

				cout << endl;
			}

			cout << "\n\n";
		}
	}
};

int main(int argc, char *argv[])
{
	srand (time(NULL));

	int total_points, dimension, K, max_iterations, has_name;

	cin >> total_points >> dimension >> K >> max_iterations >> has_name;

	vector<Point> points;
	string point_name;

	for(int i = 0; i < total_points; i++)
	{
		vector<double> values;

		for(int j = 0; j < dimension; j++)
		{
			double value;
			cin >> value;
			values.push_back(value);
		}

		if(has_name)
		{
			cin >> point_name;
			Point p(i, values, point_name);
			points.push_back(p);
		}
		else
		{
			Point p(i, values);
			points.push_back(p);
		}
	}

	KMeans kmeans(K, total_points, dimension, max_iterations);
	kmeans.run(points);

	return 0;
}
