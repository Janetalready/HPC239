// Implementation of the KMeans Algorithm
// reference: http://mnemstudio.org/clustering-k-means-example-1.htm

#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <omp.h>
#include <chrono>

using namespace std;
using namespace std::chrono;
int NUM_THREAD = 6;
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

		// #pragma omp parallel for reduction(+:sum) num_threads(NUM_THREAD)
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

		// #pragma omp parallel for num_threads(NUM_THREAD)
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

		// #pragma omp parallel for reduction(-:result) num_threads(NUM_THREAD)
		for(int i = 0; i < cluster_points.size(); i++){
			Point xk = cluster_points[i];
			result -= 2 * Point::gaussian_kernel_distance(point, xk) / cluster_size;
		}
		// for (Point xk : cluster_points) {
		// 	result -= 2 * Point::gaussian_kernel_distance(point, xk) / cluster_size;
		// }

		// #pragma omp parallel for collapse(2) reduction(+:result) num_threads(NUM_THREAD)
		for(int i = 0; i <  cluster_points.size(); i++){
			for(int j = 0; j < cluster_points.size(); j++){
				Point xk = cluster_points[i];
				Point xl = cluster_points[j];
				result += Point::gaussian_kernel_distance(xl, xk) / pow(cluster_size, 2.0);
			}
		}
		// for (Point xk : cluster_points) {
		// 	for (Point xl : cluster_points) {
		// 		result += Point::gaussian_kernel_distance(xl, xk) / pow(cluster_size, 2.0);
		// 	}
		// }

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

	void removePoint(int id_point)
	{
		int total_points = points.size();

		// #pragma omp parallel for num_threads(NUM_THREAD)
		for(int i = 0; i < total_points; i++)
		{	
			// #pragma omp critical
			if(points[i].getID() == id_point)
			{
				points.erase(points.begin() + i);
			}
		}
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

		// #pragma omp parallel for private(dist) num_threads(NUM_THREAD)
		for(int i = 1; i < K; i++)
		{
			dist = Cluster::distance_btw_cluster_and_point(clusters[i], point);
			// #pragma omp critical
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
		srand (12345);
		if(K > total_points)
			return;

		vector<int> prohibited_indexes;

		// choose K distinct values for the centers of the clusters
		// #pragma omp parallel for num_threads(NUM_THREAD)
		for(int i = 0; i < K; i++)
		{
			bool flag = true;
			while(flag)
			{
				int index_point = rand() % total_points;

				// #pragma omp critical
				if(find(prohibited_indexes.begin(), prohibited_indexes.end(),
						index_point) == prohibited_indexes.end())
				{
					prohibited_indexes.push_back(index_point);
					points[index_point].setCluster(i);
					Cluster cluster(i, points[index_point]);
					clusters.push_back(cluster);
					flag = false;
				}
			}
		}

		int iter = 1;

		vector<int> id_new_clusters;
		vector<int> id_old_clusters;

		// #pragma omp parallel for num_threads(NUM_THREAD)
		for (int i = 0; i < total_points; i++) {
			id_new_clusters.push_back(-1);
			id_old_clusters.push_back(-1);
		} 

		while(true)
		{
			bool done = true;

			// associates each point to the nearest center
			#pragma omp parallel default(shared) num_threads(NUM_THREAD)
			{
				#pragma omp for schedule(static)
				
					for(int i = 0; i < total_points; i++)
					{
						id_old_clusters[i] = points[i].getCluster();
						id_new_clusters[i] = getIDNearestCenter(points[i]);
					}
				
					// for(int i = 0; i < total_points; i++)
					// {
					// 	id_old_clusters[i] = points[i].getCluster();
					// 	id_new_clusters[i] = getIDNearestCenter(points[i]);
					// }
			}
			// #pragma omp parallel num_threads(NUM_THREAD) default(shared)
			// {
			// 	#pragma omp for schedule(static)
				
			int count_move = 0;
			for (int i = 0; i < total_points; i++) {
				int id_old_cluster = id_old_clusters[i];
				int id_new_center = id_new_clusters[i];

				if(id_old_cluster != id_new_center)
				{
					count_move++;
					if(id_old_cluster != -1) {
						clusters[id_old_cluster].removePoint(points[i].getID());
					}
						
					points[i].setCluster(id_new_center);
					clusters[id_new_center].addPoint(points[i]);
					done = false;
				}
			}
			
			float ratio = float(count_move)/total_points;
			if(done == true || iter >= max_iterations || ratio < 0.01)
			{
				cout << "Break in iteration " << iter << "\n\n";
				break;
			}

			iter++;
		}
	}

	void showClusterElements() {
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

	void showClusterLabelResults(vector<int> labels) {
		//show labels that are put into the same cluster
		for(int i=0; i<clusters.size(); i++) {
			cout << "Cluster " << i+1 << ": ";
			for(Point point : clusters[i].getPoints()) 
				cout << labels[point.getID()] << " ";
			cout << "\n";
		}
	}
};

int main(int argc, char *argv[])
{
	srand (12345);
	int max_iterations = 10, K = 5; //To change

	vector<Point> points;
	ifstream file;
	file.open("doc2vec_reviews.txt"); 
	int id = 0;
	if(file.is_open()) {
		string line;
		while(getline(file, line)) {
			vector<double> values;
    		stringstream s_stream(line); 
    		string value;
		    while(getline(s_stream, value, ',')) {
		       values.push_back(stod(value)); 
		    }
			Point p(id++, values);
			points.push_back(p);
		}
		file.close();
	}

	int total_points = points.size(); //10000
	int dimension = points[0].getTotalValues(); //300


	vector<Point> test_points = std::vector<Point>(points.begin(), points.begin() + 500);
	int test_total_points = test_points.size();
	KMeans kmeans(K, test_total_points, dimension, max_iterations);

	auto start = high_resolution_clock::now(); 
	kmeans.run(points);
	auto stop = high_resolution_clock::now(); 
	auto duration = duration_cast<microseconds>(stop - start); 
    cout << "Time taken by function: " << duration.count() << " microseconds" << endl; 

	vector<int> labels;
	file.open("doc2vec_labels.txt"); 
	if(file.is_open()) {
		string line;
		while(getline(file, line)) {
			int label = stoi(line);
			labels.push_back(label);
		}
		file.close();
	}

	kmeans.showClusterLabelResults(labels);

	return 0;
}
