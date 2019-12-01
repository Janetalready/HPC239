# include <cstdlib>
# include <ctime>
# include <iomanip>
# include <iostream>
# include <mpi.h>
#include <vector>
#include <math.h>
#include <algorithm>
#include <fstream>
#include <sstream>

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
  vector<int> points_ids;

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

  static double distance_btw_cluster_and_point_only_index(vector<int>& cluster, Point& point, vector<Point>& points) {
		vector<Point> cluster_points;
    for(int i = 0; i < cluster.size(); i++){
      int idx = cluster[i];
      cluster_points.push_back(points[idx]);
    }
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
    points_ids.push_back(point.getID());
	}

	bool removePoint(int id_point)
	{
		int total_points = points.size();

		for(int i = 0; i < total_points; i++)
		{
			if(points[i].getID() == id_point)
			{
				points.erase(points.begin() + i);
        points_ids.erase(points_ids.begin() + i);
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

  vector<int>& getPoints_ids()
	{
		return points_ids;
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

public:
  vector<Cluster> clusters;
	KMeans(int K, int total_points, int dimension, int max_iterations)
	{
		this->K = K;
		this->total_points = total_points;
		this->dimension = dimension;
		this->max_iterations = max_iterations;
	}

  void edit_point_from_cluster(int id_old_cluster, int id_new_center, Point point)
  {
    if(id_old_cluster != -1) {
      clusters[id_old_cluster].removePoint(point.getID());
    }
      
    point.setCluster(id_new_center);
    clusters[id_new_center].addPoint(point);
  }

  void add_cluster(Cluster cluster)
  {
    clusters.push_back(cluster);
  }

  vector<Cluster> get_cluster()
  {
    return clusters;
  }

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

	void showClusterResults(vector<int> labels) {
		for(int i=0; i<clusters.size(); i++) {
			cout << "Cluster " << i+1 << ": ";
			for(Point point : clusters[i].getPoints()) 
				cout << labels[point.getID()] << " ";
			cout << "\n";
		}
	}
};

int main(int argc, char** argv) {

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

	vector<Point> test_points = std::vector<Point>(points.begin(), points.begin() + 300);
	int test_total_points = test_points.size();
  KMeans kmeans(K, test_total_points, dimension, max_iterations);

    // number of sites per processor.
  int sites_per_proc = 60;
  srand (time(NULL));

  // Initial MPI and find process rank and number of processes.
  MPI_Init(NULL, NULL);
  int rank, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  // The cluster assignments for each site.
  vector<int> id_new_clusters_per_pro;
  vector<vector<int>> total_cluster;
  vector<int> cluster_size;


  for (int i = 0; i < sites_per_proc; i++) {
			id_new_clusters_per_pro.push_back(-1);
	} 

  //
  // Data structures maintained only in root process.
  //
  // Result of program: a cluster label for each site.
  vector<int> id_new_clusters;
	vector<int> id_old_clusters;
  vector<int> cluster_concate;
  int clustered_points_size;
  if (rank == 0) {
    if(K > test_total_points)
			return -1;

		vector<int> prohibited_indexes;

    for (int i = 0; i < test_total_points; i++) {
			id_new_clusters.push_back(-1);
			id_old_clusters.push_back(-1);
		} 
		// choose K distinct values for the centers of the clusters
		for(int i = 0; i < K; i++)
		{
			while(true)
			{
				int index_point = rand() % test_total_points;

				if(find(prohibited_indexes.begin(), prohibited_indexes.end(),
						index_point) == prohibited_indexes.end())
				{
					prohibited_indexes.push_back(index_point);
					vector<int> cluster;
          cluster.push_back(index_point);
          id_old_clusters[index_point] = i;
          total_cluster.push_back(cluster);
          cluster_size.push_back(1);
					break;
				}
			}
		}

    clustered_points_size = accumulate(cluster_size.begin(), cluster_size.end(), 0);
    cluster_concate.reserve(clustered_points_size);
    for(int i = 0; i < K; i++){
      cluster_concate.insert(cluster_concate.end(), total_cluster[i].begin(), total_cluster[i].end());
    }

  }

  int iter = 1;
  
  bool done = true;
  while (true) { // While they've moved...


    // Broadcast the current cluster centroids to all processes.

    MPI_Bcast(&clustered_points_size, 1, MPI_INT, 0, MPI_COMM_WORLD );
    if(rank > 0){
      cluster_concate.resize(clustered_points_size);
    }
    MPI_Bcast(&cluster_concate[0], cluster_concate.size(), MPI_INT, 0, MPI_COMM_WORLD );
    if(rank > 0){
      cluster_size.resize(5);
    }
    MPI_Bcast(&cluster_size[0], cluster_size.size(), MPI_INT, 0, MPI_COMM_WORLD );

    // Find the closest centroid to each site and assign to cluster.
    

    if(rank > 0){
      int count = 0;
      for(int i = sites_per_proc * (rank - 1); i < sites_per_proc * rank; i++)
			{
        double min_dist;
        int id_cluster_center = 0;
        Point point = test_points[i];

        vector<int> cluster_concate_1 = std::vector<int>(cluster_concate.begin(), cluster_concate.begin() + cluster_size[0]);
        int pointer = cluster_size[0];
        double dist = Cluster::distance_btw_cluster_and_point_only_index(cluster_concate_1, point, test_points);

        min_dist = dist;
        vector<int> cluster_concate_temp;

        for(int j = 1; j < K; j++)
        {
          cluster_concate_temp = std::vector<int>(cluster_concate.begin() + pointer, cluster_concate.begin() + pointer + cluster_size[j]);
          pointer = pointer + cluster_size[j];
          dist = Cluster::distance_btw_cluster_and_point_only_index(cluster_concate_temp, point, test_points);

          if(dist < min_dist)
          {
            min_dist = dist;
            id_cluster_center = j;
          }
        }
				id_new_clusters_per_pro[count] = id_cluster_center;
        count++;
			}
    }
    
    MPI_Gather(&id_new_clusters_per_pro[0], sites_per_proc, MPI_INT,
	     &id_new_clusters[0], sites_per_proc, MPI_INT, 0, MPI_COMM_WORLD);
    

    if (rank == 0) {
      for (int i = 0; i < test_total_points; i++) {
				int id_old_cluster = id_old_clusters[i];
				int id_new_center = id_new_clusters[i];

				if(id_old_cluster != id_new_center)
				{
          if(id_old_cluster == -1){
            total_cluster[id_new_center].push_back(i);
            cout << cluster_size[id_new_center] << endl;
            cluster_size[id_new_center] = cluster_size[id_new_center] + 1;
          }
          // else{
          //   total_cluster[id_new_center].push_back(i);
          //   total_cluster[id_old_cluster].erase(remove(total_cluster[id_old_cluster].begin(), total_cluster[id_old_cluster].end(), i), total_cluster[id_old_cluster].end());
          //   cluster_size[id_new_center] = cluster_size[id_new_center] + 1;
          //   cluster_size[id_old_cluster] = cluster_size[id_old_cluster] - 1;
          // }
					
					done = false;
				}
			}
      // clustered_points_size = accumulate(cluster_size.begin(), cluster_size.end(), 0);
      // cluster_concate.reserve(clustered_points_size);
      // for(int i = 0; i < K; i++){
      //   cluster_concate.insert(cluster_concate.end(), total_cluster[i].begin(), total_cluster[i].end());
      // }
      
    }
  //   // Broadcast the norm.  All processes will use this in the loop test.
  //   MPI_Bcast(&done, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  //   if(done == true || iter >= max_iterations)
  //   {
  //     cout << "Break in iteration " << iter << "\n\n";
  //     break;
  //   }

	// 	iter++;
  // }

  // if(rank == 0){
  //   vector<int> labels;
  //   file.open("doc2vec_labels.txt"); 
  //   if(file.is_open()) {
  //     string line;
  //     while(getline(file, line)) {
  //       int label = stoi(line);
  //       labels.push_back(label);
  //     }
  //     file.close();
  //   }
  //   kmeans.showClusterResults(labels);
  }
      
  MPI_Finalize();

}