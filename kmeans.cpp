#include <cstddef>
#include <cstdio>
#include <iterator>
#include <omp.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <ostream>
#include <vector>

using namespace std;

class Point
{
private:
    int pointId, clusterId;
    int dimensions;
    vector<float> values;

    vector<float> bufferToVec(vector<float> &buffer)
    {
        vector<float> values;

        for (int i = 0; i < (int)buffer.size() / 4; i++)
        {
            values.push_back((float)buffer[i]);
        }

        return values;
    }

public:
    Point(int id, vector<float> buffer)
    {
        pointId = id;
        values = bufferToVec(buffer);
        dimensions = values.size();
        clusterId = 0; // Initially not assigned to any cluster
    }

    int getDimensions() { return dimensions; }

    int getCluster() { return clusterId; }

    int getID() { return pointId; }

    void setCluster(int val) { clusterId = val; }

    double getVal(int pos) { return values[pos]; }

    vector<float> getFullVal() { return values; }

    void printFullVal() {
        for (auto x : values) {
            std::cout << x << " ";
        }
        std::cout << std::endl;
    }
};

class Cluster
{
private:
    int clusterId;
    vector<double> centroid;
    vector<Point> points;

public:
    Cluster(int clusterId, Point centroid)
    {
        this->clusterId = clusterId;
        for (int i = 0; i < centroid.getDimensions(); i++)
        {
            this->centroid.push_back(centroid.getVal(i));
        }
        this->addPoint(centroid);
    }

    void addPoint(Point p)
    {
        p.setCluster(this->clusterId);
        points.push_back(p);
    }

    bool removePoint(int pointId)
    {
        int size = points.size();

        for (int i = 0; i < size; i++)
        {
            if (points[i].getID() == pointId)
            {
                points.erase(points.begin() + i);
                return true;
            }
        }
        return false;
    }

    void removeAllPoints() { points.clear(); }

    int getId() { return clusterId; }

    Point getPoint(int pos) { return points[pos]; }

    int getSize() { return points.size(); }

    double getCentroidByPos(int pos) { return centroid[pos]; }

    void setCentroidByPos(int pos, double val) { this->centroid[pos] = val; }
};

class KMeans
{
private:
    int K, iters, dimensions, total_points;
    vector<Cluster> clusters;
    string output_dir;

    void clearClusters()
    {
        for (int i = 0; i < K; i++)
        {
            clusters[i].removeAllPoints();
        }
    }

    int getNearestClusterId(Point point)
    {
        double sum = 0.0, min_dist;
        int NearestClusterId;
        if(dimensions==1) {
            min_dist = abs(clusters[0].getCentroidByPos(0) - point.getVal(0));
        }	
        else 
        {
          for (int i = 0; i < dimensions; i++)
          {
             sum += pow(clusters[0].getCentroidByPos(i) - point.getVal(i), 2.0);
             // sum += abs(clusters[0].getCentroidByPos(i) - point.getVal(i));
          }
          min_dist = sqrt(sum);
        }
        NearestClusterId = clusters[0].getId();

        for (int i = 1; i < K; i++)
        {
            double dist;
            sum = 0.0;
            
            if(dimensions==1) {
                  dist = abs(clusters[i].getCentroidByPos(0) - point.getVal(0));
            } 
            else {
              for (int j = 0; j < dimensions; j++)
              {
                  sum += pow(clusters[i].getCentroidByPos(j) - point.getVal(j), 2.0);
                  // sum += abs(clusters[i].getCentroidByPos(j) - point.getVal(j));
              }

              dist = sqrt(sum);
              // dist = sum;
            }
            if (dist < min_dist)
            {
                min_dist = dist;
                NearestClusterId = clusters[i].getId();
            }
        }

        return NearestClusterId;
    }

public:
    KMeans(int K, int iterations, string output_dir)
    {
        this->K = K;
        this->iters = iterations;
        this->output_dir = output_dir;
    }

    void run(vector<Point> &all_points)
    {
        total_points = all_points.size();
        dimensions = all_points[0].getDimensions();

        // Initializing Clusters
        vector<int> used_pointIds;

        for (int i = 1; i <= K; i++)
        {
            while (true)
            {
                int index = rand() % total_points;

                if (find(used_pointIds.begin(), used_pointIds.end(), index) ==
                    used_pointIds.end())
                {
                    used_pointIds.push_back(index);
                    all_points[index].setCluster(i);
                    Cluster cluster(i, all_points[index]);
                    clusters.push_back(cluster);
                    break;
                }
            }
        }
        cout << "Clusters initialized = " << clusters.size() << endl
             << endl;

        cout << "Running K-Means Clustering.." << endl;

        int iter = 1;
        while (true)
        {
            cout << "Iter - " << iter << "/" << iters << endl;
            bool done = true;

            // Add all points to their nearest cluster
            #pragma omp parallel for reduction(&&: done) num_threads(16)
            for (int i = 0; i < total_points; i++)
            {
                int currentClusterId = all_points[i].getCluster();
                int nearestClusterId = getNearestClusterId(all_points[i]);

                if (currentClusterId != nearestClusterId)
                {
                    all_points[i].setCluster(nearestClusterId);
                    done = false;
                }
            }

            // clear all existing clusters
            clearClusters();

            // reassign points to their new clusters
            for (int i = 0; i < total_points; i++)
            {
                // cluster index is ID-1
                clusters[all_points[i].getCluster() - 1].addPoint(all_points[i]);
            }

            // Recalculating the center of each cluster
            for (int i = 0; i < K; i++)
            {
                int ClusterSize = clusters[i].getSize();

                for (int j = 0; j < dimensions; j++)
                {
                    double sum = 0.0;
                    if (ClusterSize > 0)
                    {
                        #pragma omp parallel for reduction(+: sum) num_threads(16)
                        for (int p = 0; p < ClusterSize; p++)
                        {
                            sum += clusters[i].getPoint(p).getVal(j);
                        }
                        clusters[i].setCentroidByPos(j, sum / ClusterSize);
                    }
                }
            }

            if (done || iter >= iters)
            {
                cout << "Clustering completed in iteration : " << iter << endl
                     << endl;
                break;
            }
            iter++;
        }

        ofstream pointsFile;
        pointsFile.open(output_dir + "/" + to_string(K) + "-points.txt", ios::out);

        for (int i = 0; i < total_points; i++)
        {
            pointsFile << all_points[i].getCluster() << endl;
        }

        pointsFile.close();

        // Write cluster centers to file
        ofstream outfile;
        outfile.open(output_dir + "/" + to_string(K) + "-clusters.txt");
        if (outfile.is_open())
        {
            for (int i = 0; i < K; i++)
            {
                cout << "Cluster " << clusters[i].getId() << " centroid : ";
                for (int j = 0; j < dimensions; j++)
                {
                    cout << clusters[i].getCentroidByPos(j) << " ";    // Output to console
                    outfile << clusters[i].getCentroidByPos(j) << " "; // Output to file
                }
                cout << endl;
                outfile << endl;
            }
            outfile.close();
        }
        else
        {
            cout << "Error: Unable to write to clusters.txt";
        }
    }
};

void load_fvec_data(char* filename, unsigned& num, unsigned& dim, vector<Point>& all_points) { 
    std::ifstream in(filename, std::ios::binary);
    if (!in.is_open()) {
        std::cout << "open file error" << std::endl;
        exit(-1);
    }
    in.read((char*)&dim, 4);	//读取向量维度
    in.seekg(0, std::ios::end);	//光标定位到文件末尾
    std::ios::pos_type ss = in.tellg();	//获取文件大小（多少字节）
    size_t fsize = (size_t)ss;
    num = (unsigned)(fsize / (dim + 1) / 4);	//数据的个数 (ID + dim)
    // char* data;
    
    int pointId = 1;
    in.seekg(0, std::ios::beg);	//光标定位到起始处
    std::cout << dim  << " "  << num << std::endl; 
    std::vector<float> buffer(dim * 4);
    for (size_t i = 0; i < num; i++) {
        
        in.seekg(4, std::ios::cur);	//光标向右移动4个字节
        in.read((char*)buffer.data(),  buffer.size());	//读取数据到一维数据data中
        Point point(pointId, buffer);
        pointId++;
        std::cout << point.getID() << " " << point.getDimensions() << std::endl;
        all_points.push_back(point);
    }

    in.close();
}

int main(int argc, char **argv)
{
    // Need 3 arguments (except filename) to run, else exit
    if (argc != 4)
    {
        cout << "Error: command-line argument count mismatch. \n ./kmeans <INPUT> <K> <OUT-DIR>" << endl;
        return 1;
    }

    string output_dir = argv[3];

    // Fetching number of clusters
    int K = atoi(argv[2]);

    // Open file for fetching points
    string filename = argv[1];
    ifstream infile(filename.c_str());

    if (!infile.is_open())
    {
        cout << "Error: Failed to open file." << endl;
        return 1;
    }

    unsigned points_num, dim;
    vector<Point> all_points;
    load_fvec_data(argv[1], points_num, dim, all_points);
    std::cout << std::endl << "points_num: "<< points_num << std::endl << "data dimension: " << dim << std::endl;

    cout << "\nData fetched successfully!" << endl;

    // Return if number of clusters > number of points
    if ((int)all_points.size() < K)
    {
        cout << "Error: Number of clusters greater than number of points." << endl;
        return 1;
    }

    // Running K-Means Clustering
    int iters = 100;

    KMeans kmeans(K, iters, output_dir);
    kmeans.run(all_points);

    return 0;
}
