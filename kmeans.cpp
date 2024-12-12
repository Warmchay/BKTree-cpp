#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <math.h>
#include <omp.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <ostream>
#include <stdlib.h>
#include <vector>
#include <queue>
#include <map>
#include <utility>
#include <vector>

using namespace std;

static omp_lock_t lock;

typedef struct {
    int cluster_id;
    double dist;
} cluster_with_distance;

class Cluster_centroid {
private: 
    int cluster_id_;
    int level_;
    vector<double> cluster_centroid_;
public:
    Cluster_centroid() {
        cluster_id_ = 0;
        level_ = 0;
        cluster_centroid_ = vector<double> {};
    }

    Cluster_centroid(int cluster_id, int level, vector<double> cluster_centroid) {
        this->cluster_id_ = cluster_id;
        this->level_ = level;
        this->cluster_centroid_ = cluster_centroid;
    }

    int get_cluster_id() { return cluster_id_; }
    int get_cluster_level() { return level_; }
    vector<double> get_cluster_centroid() { return cluster_centroid_; }

};

typedef struct {
    int level_;
    vector<Cluster_centroid> cluster_centroid_;
} cluster_centroid_with_level;

struct compare {
    bool operator()(cluster_with_distance a, cluster_with_distance b) {
        return a.dist < b.dist;
    }
};

class Point
{
private:
    int pointId, clusterId, levelID;
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

    void setCluster(int val) { this->clusterId = val;}

    void setLevel(int val) { this->levelID = val; }

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
    int clusterId, levelId;
    vector<double> centroid;
    vector<Point> points;

public:
    Cluster() {
        this->clusterId = 0;
        this->levelId = 0;
        this->points.clear();
    }

    Cluster(int clusterId, Point centroid)
    {
        this->clusterId = clusterId;
        for (int i = 0; i < centroid.getDimensions(); i++)
        {
            this->centroid.push_back(centroid.getVal(i));
        }
        this->addPoint(centroid);
    }

    Cluster(int clusterId, int levelId, Point centroid)
    {
        this->clusterId = clusterId;
        this->levelId = levelId;
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

    KMeans(int K, int iterations) {
        this->K = K;
        this->iters = iterations;
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

// 层内 kmeans
class BKmeans
{
private:
    int level_, dim_, total_points_, k_hortizontal_num_, k_vertical_num_;
    int inner_max_size_, iters_;
    int k_index_;
    // vector<Cluster> clusters_;
    map<int, Cluster> clusters_; // (cluster id, Cluster)
    
public:
    BKmeans(int dim, int level, int total_points, int iterations, int k_vertical_num, int k_index) {
        this->dim_ = dim;
        this->level_ = level;
        this->total_points_ = total_points;
        this->iters_ = iterations;
        this->k_index_ = k_index;
        this->k_hortizontal_num_ = k_index + k_vertical_num;
        this->k_vertical_num_ = k_vertical_num;
        this->inner_max_size_ = ceil(total_points / (double)this->k_vertical_num_);
    }

    void run_bkmeans(vector<Point> &all_points, int level) {
        
        // 1. 选该层所有聚类的质心向量
        vector<int> use_ids;
        cout << "L" << level << " all points size in bkmeans: " << all_points.size() << endl;
        for (int i = this->k_index_; i < this->k_hortizontal_num_; ++i) {
            while (true) {
                int centroid_id = rand() % total_points_;
                if ( find(use_ids.begin(), use_ids.end(), centroid_id) 
                     == use_ids.end() ) 
                {   
                    Point point = all_points[centroid_id];
                    use_ids.push_back(centroid_id);
                    Cluster tmp_cluster{i, level, point};
                    // clusters_[i] = tmp_cluster;
                    if (clusters_.find(i) == clusters_.end()) {
                        clusters_.emplace(i, tmp_cluster);
                    }
                    // clusters_.push_back(tmp_clustler);
                    break;
                }
            }
        }


        // 2. 进行 balanced kmeans
        // 选好了该层的质心向量后，我们需要将这些 points 分配到相应的聚类当中，并且限制数量
        
        for (auto& p : all_points) {
            int best_cluster_id = order_nearest_cluster(p);
            
            p.setCluster(best_cluster_id);

            if (find(use_ids.begin(), use_ids.end(), p.getID()) != use_ids.end()) {
                continue;
            }
            clusters_[best_cluster_id].addPoint(p);
            // cout << best_cluster_id << endl;
        }

        for (auto c : clusters_) {
            cout << "BKmeans Cluster ID(L0): " << c.second.getId() << " Cluster size: " << c.second.getSize() << endl;
        }        
    }

    void run_bkmeans(vector<Point> &all_points, int level, vector<Point> &ori_points) {
        
        // 1. 选该层所有聚类的质心向量
        vector<int> use_ids;
        cout << "L" << level << " all points size in bkmeans: " << all_points.size() << endl;
        for (int i = this->k_index_; i < this->k_hortizontal_num_; ++i) {
            while (true) {
                int centroid_id = rand() % total_points_;
                if ( find(use_ids.begin(), use_ids.end(), centroid_id) 
                     == use_ids.end() ) 
                {   
                    Point point = all_points[centroid_id];
                    use_ids.push_back(centroid_id);
                    Cluster tmp_cluster{i, level, point};
                    if (clusters_.find(i) == clusters_.end()) {
                        clusters_.emplace(i, tmp_cluster);
                    }
                    break;
                }
            }
        }

        // 2. 进行 balanced kmeans
        // 选好了该层的质心向量后，我们需要将这些 points 分配到相应的聚类当中，并且限制数量
        for (auto& p : all_points) {
            int best_cluster_id = order_nearest_cluster(p);
            p.setCluster(best_cluster_id);
            if(find(use_ids.begin(), use_ids.end(), p.getID()) == use_ids.end()) {
                clusters_[best_cluster_id].addPoint(p);
                ori_points[p.getID()].setCluster(best_cluster_id);
            }
        }
        int delete_key = 0;
        bool flag = false;
        for (auto m : clusters_) {
            if (m.first < this->k_index_) {
                delete_key = m.first;
                flag = true;
            }
        }
        if (flag) {
            clusters_.erase(delete_key);
        }
        
        for (auto c : clusters_) {
            cout << "BKMeans Cluster ID: " << c.second.getId() << " Cluster size: " << c.second.getSize() << endl;
        }
    }

    int order_nearest_cluster(Point p) {
        // 返回没有满的 clusters，且距离最近的质心向量
        double min_dist = INT32_MAX, sum = 0.0, dist = 0.0;
        int best_cluster_id = 0;
        for (auto c : clusters_) {
            sum = 0.0;
            if (dim_ == 1) {
                dist = abs(c.second.getCentroidByPos(0) - p.getVal(0));
            } else {
                for (int i = 0; i < dim_; ++i) {
                    sum += pow(c.second.getCentroidByPos(i) - p.getVal(i), 2.0);
                }
                dist = sqrt(sum);
            }
            cluster_with_distance tmp = {c.second.getId(), dist};
            // omp_set_lock(&lock);
            if (c.second.getSize() < inner_max_size_ && dist < min_dist) {
                min_dist = dist;
                best_cluster_id = c.second.getId();
            }
            // omp_unset_lock(&lock);
        }
        return best_cluster_id;
    }

    vector<double> get_cluster_centroid_by_id(int cluster_id) {
        vector<double> ret;
        for (int i = 0; i < dim_; ++i) {
            ret.push_back(clusters_[cluster_id].getCentroidByPos(i));
        }
        return ret;
    }

    int get_vector_num_by_cluster_id(int cluster_id) { return clusters_[cluster_id].getSize(); }
    
};

class HBKmeans
{
private:
    vector<BKmeans> balanced_kmeans_;
    int level_num_, dim_, total_points_;
    int inner_max_size_, iters_, k_vertical_num_;
    vector<cluster_centroid_with_level> cluster_centroid_with_level_; 
    map<int, vector<Cluster_centroid>> cluster_with_level_;

public:
    // 对于用户而言，最方便的 hkbmeans 方式应该是输入每个子聚类应该有的向量数目
    // 但那样不好规定层高和每一层有多少类，也可以算出应该有的向量数目，只需要算最后一层就好
    HBKmeans(int dim, int level_num, int total_points, int iterations, int k_vertical_num) {
        this->dim_ = dim;
        this->level_num_ = level_num;
        this->total_points_ = total_points;
        this->iters_ = iterations;
        this->k_vertical_num_ = k_vertical_num;
    }

    void print_vector(vector<double> vec) {
        cout << "[";
        for (auto v : vec) {
            cout << v << " ";
        }
        cout << "]" << endl;
    }

    void run_hbkmeans(vector<Point> &all_points) {
        vector<Cluster_centroid> cluster_centroid;
        for (int i = 1; i <= level_num_; ++i) {
            std::cout << "Level " << i << " BKMeans start" << std::endl;
            if (i == 1) {
                // 1. 执行 bkmeans
                BKmeans tmp_bkmeans = {dim_, i, total_points_, iters_, k_vertical_num_, 0};
                tmp_bkmeans.run_bkmeans(all_points, i);
                
                // std::cout << "Inject bkmeans centroid start" << std::endl;
                for (int j = 0; j < pow(2, i); ++j) {
                    int cluster_id = j;
                    vector<double> tmp_centroid = tmp_bkmeans.get_cluster_centroid_by_id(j);
                    Cluster_centroid tmp_cluster_centroid = {j, i, tmp_centroid};
                    cluster_with_level_[i].push_back(tmp_cluster_centroid);
                    // cout << "Centroid " << j << ":";
                    // print_vector(tmp_centroid);
                }

            }
            else {
                // 2. 根据上一轮的 cluster id，再次 kmeans 不同 cluster 内的向量
                // 2.1 将 points 根据不同 cluster id 划分好
                map<int, vector<Point>> cluster_with_points;
                for (auto& p : all_points) {
                    cluster_with_points[p.getCluster()].push_back(p);
                }
                // 2.2 训练不同 cluster 内的 points
                int round = 0;
                for (auto m : cluster_with_points) {
                    cout << "Cluster ID:" << m.first << " size: " << m.second.size() << endl;
                    BKmeans tmp_bkmeans = {dim_, i, (int)m.second.size(), iters_, k_vertical_num_, 2 * round};
                    tmp_bkmeans.run_bkmeans(m.second, i, all_points);
                    // 2.3 将不同 cluster 的质心向量记录到 hbkmeans 的队列中
                    for (int j = 2 * round; j < (2 * round + this->k_vertical_num_); ++j) {
                        int cluster_id = j;
                        vector<double> tmp_centroid = tmp_bkmeans.get_cluster_centroid_by_id(j);
                        Cluster_centroid tmp_cluster_centroid = {j, i, tmp_centroid};
                        cluster_with_level_[i].push_back(tmp_cluster_centroid);
                        // cout << "Centroid " << j << ":";
                        // print_vector(tmp_centroid);
                    }
                    round++;
                }
            }
        }
        // 3. 运行结束, 打印层级和每层的 cluster id
        for (auto it = cluster_with_level_.begin(); it != cluster_with_level_.end(); it++) {
            cout << "Level " << it->first << ": ";
            for (auto cluster : it->second) {
                cout << cluster.get_cluster_id() << " ";
            }
            cout << endl;
        }

        cout << "Print each vector belong clusters:" << endl;
        // for (auto p : all_points) {
        //     cout << "Point ID: " << p.getID() << " Cluster ID: " << p.getCluster() << endl;
        // }
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
    for (size_t i = 0; i < 500; i++) {
        
        in.seekg(4, std::ios::cur);	//光标向右移动4个字节
        in.read((char*)buffer.data(),  buffer.size());	//读取数据到一维数据data中
        Point point(pointId, buffer);
        pointId++;
        //std::cout << point.getID() << " " << point.getDimensions() << std::endl;
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

    // KMeans kmeans(K, iters, output_dir);
    // kmeans.run(all_points);

    HBKmeans hbkmeans(dim, 3, all_points.size(), iters, 2);
    hbkmeans.run_hbkmeans(all_points);
    return 0;
}
