#ifndef __Cluster__
#define __Cluster__

#include "Utilities.h"
#include "Dataset.h"
#include "Dataset.cpp"
#include "LSH.h"
#include "LSH.cpp"
#include "Hypercube.h"
#include "Hypercube.cpp"
#include "BinaryTree.h"
#include "BinaryTree.cpp"


#define ITERATIONS 20


template <typename T>
struct node_cluster {
    int centroid;               // centroid belonging to the vectors (using mapping)
    std::vector<struct dataNode<T>*>* vectors;
};


template <typename T>
class Cluster {
    private:
        std::string represantation;
        std::string assignment;
        std::string update_method;

        Dataset<T>* dataset;
        int number_of_clusters;      
        int number_of_vector_hash_tables;
        int number_of_vector_hash_functions;
        int max_number_M_hypercube;
        int number_of_hypercube_dimensions;
        int number_of_probes;

        std::vector<std::vector<T>> centroids;
        std::vector<node_cluster<T>*> clusters;    //keep vectors for every centroid

        LSH<T>* lsh;
        Hypercube<T>* hypercube;
        


    public:
        Cluster(Dataset<T>* data, std::string assignment, std::string represantation, std::string conf_file, std::string update_method);    // in conf file are the arguments
        ~Cluster();
        void get_arguments(std::string file);                // get the arguments from file.conf
        void K_means_plusplus(int number_of_clusters);        // initializing cluster with K-means++
        void K_means_plusplus_2(int number_of_clusters);
        std::vector<std::vector<T>> update_centroids_Mean_Vectors();     // mean of all vectors in cluster
        std::vector<std::vector<T>> update_centroids_Mean_Frechet_2();     // mean of all curve in cluster with frechet
        std::vector<std::vector<T>> update_centroids_Mean_Frechet();     // mean of all curve in cluster with frechet
        void Lloyds();
        void Silhouette(std::ofstream& output_file);
        void Range_Search(std::string method);              // range search based on method (LSH or hypercube)
        void print(std::ofstream& output_file, bool complete);
};




#endif