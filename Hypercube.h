#ifndef __HYPERCUBE___
#define __HYPERCUBE___

#include "HashTable.h"
#include "HashTable.cpp"


template <typename T>
class Hypercube {
    private:
        HashTableHC<T>* hashTable;
        int k;
        int probes;     //maximum number of vertices for check
        int m;     //maximum number of vectors for check
        int num_vectors; 

        // this vector will keep all nearest vectors for query vector
        // this vector will update for each query
        // std::vector<vector_distances> all_vectors_for_query;    
    
    public: 
        Hypercube(int k, int w, int m, int probes, int dimensions);
        ~Hypercube();
        int insert_vector(struct dataNode<T>* vec);    // return the bucket belonging to the vector
        std::vector<vector_distances*> findNN_Hypercube(std::vector<T>* vec_query, int n); // find n nearest neighbors with algorithm Hypercube
        std::vector<struct dataNode<T>*> findNN_range_search(std::vector<T>* vec, float radius);     //find all nearest neighbors inside radius
        std::vector<struct dataNode<T>>* get_M_vectors_from_buckets(std::vector<T>* vec_query);
        int get_num_vectors();
        int get_num_of_h() { return k; }
        inline void set_num_of_h(int k) { this->k = k; }
        int max_num_of_vectors() { return m; }
        int max_num_of_vertices() { return probes; }
        void print();
};

#endif // __HYPERCUBE__
