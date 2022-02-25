#ifndef __LSH__
#define __LSH__

#include "Utilities.h"
#include "HashTable.h"
#include "HashTable.cpp"



template <typename T>
class LSH {
    private:
        std::vector<HashTableLSH<T>*> hashTables;
        int k;                  // number of hi functions and v vectors
        int l;                  // number of hash tables
        int w;                  
        int TableSize;          // number of buckets for each hash table
        int num_vectors;
        int dimensions;
        std::string represantation;
        std::string metric; 


    public:
        LSH(int buckets, int k, int l, int w, int dimensions, double delta, std::string representation, std::string metric);    // representation can be vector or curve
        ~LSH();
        std::vector<int> insert_vector(struct dataNode<T>* vec);        //insert vector in L hash tables, return the buckets belonging to the vector
        std::vector<struct dataNode<T>>* get_all_vectors_from_L_buckets(std::vector<T>* vector_query);    // get vectors from L buckets for query
        std::vector<vector_distances*> findNN_LSH(std::vector<T>* vec, int n);  
        std::vector<struct dataNode<T>*> findNN_range_search(std::vector<T>* vec_query, float r);
        int get_num_vectors();
        int get_num_of_h() { return k; }
        int get_num_hashTables() { return l; }

        void print();       //for debuging
};

#endif // __LSH__