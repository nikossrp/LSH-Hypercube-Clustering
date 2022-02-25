#ifndef __HashTable__
#define __HashTable__
#include "Utilities.h"
#include "Grid.h"
#include "Grid.cpp"


template<typename T>
class HashTable {
    protected:
        std::list<struct Node<T>*>* table;  // each node  = {item_id, vector, ID}
        int TableSize;      
        int k;                              // number of hi functions and v vectors
        int num_vectors;    
        std::vector<std::vector<float>*> v;  //the random vectors v in function h
        std::vector<float> t;               // the random numbers t for every hi
        int w;
        int dimensions;             // number of dimensions for the vectors
        std::string represantation;
        std::string metric;

        int h(std::vector<T>* p, int i);   // calculation of hash function: h = floor[(p v + t) / w] 

    public:
        HashTable(int buckets, int k, int w, int dimensions, std::string represantation, std::string metric);
        ~HashTable();
        std::vector<struct Node<T>*>* get_vectors_inBucket(int bucket);  //get all vectors in bucket
        int get_num_vectors() { return this->num_vectors; };
};


// hash table for LSH algorithm
template <typename T>
class HashTableLSH: public HashTable<T> {
    private:
        std::vector<int> r;              // the random numbers r in the fromula g(p) = [r1h1(p) + r2h2(p) + ... + rkhk (p) mod M] mod TableSize (is the same for each g)
        Grid<T>* grid;

    public:
        HashTableLSH(int buckets, std::vector<int> r, int k, int w, int dimensions, double delta, std::string represantation, std::string metric);
        ~HashTableLSH();
        int hashFunction_g(long int ID) { return ID % this->TableSize; };
        long int ID(std::vector<T>* p);
        int insert_vector(struct dataNode<T>* vec);
        long int get_ID(std::string item_id, int bucket);   // get the ID from a given item in a given bucket
        std::vector<T>* get_snapped_curve_1D(std::vector<T>* curve);  
        std::vector<T>* get_snapped_curve_2D(std::vector<T>* curve);  

        void print();   //for debuging
};



// hash table for Hypercube algorithm
template <typename T>
class HashTableHC: public HashTable<T> {
    private:
        std::unordered_map<std::string, int> fi;     //map each hi value to {0, 1}
        int num_vectors;
    
    public:
        HashTableHC(int buckets, int k, int w, int dimensions);
        ~HashTableHC() { };
        int insert_vector(struct dataNode<T>* vec);
        int get_hashValue(std::vector<T>* vec);         
        int check_map(std::string i, std::string hi);   // search in fi map the fi(hi) value if not exist create one
        int get_TableSize() { return this->TableSize; }
        int get_num_vectors() { return this->num_vectors; }

        void print(); //for debuging
};










#endif // __HashTable__
