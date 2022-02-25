#ifndef __DATASET__
#define __DATASET__


#include "Utilities.h"

#define SUBSET 100


template<typename T>
class Dataset {

    private:
        int num_vectors;
        std::vector<struct dataNode<T>*> items;     //array of {item_ids, vectors}
        std::unordered_map<std::string, struct dataNode<T>*> hash_items;
        int dimensions;
        std::string metric;
        std::string represantation;

    public:
        Dataset(std::string represantation, std::string metric);
        ~Dataset();
        void insert(std::string  line);
        struct dataNode<T>* get_vector(std::string item_id);
        struct dataNode<T>* get_i_vector(int i);
        int get_num_vectors() { return num_vectors; }
        int get_dimensions() { return dimensions; }
        std::vector<vector_distances*> findNN_true(std::vector<T>* vec_query, int number_of_NN);    //return a vecotr with N items
        std::vector<std::vector<T>>  k_means_plusplus(int k, std::string represantation);
        std::vector<std::vector<T>>  k_means_plusplus_2(int k, std::string represantation);
        int calculate_W();      // calcuate the W from a data subset
        void print(int num_of_col_print);   // debuging


};



#endif // __DATASET__