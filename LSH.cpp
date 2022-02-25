#include <iostream>
#include "LSH.h"

using namespace std;

template <typename T>
LSH<T>::LSH(int buckets, int k, int l, int w, int dimensions, double delta, string representation, string metric): TableSize(buckets)
{
    this->k = ((k <= 0) ? DEFAULT_K_LSH : k);   //check if user gave you varaible k otherwise use the default

    this->l = ((l <= 0) ? DEFAULT_L : l);
    this->w = w;
    this->represantation = representation;
    this->metric = metric;
    this->dimensions = dimensions;

    // Generate random numbers r, k times   (vecot r will be the same for every g function)
    vector<int> r;
    random_device rd;
	mt19937 gen(rd()); // for random numbers r 
	uniform_int_distribution<int> distrib(0, INT32_MAX);
    for (int i = 0; i < this->k; i++)
    {
        r.push_back(distrib(gen));
    }

    for (int i = 0; i < this->l; i++) 
    {
        this->hashTables.push_back(new HashTableLSH<T>(this->TableSize, r, this->k, this->w, dimensions, delta, representation, metric));
    }
}

template<typename T>
LSH<T>::~LSH()
{
    for(auto curr_table: this->hashTables) {
        delete curr_table;
    }
}


//insert the vector at the L hash tables
template<typename T>
vector<int> LSH<T>::insert_vector(struct dataNode<T>* vec)  
{
    HashTableLSH<T>* curr_table;
    vector<int> buckets;
    int bucket = 0;

    for (int i = 0; i < this->hashTables.size(); i++)
    {
        curr_table = this->hashTables[i];
        bucket = curr_table->insert_vector(vec);
        buckets.push_back(bucket);
    }

    return buckets;
}


/*
    1) use the function hashTable->ID, find the ID for query vector
    2) find the bucket which belong the query vector
    2) go to all hash tables in the corresponding buckets and search for vectors with the same ID
    3) if the number of vectors !=  N, search with L2 / Discrete_Frechet / Continuous Frechet
*/
template<typename T>
std::vector<vector_distances*> LSH<T>::findNN_LSH(vector<T>* vec, int n)
{
    n = ((n <= 0) ? N : n);      //check the number of nearest neighbors
    
    int bucket;
    vector<int> buckets;
    HashTableLSH<T>* curr_hashTable;
    vector<long int> query_ID;
    vector<struct Node<T>*>* vectors_in_bucket;
    struct Node<T>* curr_vector;
    vector<struct Node <T>*> vectors_with_same_ID;       // keep only the nearest vectors with query
    vector<struct Node <T>*> all_nearest_vectors;       // keep only the "nearest" vectors from all buckets
    vector<vector_distances> all_vectors_for_query;     //do the same but compute the distances with query 
    vector<vector_distances*> NN_vectors;
    vector<vector_distances*> null;
    vector<T>* snapped_curve;



    //find the bucket for the vector_query in all hash Tables 
    for (int i = 0; i < this->l; i++) 
    {
        long int id; // = this->hashTables[i]->ID(vec);

        if  (!this->represantation.compare("vector")) {     // lsh with L2 when curve is like vector
            id = this->hashTables[i]->ID(vec);
        }
        else if (!this->represantation.compare("curve")) {  // lsh with Frechet curve is a polygynous curve

            if (!this->metric.compare("discrete"))    // curve with discrete  (Aii)
            {
                snapped_curve = this->hashTables[i]->get_snapped_curve_2D(vec); 
                padding(snapped_curve, 2*this->dimensions);     // padding the snapped curve with the max coordinate in dataset
                id = this->hashTables[i]->ID(snapped_curve);    // find the ID for the snapped curve

                delete snapped_curve;
            }
            else if (!this->metric.compare("continuous")) // curve with continuous    (Aiii)
            {
                vector<T>* vec_copy = new vector<T>(*vec);
                // filtering(vec_copy);
                snapped_curve = this->hashTables[i]->get_snapped_curve_1D(vec_copy);
                padding(snapped_curve, this->dimensions);
                id = this->hashTables[i]->ID(snapped_curve); 
                delete vec_copy;
            }

        }

        int hash_value = this->hashTables[i]->hashFunction_g(id);
        buckets.push_back(hash_value);
        query_ID.push_back(id);
    }



    // for each bucket of the corresponding hash table we must find the nearest vectors
    // theoretically nearest vectors is vectors with same ID with query
    for (int i = 0; i < this->l; i++)    
    {
        bucket = buckets[i];    // get the bucket which belong the vector for query
        curr_hashTable = this->hashTables[i];

        vectors_in_bucket = curr_hashTable->get_vectors_inBucket(bucket);

        // fill the vecotr of nearest with query vectors
        for (int j = 0; j < vectors_in_bucket->size(); j++) 
        {
            curr_vector = vectors_in_bucket->at(j);

            all_nearest_vectors.push_back(curr_vector); 

            if ((curr_vector->ID == query_ID[i]))
            {
                vectors_with_same_ID.push_back(curr_vector);
            }
        }

        delete vectors_in_bucket;
    }
    
    struct Node<T>* curr_node;

    if (vectors_with_same_ID.size() >= n)
    {
        // for all vectors in buckets find the dinstace from query vector 
        // only for these vectors which have same id with query
        for (int i = 0; i < vectors_with_same_ID.size(); i++)
        {
            vector_distances new_node;
            curr_node = vectors_with_same_ID.at(i);
            new_node.item_id = curr_node->item_id;
            new_node.vec = curr_node->vec;


            // compute the distance based on represantation and metric
           if (!this->represantation.compare("vector")) {  // if represantation is with vectors

                new_node.distance_with_query = euclidean_distance(vec, new_node.vec);
            
            }
            else if (!this->represantation.compare("curve") && !this->metric.compare("discrete")) { // if represantations is with discrete curves
                vector<T> mean_curve;
                new_node.distance_with_query = Discrete_Frechet_distance(vec, new_node.vec, mean_curve); // function wont compute mean curve, to save time

            }
            else if (!this->represantation.compare("curve") && !this->metric.compare("continuous")) { // if represantations is with continuous curves

                new_node.distance_with_query = Continuous_Frechet_distance(vec, new_node.vec);   // we will take the code fromm github for this (but.. how?)

            }


            all_vectors_for_query.push_back(new_node);
        }
    }
    else {
        // for all vectors in buckets find the dinstace from query vector
        for (int i = 0; i < all_nearest_vectors.size(); i++)
        {
            vector_distances new_node;
            curr_node = all_nearest_vectors.at(i);
            new_node.item_id = curr_node->item_id;
            new_node.vec = curr_node->vec;

            // compute the distance based on represantation and metric
            if (!this->represantation.compare("vector")) {  // if represantation is with vectors

                new_node.distance_with_query = euclidean_distance(vec, new_node.vec);
            
            }
            else if (!this->represantation.compare("curve") && !this->metric.compare("discrete")) { // if represantations is with discrete curves
                vector<T> mean_curve;
                new_node.distance_with_query = Discrete_Frechet_distance(vec, new_node.vec, mean_curve);

            }
            else if (!this->represantation.compare("curve") && !this->metric.compare("continuous")) { // if represantations is with continuous curves

                new_node.distance_with_query = Continuous_Frechet_distance(vec, new_node.vec);   // we will take the code fromm github for this

            }

            all_vectors_for_query.push_back(new_node);
        }
    }


    // cout << "After computed the distance all_vecotrs_for_query = " << all_vectors_for_query.size() << "\n";

    // if we haven't vectors for query just return
    if (all_vectors_for_query.size() == 0) {
        return null;
    }

    // sort based on distance
    sort(all_vectors_for_query.begin(), all_vectors_for_query.end(), [](vector_distances a, vector_distances b) {
        return a.distance_with_query < b.distance_with_query;
    });

    // remove duplicates vectors (in case 2 hashtables has the same vector in the same bucket)
    all_vectors_for_query.erase(std::unique(all_vectors_for_query.begin(), all_vectors_for_query.end(), [](vector_distances a, vector_distances b)
    {
        return !(*a.item_id).compare(*b.item_id);
    }
    ), all_vectors_for_query.end());

    vector_distances* node;

    /* Get only the first n Nearest neighbors from the vector: all_vectors_for_query */
    for (int i = 0; i < n; i++)  
    {
        node = (vector_distances*)malloc(sizeof(vector_distances));
        node->item_id = all_vectors_for_query.at(i).item_id;
        node->distance_with_query = all_vectors_for_query.at(i).distance_with_query;
        node->vec = all_vectors_for_query.at(i).vec;
        NN_vectors.push_back(node);
    }

    return NN_vectors;
}



/* Range search */
// we have already find all vectors in bucket (see all_vectors_for_query), 
//so now we have to see which vector is inside the radius
template<typename T>
vector<struct dataNode<T>*> LSH<T>::findNN_range_search(vector<T>* vec_query, float radius)
{
    vector <struct dataNode<T>*> R_near;
    struct dataNode<T>* node;

    radius = ((radius <= 0) ? R : radius);
    int bucket;
    vector<int> buckets;
    HashTableLSH<T>* curr_hashTable;
    vector<long int> query_ID;
    vector<struct Node<T>*>* vectors_in_bucket;
    struct Node<T>* curr_vector;
    vector<struct Node <T>*> vectors_with_same_ID;       // keep only the nearest vectors with query
    vector<struct Node <T>*> all_nearest_vectors;       // keep only the nearest vectors from all buckets
    vector<vector_distances> all_vectors_for_query;     //do the same but compute the distances with query 
    vector<vector_distances*> NN_vectors;
    long int ID;
    
    
    //check what is the bucket for the vector_query in all hash Tables 
    for (int i = 0; i < this->l; i++) 
    {
        ID = this->hashTables[i]->ID(vec_query);
        if  (!this->represantation.compare("vector")) {     // lsh with L2 when curve is like vectors
            ID = this->ID(vec_query);
        }
        else if (!this->represantation.compare("curve")) {  // lsh with Frechet curve is a polygynous curve
            vector<T>* snapped_curve;
            if (!this->metric.compare("discrete"))    // discrete curve
            {
                snapped_curve = this->hashTables[i]->get_snapped_curve_2D(vec_query); 
                padding(snapped_curve, 2*this->dimensions);   // padding the snapped curve with the max_point in dataset
                ID = this->ID(snapped_curve);               // find the ID for the snapped curve
            }
        }
        int hash_value = this->hashTables[i]->hashFunction_g(ID);
        buckets.push_back(hash_value);
        query_ID.push_back(ID);
    }



    // for each bucket of the corresponding hash table we must find the nearest vectors
    // nearest vectors is vectors with same ID with query
    for (int i = 0; i < this->l; i++)    
    {
        bucket = buckets[i];    //get the bucket which belong the vector for query
        curr_hashTable = this->hashTables[i];

        vectors_in_bucket = curr_hashTable->get_vectors_inBucket(bucket);

        // fill the vecotr of nearest with query vectors
        for (int j = 0; j < vectors_in_bucket->size(); j++) 
        {
            curr_vector = vectors_in_bucket->at(j);

            all_nearest_vectors.push_back(curr_vector);    

            if ((curr_vector->ID == query_ID[i]))
            {
                vectors_with_same_ID.push_back(curr_vector);
            }
        }

        delete vectors_in_bucket;
    }

    struct Node<T>* curr_node;

    // for all vectors in buckets find the dinstace from query vector
    for (int i = 0; i < all_nearest_vectors.size(); i++)
    {
        vector_distances new_node;
        curr_node = all_nearest_vectors.at(i);
        new_node.item_id = curr_node->item_id;
        new_node.vec = curr_node->vec;

        // compute the distance based on represantation and metric
        if (!this->represantation.compare("vector")) {  // if represantation is with vectors

            new_node.distance_with_query = euclidean_distance(vec_query, new_node.vec);
        
        }
        else if (!this->represantation.compare("curve") && !this->metric.compare("discrete")) { // if represantations is with discrete curves
            vector<T> mean_curve;
            new_node.distance_with_query = Discrete_Frechet_distance(vec_query, new_node.vec, mean_curve);

        }
        else if (!this->represantation.compare("curve") && !this->metric.compare("continuous")) { // if represantations is with continuous curves

            new_node.distance_with_query = Continuous_Frechet_distance(vec_query, new_node.vec);   // we will take the code fromm github for this

        }

        all_vectors_for_query.push_back(new_node);
    }

    
    // sort based on distance
    sort(all_vectors_for_query.begin(), all_vectors_for_query.end(), [](vector_distances a, vector_distances b) {
        return a.distance_with_query < b.distance_with_query;
    });

    // remove duplicates vectors (in case 2 hashtables has the same vector in the same bucket)
    all_vectors_for_query.erase(std::unique(all_vectors_for_query.begin(), all_vectors_for_query.end(), [](vector_distances a, vector_distances b) {
        return *a.item_id == *b.item_id;
    }
    ), all_vectors_for_query.end());
    

    for (int i = 0; i < all_vectors_for_query.size(); i++)
    {
        // check if the distance from q is inside the R
        // until we find the first item out of radius (items are sort by distance)
        if (all_vectors_for_query.at(i).distance_with_query <= radius)
        {
            node = (struct dataNode<T>*)malloc(sizeof(struct dataNode<T>));
            node->item_id = all_vectors_for_query.at(i).item_id;
            node->vec = all_vectors_for_query.at(i).vec;
            R_near.push_back(node);
        }
        else 
            break;  
    }

    return R_near;  
}



template <typename T>
vector<struct dataNode<T>>* LSH<T>::get_all_vectors_from_L_buckets(vector<T>* vec_query)
{
    // vector<vector_distances*> R_near;   //return a defferent vector with nearest neighbors in ra
    vector <struct dataNode<T>*> R_near;
    vector_distances* new_node;
    struct dataNode<T>* node;

    int bucket;
    vector<int> buckets;
    HashTableLSH<T>* curr_hashTable;
    vector<long int> query_ID;
    vector<struct Node<T>*>* vectors_in_bucket;
    struct Node<T>* curr_vector;
    vector<struct Node <T>*> vectors_with_same_ID;       // keep only the nearest vectors with query
    vector<struct Node <T>*> all_nearest_vectors;       // keep only the nearest vectors from all buckets
    vector<struct dataNode<T>>* all_vectors_for_query = new vector<struct dataNode<T>>();     //do the same but compute the distances with query 
    vector<vector_distances*> NN_vectors;
    vector<T>* snapped_curve;



    //check what is the bucket for the vector_query in all hash Tables 
    for (int i = 0; i < this->l; i++) 
    {
        long int ID;

        if  (!this->represantation.compare("vector")) {     // lsh with L2 when curve is like vectors
            ID = this->hashTables[i]->ID(vec_query);
        }
        else if (!this->represantation.compare("curve")) {  // lsh with Frechet curve is a polygynous curve

            if (!this->metric.compare("discrete"))    // curve with discrete  (Aii)
            {
                snapped_curve = this->hashTables[i]->get_snapped_curve_2D(vec_query); 

                padding(snapped_curve, 2*this->dimensions);   // padding the snapped curve with the DBL_MAX
            
                ID = this->hashTables[i]->ID(snapped_curve);               // find the ID for the snapped curve

                delete snapped_curve;
            }

        }
        int hash_value = this->hashTables[i]->hashFunction_g(ID);
        buckets.push_back(hash_value);
        query_ID.push_back(ID);
    }

    // for each bucket of the corresponding hash table we must find the nearest vectors
    // nearest vectors is vectors with same ID with query
    for (int i = 0; i < this->l; i++)    
    {
        bucket = buckets[i];    //get the bucket belonging on the vector for query
        curr_hashTable = this->hashTables[i];
        vectors_in_bucket = curr_hashTable->get_vectors_inBucket(bucket);

        // fill the vecotr of nearest with query vectors
        for (int j = 0; j < vectors_in_bucket->size(); j++) 
        {
            curr_vector = vectors_in_bucket->at(j);

            all_nearest_vectors.push_back(curr_vector);    
        }

        delete vectors_in_bucket;
    }

    struct Node<T>* curr_node;

    // for all vectors in buckets find the dinstace from query vector
    for (int i = 0; i < all_nearest_vectors.size(); i++)
    {
        struct dataNode<T> new_node;
        curr_node = all_nearest_vectors.at(i);
        new_node.item_id = curr_node->item_id;
        new_node.vec = curr_node->vec;

        all_vectors_for_query->push_back(new_node);
    }

    
    // remove duplicates vectors (in case 2 hashtables has the same vector in the same bucket)
    all_vectors_for_query->erase(std::unique(all_vectors_for_query->begin(), all_vectors_for_query->end(), [] (struct dataNode<T> a, struct dataNode<T> b) {
        return *a.item_id == *b.item_id;
    }
    ), all_vectors_for_query->end());


    return all_vectors_for_query;
}


template<typename T>
int LSH<T>::get_num_vectors()
{
    int counter = 0;

    for (int i = 0; i < this->hashTables.size(); i++)
    {
        // cout << "In hash table: " << i << " " << this->hashTables.at(i)->get_num_vectors() << endl;
        counter += this->hashTables.at(i)->get_num_vectors();
    }

    return counter;
}



template<typename T>
void LSH<T>::print()
{
    for (int i = 0; i < this->hashTables.size(); i++) {
        cout << "########## Hash Table: " << i << " ################" << endl;
        this->hashTables[i]->print();
    }
}
