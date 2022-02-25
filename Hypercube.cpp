#include <iostream>
#include "Hypercube.h"

using namespace std;

template <typename T>
Hypercube<T>::Hypercube(int k, int w, int m, int probes, int dimensions)
{
    this->m = ((m <= 0) ? DEFAULT_M : m);    
    this->probes = ((probes <= 0)? DEFAULT_PROBES : probes);  // maximun number of vertices to check
    
    this->k = ((k <= 0) ? DEFAULT_K_HYPERCUBE : k);      //if not give k take the default value
    int num_buckets = pow(2, this->k);  // table size

    if (this->probes > num_buckets) 
    {
        cout << "you gave too many vertices. Maximum vertices: " << num_buckets << endl;
        exit(EXIT_FAILURE);
    }

    this->hashTable = new HashTableHC<T>(num_buckets, this->k, w, dimensions);
    this->num_vectors = 0;
}


template <typename T>
Hypercube<T>::~Hypercube()
{
    delete this->hashTable;
}


template <typename T>
int Hypercube<T>::insert_vector(struct dataNode<T>* vec)    // return the buckets belonging to the vector
{
    int bucket = this->hashTable->insert_vector(vec);
    return bucket;
}


template <typename T>
std::vector<vector_distances*> Hypercube<T>::findNN_Hypercube(std::vector<T>* vec_query, int n) 
{
    n = ((n <= 0) ? N : n);   
    int probes_temp = this->probes;
    int m_temp = this->m;   


    vector<struct Node<T>*>* vectors_in_bucket;
    vector<struct Node<T>*> all_vectors_in_buckets;
    std::vector<vector_distances> all_vectors_for_query;    
    vector<int> vertices_checked;
    int vertex = 0;
    struct Node<T>* curr_vector;
    vector<vector_distances*> NN_vectors;
    vector<vector_distances*> null;     // in case we dont find nearest neighbor


    // get the bucket which should go the query
    int bucket = this->hashTable->get_hashValue(vec_query);


    vectors_in_bucket = this->hashTable->get_vectors_inBucket(bucket);


    for (int i = 0; i < vectors_in_bucket->size(); i++) //gather all vectors in a std::vector
    {
        curr_vector = vectors_in_bucket->at(i);
        if (all_vectors_in_buckets.size() >= m_temp)     //dont superpass the limit of check vectors
            break;
        all_vectors_in_buckets.push_back(curr_vector);
    }
    delete vectors_in_bucket;

    vector<string> all_nearest_vertices;    // here we will keep all nearest vertices
    string s_str = toBinary(bucket);

    while(s_str.length() < this->k)
        s_str.insert(0, "0");

    char s_c [s_str.length() + 1];

    // copying the contents of the
    // string to char array, we will use the function magic
    strcpy(s_c, s_str.c_str());

    int hamming_distance = 1;   // starting with hamming distance = 1

    // if we dont have already check the maximum vectors go to the nearest vertices until consume all m or probes
    while((all_vectors_in_buckets.size() < m_temp ) && (probes_temp > 0))  
    {
        get_vertices(s_c, strlen(s_c)-1, hamming_distance++, all_nearest_vertices);
        all_nearest_vertices.erase(all_nearest_vertices.begin());   //remove the first bucket, it is the query bucket there
        
        for (int i = 0; i < all_nearest_vertices.size(); i++)   // go to each vertex
        {
            if (probes_temp <= 0) //dont superpass the maximum vertex
                break; 

            string bucket =  all_nearest_vertices[i];
            vertex = stoi(bucket, 0, 2);
            vectors_in_bucket = this->hashTable->get_vectors_inBucket(vertex);  //for every vertex get the vectors
            for (int j = 0; j < vectors_in_bucket->size(); j++) 
            {
                if (all_vectors_in_buckets.size() >= m )    //dont superpass the maximum number of points for check
                    break;
                curr_vector = vectors_in_bucket->at(j);
                all_vectors_in_buckets.push_back(curr_vector);      // insert all vectors in keeper
            }
            probes_temp--;       // we checked a vertex
            delete vectors_in_bucket;
        }
        all_nearest_vertices.clear();   // go to the next vertex
    }




    struct Node<T>* curr_node;
    // for all vectors in buckets find the dinstace from query vector
    for (int i = 0; i < all_vectors_in_buckets.size(); i++)
    {
        vector_distances new_node;
        curr_node = all_vectors_in_buckets.at(i);
        new_node.item_id = curr_node->item_id;
        new_node.vec = curr_node->vec;
        new_node.distance_with_query = euclidean_distance(vec_query, new_node.vec);
        all_vectors_for_query.push_back(new_node);
    }

    // sort the vectors based on distance between query
    sort(all_vectors_for_query.begin(), all_vectors_for_query.end(), [](vector_distances a, vector_distances b) {
        return a.distance_with_query < b.distance_with_query;
    });

    // if we haven't vectors for query just return
    if (all_vectors_for_query.size() == 0) {
        return null;
    }

    vector_distances* node;

    if (all_vectors_for_query.size() < n)   //Take whatever you have already find
        n = all_vectors_for_query.size();

    /* Get only the first n Nearest neighbors */
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


// we have already find all vectors < M  so now we have to see which vector is inside radius
template <typename T>
std::vector<struct dataNode<T>*> Hypercube<T>::findNN_range_search(vector<T>* vec_query, float radius)     //find all nearest neighbors inside radius
{

    vector<struct dataNode<T>*> R_near;   //return a defferent vector with nearest neighbors in ra
    struct dataNode<T>* new_node;

    radius = ((radius <= 0) ? R : radius);

    int probes_temp = this->probes;
    int m_temp = this->m;   


    struct Node<T>* curr_vector;

    vector<struct Node<T>*>* vectors_in_bucket;
    vector<struct Node<T>*> all_vectors_in_buckets;
    std::vector<vector_distances>* all_vectors_for_query = new vector<vector_distances>();    
    vector<int> vertices_checked;
    int vertex = 0;
    vector<vector_distances*> NN_vectors;


    // get the bucket which should go the query
    int bucket = this->hashTable->get_hashValue(vec_query);
    
    // get the vectors in this bucket
    vectors_in_bucket = this->hashTable->get_vectors_inBucket(bucket);

    for (int i = 0; i < vectors_in_bucket->size(); i++) //gather all vectors in a std::vector
    {
        curr_vector = vectors_in_bucket->at(i);
        if (all_vectors_in_buckets.size() >= m_temp)     //dont superpass the limit of check vectors
            break;
        all_vectors_in_buckets.push_back(curr_vector);
    }

    delete vectors_in_bucket;
    vector<string> all_nearest_vertices;
    string s_str = toBinary(bucket);

    while(s_str.length() < this->k) // make the string length equal to k
        s_str.insert(0, "0");

    char s_c [s_str.length() + 1];
    strcpy(s_c, s_str.c_str());

    // find M vectors
    int hamming_distance = 1;

    // if we dont have already check the maximum number of vectors (M) go to the nearest vertices until consume all m or probes
    while((all_vectors_in_buckets.size() < m_temp ) && (probes_temp > 0))  
    {
        get_vertices(s_c, strlen(s_c)-1, hamming_distance++, all_nearest_vertices);
        all_nearest_vertices.erase(all_nearest_vertices.begin());   //remove the first bucket, it is the query bucket there
        

        for (int i = 0; i < all_nearest_vertices.size(); i++)   // go to each vertex
        {
            if (probes_temp <= 0) //dont superpass the maximum number of vertices
                break; 

            string bucket =  all_nearest_vertices[i];
            vertex = stoi(bucket, 0, 2);
            vectors_in_bucket = this->hashTable->get_vectors_inBucket(vertex);  //for every vertex get the vectors
            for (int j = 0; j < vectors_in_bucket->size(); j++) 
            {
                if (all_vectors_in_buckets.size() >= m )    //dont superpass the maximum number of points for check
                    break;
                curr_vector = vectors_in_bucket->at(j);
                all_vectors_in_buckets.push_back(curr_vector);      // insert all vectors in keeper
            }
            probes_temp--;       // we checked a vertex
            delete vectors_in_bucket;
        }
        all_nearest_vertices.clear();   // go to the next vertex
    }


    struct Node<T>* curr_node;

    // for all vectors in buckets find the dinstace from query vector
    for (int i = 0; i < all_vectors_in_buckets.size(); i++)
    {
        vector_distances new_node;
        curr_node = all_vectors_in_buckets.at(i);
        new_node.item_id = curr_node->item_id;
        new_node.vec = curr_node->vec;
        new_node.distance_with_query = euclidean_distance(vec_query, new_node.vec);
        all_vectors_for_query->push_back(new_node);
    }

    // sort the vectors based on distance between query
    sort(all_vectors_for_query->begin(), all_vectors_for_query->end(), [] (vector_distances a, vector_distances b) {
        return a.distance_with_query < b.distance_with_query;
    });

    for (int i = 0; i < all_vectors_for_query->size(); i++)
    {
        // check if the distance from q is inside the R
        //untill we find a item out of radius (items are sorted by distance)
        if (all_vectors_for_query->at(i).distance_with_query <= radius)
        {
            new_node = (struct dataNode<T>*)malloc(sizeof(struct dataNode<T>));
            new_node->item_id = all_vectors_for_query->at(i).item_id;
            new_node->vec = all_vectors_for_query->at(i).vec;
            R_near.push_back(new_node);
        }
        else 
            break;
    }
    delete all_vectors_for_query;

    return R_near;  
}



template<typename T>
vector<struct dataNode<T>>* Hypercube<T>::get_M_vectors_from_buckets(vector<T>* vec_query)
{
    int probes_temp = this->probes;
    int m_temp = this->m;   


    vector<struct Node<T>*>* vectors_in_bucket;
    vector<struct Node<T>*> all_vectors_in_buckets;
    std::vector<struct dataNode<T>>* all_vectors_for_query = new vector<struct dataNode<T>>();    
    vector<int> vertices_checked;
    int vertex = 0;
    struct Node<T>* curr_vector;


    // get the bucket which should go the query
    int bucket = this->hashTable->get_hashValue(vec_query);
    
    // get the vectors in this bucket
    vectors_in_bucket = this->hashTable->get_vectors_inBucket(bucket);

    for (int i = 0; i < vectors_in_bucket->size(); i++) //gather all vectors in a std::vector
    {
        curr_vector = vectors_in_bucket->at(i);
        if (all_vectors_in_buckets.size() >= m_temp)     //dont superpass the limit of check vectors
            break;
        all_vectors_in_buckets.push_back(curr_vector);
    }

    delete vectors_in_bucket;
    vector<string> all_nearest_vertices;
    string s_str = toBinary(bucket);

    while(s_str.length() < this->k) // make the string length equal to k
        s_str.insert(0, "0");

    char s_c [s_str.length() + 1];
    strcpy(s_c, s_str.c_str());

    // find M vectors
    int hamming_distance = 1;

    // if we dont have already check the maximum number of vectors (M) go to the nearest vertices until consume all m or probes
    while((all_vectors_in_buckets.size() < m_temp ) && (probes_temp > 0))  
    {
        get_vertices(s_c, strlen(s_c)-1, hamming_distance++, all_nearest_vertices);
        all_nearest_vertices.erase(all_nearest_vertices.begin());   //remove the first bucket, it is the query bucket there
        

        for (int i = 0; i < all_nearest_vertices.size(); i++)   // go to each vertex
        {
            if (probes_temp <= 0) //dont superpass the maximum number of vertices
                break; 

            string bucket =  all_nearest_vertices[i];
            vertex = stoi(bucket, 0, 2);
            vectors_in_bucket = this->hashTable->get_vectors_inBucket(vertex);  //for every vertex get the vectors
            for (int j = 0; j < vectors_in_bucket->size(); j++) 
            {
                if (all_vectors_in_buckets.size() >= m )    //dont superpass the maximum number of points for check
                    break;
                curr_vector = vectors_in_bucket->at(j);
                all_vectors_in_buckets.push_back(curr_vector);      // insert all vectors in keeper
            }
            probes_temp--;       // we checked a vertex
            delete vectors_in_bucket;
        }
        all_nearest_vertices.clear();   // go to the next vertex
    }


    struct Node<T>* curr_node;

    // for all vectors in buckets find the dinstace from query vector
    for (int i = 0; i < all_vectors_in_buckets.size(); i++)
    {
        struct dataNode<T> new_node;
        curr_node = all_vectors_in_buckets.at(i);
        new_node.item_id = curr_node->item_id;
        new_node.vec = curr_node->vec;
        all_vectors_for_query->push_back(new_node);
    }

    return all_vectors_for_query;
}



template<typename T>
void Hypercube<T>::print()
{
    this->hashTable->print();
}


template<typename T>
int Hypercube<T>::get_num_vectors()
{
    return this->hashTable->get_num_vectors();
}

