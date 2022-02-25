#pragma once
#include <iostream>
#include "HashTable.h"

using namespace std;


/* ================================ HashTable =============================*/
template<typename T>
HashTable<T>::HashTable(int num_buckets, int k, int w, int dimensions, std::string represantation, std::string metric)
{
    this->TableSize = num_buckets;
    this->table = new list<struct Node<T>*>[num_buckets];
    this->k = k;
    this->w = ((w <= 0) ? W : w);
    this->represantation = represantation;
    this->metric = metric;
    this->dimensions = dimensions;

    if (!this->represantation.compare("vector") || !this->metric.compare("continuous")) // for represantation continuous curve and vectors
    {
        this->v = vectors_v(this->k, dimensions);
    }
    else if (!this->represantation.compare("curve") && !this->metric.compare("discrete"))   // for represantation discrete curve
    {
        // get the random vectors v for every hi function 
        // (in discrete frechet we have 2*D coordinates, we flattening the coordinates 
        // e.g: curve: {(x1, y1), (x2, y2)} -- for hashing --> {x1, y1, x2, y2} )
        this->v = vectors_v(this->k, 2*dimensions);   
    }

    this->t = vector_t(this->k);
    this->num_vectors = 0;
}


template<typename T>
HashTable<T>::~HashTable()
{
    for (int i = 0; i < this->TableSize; i++)   //free the table
    {
        while(!table[i].empty()) {
            free(table[i].front());
            table[i].pop_front();
        }
    }
    
    for (auto p : this->v)  //delete random vectors v
    {
        delete p;
    } 

    this->v.clear();

    delete [] table;
}




// calculation of hash function: h = floor[(p v + t) / w]
template<typename T>
int HashTable<T>::h(vector<T>* p, int i)
{

    int result = 0;
    float pv = 0.0;
    float t = this->t[i];
    vector<float>* vect = this->v[i];

    pv = inner_product(begin(*p), std::end(*p), std::begin(*vect), 0.0);
    result = floor((pv + t) / (float)this->w);

    return (int)result;
}


template<typename T>
vector<struct Node<T>*>* HashTable<T>::get_vectors_inBucket(int bucket)
{
    return  new vector<struct Node<T>*>(this->table[bucket].begin(), this->table[bucket].end());
}



/* =============================  HashTableLSH ===========================*/

template <typename T>
HashTableLSH<T>::HashTableLSH(int num_buckets, vector<int> r, int k, int w, int dimensions, double delta, string represantation, string metric) :
 HashTable<T>(num_buckets, k, w, dimensions, represantation, metric)
{
    this->r = r;        // take the random vector r for the calculatation g(p) = [r1h1(p) + r2h2(p) + ... + rkhk (p) mod M] mod TableSize
    this->represantation = represantation;
    this->metric = metric;

    if (!represantation.compare("curve")) 
    {
        this->grid = new Grid<T>(delta);    // grid for curves
    }
}


template <typename T>
HashTableLSH<T>::~HashTableLSH()
{
    if(!this->represantation.compare("curve"))
        delete grid;
}




template <typename T>
int HashTableLSH<T>::insert_vector(struct dataNode<T>* vec)
{
    int hash_value = 0;
    struct Node <T>* node = (struct Node<T>*) malloc(sizeof(struct Node<T>));
    long int id = 0;
    vector<T>* snapped_curve;

    if  (!this->represantation.compare("vector")) {     // lsh with L2 when curve is like vectors
        id = this->ID(vec->vec);
    }
    else if (!this->represantation.compare("curve")) {  // lsh with Frechet curve is a polygynous curve

        if (!this->metric.compare("discrete"))    // curve with discrete  (Aii)
        {
            snapped_curve = this->grid->get_snapped_curve_2D(vec->vec); 

            padding(snapped_curve, 2*this->dimensions);   // padding the snapped curve with the DBL_MAX
           
            id = this->ID(snapped_curve);               // find the ID for the snapped curve
            delete snapped_curve;
        }
        else if (!this->metric.compare("continuous")) // curve with continuous    (Aiii)
        {
            vector<T>* vec_copy = new vector<T>(*(vec->vec));
            filtering(vec_copy);
            snapped_curve = this->grid->get_snapped_curve_1D(vec_copy);
            padding(snapped_curve, this->dimensions);
            id = this->ID(snapped_curve);               // get the right id
            delete vec_copy;                            // delete vec copy, no need any more
        }

    }

    node->ID = id;
    node->vec = vec->vec;
    node->item_id = vec->item_id;

    // calculation of g(p) = [r1h1(p) + r2h2(p) + ... + rkhk (p) mod M] mod TableSize
    // => g(p) = ID % TableSize;
    hash_value = this->hashFunction_g(id);

    this->table[hash_value].push_back(node);

    this->num_vectors++;

    return hash_value;      //return the bucket (we need that for query)
}


// calculation of ID(p) = r1h1(p) + r2h2(p) + ... + rkhk (p) mod M
template<typename T>
long int HashTableLSH<T>::ID(vector<T>* p)
{
    long int M = pow(2, 32) - 5;
    long int result = 0;
    int h = 0;

    for (int i = 0; i < this->k; i++) {
        h = this->h(p, i);
        result += mod(r[i] * h, M);
    }
    result = mod(result, M);

    return result;
}


template <typename T>
long int HashTableLSH<T>::get_ID(string item_id, int bucket)
{
    for (const auto& node: this->table[bucket])
    {
        if (*node->item_id == item_id)
            return node->ID;
    }
    return 0;
}


template <typename T>
std::vector<T>* HashTableLSH<T>::get_snapped_curve_1D(std::vector<T>* curve)
{
    return this->grid->get_snapped_curve_1D(curve);
}  


template <typename T>
std::vector<T>* HashTableLSH<T>::get_snapped_curve_2D(std::vector<T>* curve)
{
    return this->grid->get_snapped_curve_2D(curve);
}  


template<typename T>
void HashTableLSH<T>::print()
{
    struct Node<T>* temp;

    for (int i = 0; i < this->TableSize; i++)
    {
        if (this->table[i].size() == 0)
            continue;
        cout << "bucket: " << i << endl;
        cout << "Number of vectors: " << this->table[i].size() << endl;
        // for (auto x : this->table[i]) {

        //     cout << " ( ID: " << x->ID << ", item_id: " << *(x->item_id) <<" ) \n" ;

        // }
        cout << endl;
    }
}



/* =============================  HashTable Hypercube ===========================*/

template <typename T>
HashTableHC<T>::HashTableHC(int num_buckets, int k, int w, int dimensions) : HashTable<T>(num_buckets, k, w, dimensions, "vector", "")
{
    this->num_vectors = 0;
}

template <typename T>
int HashTableHC<T>::insert_vector(struct dataNode<T>* vec)
{
    int hash_value = 0;
    struct Node <T>* node = (struct Node<T>*) malloc(sizeof(struct Node<T>));
    node->ID = -1;              //ID is for lsh no need for hypercube
    node->vec = vec->vec;
    node->item_id = vec->item_id;
    hash_value = get_hashValue(vec->vec);
    this->table[hash_value].push_back(node);

    this->num_vectors++;

    return hash_value;  //return the bucket which inserted the vector
}


// fill the fi table with 0-1, convert the string of 0-1 to int
template <typename T>
int HashTableHC<T>::get_hashValue(vector<T>* p)
{
    vector<bool> hashValue_binary; //keep the calculation of p -> [f1(h1(p)); : : : ; fd0(hd0(p))]
    int hashValue;
    int h;

    for (int i = 0; i < this->k; i++)
    {
        h = this->h(p, i);
        int fi = check_map(to_string(i), to_string(h));
        hashValue_binary.push_back(fi);
    }

    // convert the binary number to decimal ( this is the bucket that we must save the vector)
    hashValue = accumulate(hashValue_binary.rbegin(), hashValue_binary.rend(), 0, [](int x, int y) { return (x << 1) + y; });

    return hashValue;
}


template<typename T>
int HashTableHC<T>::check_map(string i, string hi)
{
    // the key format: i / hi
    unordered_map<string, int>::iterator got = fi.find (i + "/" + hi); // search in map to find the value with key hi
    bool not_founded = (got == fi.end()); 

    if (not_founded == true) {    // if we dont find it, we calculate the f(hi) -> {0, 1}, also we save it in unordered map
        int f_hi;
        random_device rd;
        mt19937 gen(rd()); // for random numbers fi->{0, 1}
        uniform_int_distribution<int> distrib(0, 1);
        f_hi = distrib(gen);

        fi[i + "/" + hi] = f_hi;        // save for the hi / value the random number {0, 1}    

        return f_hi;
    }
    else {  // else we found the value, so we return the f(hi) - > {0, 1}

        return fi[i + "/" + hi];
    }

    cout << "Something gone wrong in function HashTableHC::check_map\n";
    exit(EXIT_FAILURE);
}



template<typename T>
void HashTableHC<T>::print() //for debuging
{
    struct Node<T>* temp;

    for (int i = 0; i < 2/*this->TableSize*/; i++)
    {
        cout << "bucket: " << i << endl;

        for (auto x : this->table[i]) {
            // for (int i = 0; i < (x->vec)->size(); i++) {
            //     cout << (x->vec)->at(i) << ", "; 
            // }
            cout << "item_id: " << *(x->item_id) <<"\n" ;

            // cout << " ----> ";
        }
        cout << endl;
    }
}


