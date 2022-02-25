#ifndef __UTILITIES__
#define __UTILITIES__

// libraries
#include <math.h>   //sqrt, pow for L2 
#include <vector>
#include <list>
#include <unordered_map>
#include <string>
#include <algorithm>
#include <numeric>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <string.h>
#include <stdio.h>
#include <set>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <sys/stat.h>
#include <iomanip>
#include <cstddef>
#include <limits>
#include <float.h>

// include the necessary files to compute the continuous frechet (i just changed a little bit this repo https://github.com/derohde/Fred )
#include "Continuous_Frechet_github/curve.hpp"
#include "Continuous_Frechet_github/frechet.hpp"


// #define D 128 // dimentions
// #pragma once

// defaults input for LSH
#define W 2
#define DEFAULT_K_LSH 4     // number of hi functions
#define DEFAULT_L 5     // number of hash tables 

// defaults input for Hypercube
#define DEFAULT_K_HYPERCUBE 14       // dimension for the vectors 
#define DEFAULT_M 10       // Maximum number of vectors should be check on hypercube
#define DEFAULT_PROBES 2    // Maximum number of vertices should be check 

#define N 1     // number of nearest neighbors
#define R 10000 //for radius search

// defaults input for cluster
#define NUMBER_OF_VECTOR_HASH_TABLES 3
#define NUMBER_OF_VECTOR_HASH_FUNCTION 4
#define MAX_NUMBER_M_HYPERCUBE 10
#define NUAMBER_OF_HYPERCUBE_DIMENSIONS 3
#define NUMBER_OF_PROBES 2

#define MAX_DIMENSION 200   // max number of coordinates for the centroids, we need this for clustering with frechet discrete

// structs
typedef double Type;

typedef struct _vector {    // for the function convert_string_to_vector
    std::vector<Type>* vec;
    int d;                  //dimensions
}_vector;


template <typename T>   // for hash Tables in lsh
struct Node {
    long int ID;                 
    std::vector<T>* vec;
    std::string* item_id;
};


template<typename T>     // for the dataset
struct dataNode {
    std::string* item_id;
    std::vector<T>* vec;
};


typedef struct vector_distances {   // for find nearest neighbors
    std::vector<Type>* vec;
    std::string* item_id;
    float distance_with_query;
    
}vector_distances;


// functions

// Euclidean distance between 2 vectors
double euclidean_distance(std::vector<Type>* vec1, std::vector<Type>* vec2);  
double euclidean_distance(std::vector<Type>& vec1, std::vector<Type>& vec2);

std::vector<std::vector<float>*> vectors_v(int k, int D);

std::vector<float> vector_t(int k);

_vector* convert_string_to_vector(std::string vec_string);


// make the negative mod positive
long int mod(long int x, long int y);

std::vector<int> vector_r();    // generate the random numbers r for the hash function h

std::string toBinary(int n);    // convert a number to binary, save it as a string

size_t count (const std::string & src, const std::string & str);

//fill the vector with all nearest vertices   equal to given hamming distance
// we need that for hypercube
void get_vertices(char* str, int i, int hamming_distance,  std::vector<std::string>& all_nearest_vertices);

// return the index from the smallest element in array from a given element ( searching with binary search)
// we need that for K-means++
int binary(std::vector<float> arr, float to_search, int arr_size);

// return the index of the second smallest number in vector
int index_of_second_smallest_number(std::vector<float> v);

//get the max number from 2 numbers
float max_number(float x, float y);

// reverse the words in a string e.g str = "Hello world" after reverse str = "world hello" 
std::string reverseString(std::string str);

// check if a file exist
bool exists_file (const std::string& name);


// compute the Discrete Frechet distance between 2 curves (see pg: 9 2.curves.pdf), also fill the vector for optimal traversal
double Discrete_Frechet_distance(std::vector<Type>* curve1, std::vector<Type>* curve2, std::vector<Type>& mean_curve);

// compute Continuous Frechet distance based on repository (https://github.com/derohde/Fred)
double Continuous_Frechet_distance(std::vector<Type>* curve1, std::vector<Type>* curve2);

// padding and create the final vector for the snapped curve
void padding(std::vector<Type>* snapped_curve, int dimensions);       // padding and create the final vector for the snapped curve

// keep the value for padding in stack
double value_for_padding(double M);

// remove the middle point between 3 numbers a < b < c, if |a-b| < e and |b-c| < e
void filtering(std::vector<Type>* curve);

void remove_potential_centroid(std::vector<struct dataNode<Type>*>& potential_centroids, std::string item_id);


#endif /* __UTILITIES__ */
