#include <iostream>
#include "Dataset.h"


using namespace std;


template <typename T>
Dataset<T>::Dataset(string represantation, string metric) 
{
    this->represantation = represantation;
    this->metric = metric;
    this->num_vectors = 0;
    this->dimensions = 0;
}


template<typename T>
Dataset<T>::~Dataset() 
{
    for (int i = 0; i < this->items.size(); i++) {
        delete items[i]->item_id;
        delete items[i]->vec;
        free(items[i]);
    }
}


template<typename T>
void Dataset<T>::insert(string line) 
{
    int max_y, min_y;

    struct dataNode<T>* node = (struct dataNode<T>*)malloc(sizeof(struct dataNode<T>));

    string item_id (line.substr(0, line.find_first_of(" \t"))); // find first of space or tab

    line = line.substr(line.find_first_of(" \t")+1);
    
    node->item_id = new string(item_id);
    _vector* vec = convert_string_to_vector(line);  // get the vector

    if ((dimensions != 0) && (dimensions != vec->d)) {  // check if vector has same dimensions with otehr vectors
        cout << "Warning: input file contains items with different dimensions\n";
    }
    dimensions = vec->d;

    // When we running continuous frechet we must filter the curves to reduce complexity 
    if (!this->represantation.compare("curve") && !this->metric.compare("continuous")) {
        filtering(vec->vec);
    }

    node->vec = vec->vec;

    
    free(vec);
    this->items.push_back(node);        // insert vector in std::vector<vectors>
    this->hash_items[item_id] = node;   // insert vector in unordered_map
    this->num_vectors++;
}


template<typename T>
dataNode<T>* Dataset<T>::get_i_vector(int i)   // use the std::vector
{
    return this->items[i];
}


template<typename T>
dataNode<T>* Dataset<T>::get_vector(string id)  //use the unordered map,  time complexity O(1)
{
    return this->hash_items[id];
}




// with brute force find the N nearest neighbors from query
template <typename T>
vector<vector_distances*> Dataset<T>::findNN_true(vector<T>* vec_query, int number_of_NN)
{
    // cout << "Distance True\n";
    number_of_NN = (number_of_NN <= 0) ? N : number_of_NN;
    
    vector<T> query_vector;
    vector<vector_distances*> NN_vectors;
    vector<vector_distances*> all_vectors;
    struct dataNode<T>* curr_node;
    vector_distances* new_node;


    for(int i = 0; i < this->num_vectors; i++)
    {
        curr_node = items.at(i);

        new_node = (vector_distances*) malloc(sizeof(vector_distances));
        new_node->vec = curr_node->vec;
        new_node->item_id = curr_node->item_id;


        // compute the distance based on represantation and metric
        if (!this->represantation.compare("vector")) {  // if represantation is vector

            new_node->distance_with_query = euclidean_distance(vec_query, new_node->vec);
        
        }
        else if (!this->represantation.compare("curve") && !this->metric.compare("discrete")) { // if represantations is with discrete curves
            vector<T> mean_curve;
            new_node->distance_with_query = Discrete_Frechet_distance(vec_query, new_node->vec, mean_curve);
        }
        else if (!this->represantation.compare("curve") && !this->metric.compare("continuous")) { // if represantations is with continuous curves

            new_node->distance_with_query = Continuous_Frechet_distance(vec_query, new_node->vec);   // we will take the code fromm github for this (but.. how?)

        }

        all_vectors.push_back(new_node);
    }


    // sort vectors based on distance from query
    sort(all_vectors.begin(), all_vectors.end(), [](vector_distances* a, vector_distances* b) {
        return a->distance_with_query < b->distance_with_query;
    });


    for (int i = 0; i < number_of_NN; i++)  //get only the first N Nearest neighbors
    {
        NN_vectors.push_back(all_vectors.at(i));
    }

    for(int i = number_of_NN; i < all_vectors.size(); i++)
    {
        free(all_vectors.at(i));
    }

    return NN_vectors;
}


template<typename T>
int Dataset<T>::calculate_W()   
{
    long int w = 0;
    int num_vectors_for_w = (this->num_vectors < SUBSET ? this->num_vectors : SUBSET);   //take the subset
    vector<vector_distances*> distances;
    vector_distances* new_node = NULL;
    
    // for a point in the dataset calculate the distance between the others points
    // take the avarega distance between them
    // Do it for every point in the subset
    for (int i = 0; i < num_vectors_for_w; i++)
    {
        new_node = (vector_distances*) malloc(sizeof(vector_distances));
        new_node->item_id = items.at(i)->item_id;
        new_node->vec = items.at(i)->vec;
        new_node->distance_with_query = 0;

        for (int j = 0; j < num_vectors_for_w; j++)
        {
            if (*items.at(i)->item_id != *items.at(j)->item_id) 
            {
                new_node->distance_with_query += euclidean_distance(items.at(i)->vec, items.at(j)->vec);
                
            }
        }

        new_node->distance_with_query = new_node->distance_with_query / (num_vectors_for_w-1);   // we have subset-1 distances for each point
        distances.push_back(new_node);
    }


    // summ all distances
    for (int i = 0; i < distances.size(); i++)
    {
        w += distances.at(i)->distance_with_query;
    }

    //find the avarage for all distances
    w = w / num_vectors_for_w;

    for (int i = 0; i < num_vectors_for_w; i++)
    {
        free(distances.at(i));
    }

    return w;
}


//using this struct for k-means++
typedef struct node_P {

    string* vector_id;
    vector<float>* distances;    //all distances between centroids (we need to take the min distance)

}distances_from_centroids;



template <typename T>
bool contain_centroid(vector<struct dataNode<T>*> centroids, string* item)
{
    for (int i = 0; i < centroids.size(); i++) {
        if (*centroids[i]->item_id == *item) {
            // cout << "Centroids\n";
            // for (auto x : centroids)
            //     cout << *x->item_id << endl;
            // cout << "Current centroid: " << *item << endl;
            return true;
        }
    }
    return false;
}


// return k centroids, based on lexture
template <typename T>
vector<vector<T>> Dataset<T>::k_means_plusplus_2(int k, string represantation)   
{   
    random_device rd;
	mt19937 gen(rd()); // for random numbers r 
	uniform_int_distribution<int> distrib(0, this->num_vectors-1);
    int random_centroid = distrib(gen);     // get the first centroid randomly

    vector<vector<T>> centroids_vectors;
    vector<struct dataNode<T>*> centroids;

    vector<distances_from_centroids*> dist_from_allcentroids; //keep all distances from all centroids (we need to take the min distance)
    distances_from_centroids* new_node = NULL;
   

    centroids.push_back(this->items.at(random_centroid));   // here we have the {id, vector}

    float distance = 0;

    vector<T>* centroid = this->items.at(random_centroid)->vec; //get the vector from first random centroid

    // find the distance from all elements from the first centroid keep it in the struct
    for (int i = 0; i < this->num_vectors; i++)
    {
        distance = euclidean_distance(this->items.at(i)->vec, centroid);
        if (!this->represantation.compare("vector")) {
            distance = euclidean_distance(this->items.at(i)->vec, centroid);
        }
        else if (!this->represantation.compare("curve")) {
            vector<T> mean_curve;
            distance = Discrete_Frechet_distance(this->items.at(i)->vec, centroid, mean_curve);
        }
        else {
            cout << "Cann't run k-means++ with represantation: " << this->represantation << " and metric " << this->metric << endl;
            exit(EXIT_FAILURE);
        }
        new_node = (distances_from_centroids*)malloc(sizeof(distances_from_centroids));
        new_node->vector_id = this->items.at(i)->item_id;
        new_node->distances = new vector<float>();
        new_node->distances->push_back(distance);

        dist_from_allcentroids.push_back(new_node);     //keep the {id, <distance1> }
    }



    vector<float> Pi;
    float pi = 0.0, Di = 0.0, x = 0.0;
    float max_Di = 0.0;
    int index = 0;
    bool first_centroid = true;
    int t = 0;      //we have already the first centroid at position 0 in vector : centroids
    

    while(centroids.size() < k)     // keep searching centroids until to find all of them
    {
        if (first_centroid == false) {
            for (int i = 0; i < dist_from_allcentroids.size(); i++) 
            {
                if (contain_centroid(centroids, this->items.at(i)->item_id)) {
                    continue;   // if the point is a centroid skip
                }

                //we calculate the distance from the 't' centroid for every point
                if (!this->represantation.compare("vector")) {
                    distance = euclidean_distance(this->items.at(i)->vec, centroids[t]->vec);
                }
                else if (!this->represantation.compare("curve")) {
                    vector<T> mean_curve;
                    distance = Discrete_Frechet_distance(this->items.at(i)->vec, centroids[t]->vec, mean_curve);
                }
                else {
                    cout << "Cann't run k-means++ with represantation: " << this->represantation << " and metric " << this->metric << endl;
                    exit(EXIT_FAILURE);
                }

                dist_from_allcentroids[i]->distances->push_back(distance);
            }
        }

            
        for (int i = 0; i < dist_from_allcentroids.size(); i++)
        {
            if (contain_centroid(centroids, dist_from_allcentroids.at(i)->vector_id)) {
                continue;   // if the point is a centroid skip
            }

            vector<float> vector_i = *dist_from_allcentroids.at(i)->distances;

            // normalize all D(i)’s by dividing them by maxi D(i)
            max_Di = *max_element(vector_i.begin(), vector_i.end());
            transform(vector_i.begin(), vector_i.end(), vector_i.begin(), [max_Di](float &Di){ return Di/max_Di; });

            Di = *min_element(vector_i.begin(), vector_i.end());    // get the min distance from centroids
            Di = pow(Di, 2);   
            pi += Di;
            Pi.push_back(pi);   // P(r) = Sum_i_r ( D(i) ^ 2 );
        }



        std::random_device rd1;  // Will be used to obtain a seed for the random number engine
        std::mt19937 gen1(rd1()); // Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> dis(0, *max_element(Pi.begin(), Pi.end()));

        x = dis(gen1);      // take a random number in range [0, Pi]


        index = binary(Pi, x, Pi.size());   //get the position of the x if there isn't this nubmer get the smaller number
        if (index > Pi.size()) {
            index = Pi.size() - 1;
        }

        // we dont want to have the same centroids
        while (contain_centroid(centroids, this->items.at(index+1)->item_id))    
            index++;

        // with the previous loop there is a chance to go out of range of the vector so 
        //we should go back until we find the first item wich isn't centroid
        while(index > this->items.size() && contain_centroid(centroids, this->items.at(index+1)->item_id))
            index--;

        centroids.push_back(this->items.at(index+1));     //adding 1 we have the next centroid
        t++;
        Pi.clear();

        first_centroid = false;
    }

    for (int i = 0; i < dist_from_allcentroids.size(); i++)
    {
        new_node = dist_from_allcentroids.at(i);
        delete new_node->distances;
        free(new_node);
    }


    for (int i = 0; i < centroids.size(); i++)
    {
        centroids_vectors.push_back(*centroids.at(i)->vec);
    }

    return centroids_vectors;
}


// based on probability
template <typename T>
vector<vector<T>> Dataset<T>::k_means_plusplus(int k, string represantation)
{
    vector < pair<string, vector<double> > > dist_from_allcentroids;
    vector<pair<string, double> > prob;
    vector<struct dataNode<T>*> vectors = this->items;  // copy dataset
    vector< vector<T> > centroids_vectors;      // centroid as a vectors
    vector<struct dataNode<T>*> centroids;      // centroids as a node from dataset
    vector<double> vector_with_distances;

    vector <double> Pi;
    double pi = 0.0, Di = 0.0, x = 0.0, sum_of_Di = 0.0;
    double max_Di = 0.0;
    double distance = 0.0;


    random_device rd;
	mt19937 gen(rd()); // for random numbers r 
	uniform_int_distribution<int> distrib(0, this->num_vectors-1);
    int random_centroid = distrib(gen);     // get the first centroid randomly

    centroids.push_back(vectors.at(random_centroid));   // here we have the {id, vector}
    vector<T>* curr_centroid = centroids.at(0)->vec; 
    string curr_centroid_id = *centroids.at(0)->item_id;

    // remove from vectors the centroid
    vector<struct dataNode<double>*>::iterator it = vectors.begin();

    while(it != vectors.end())  
    {
        if (!((*it)->item_id->compare(curr_centroid_id)))    
        {   
            vectors.erase(it);
            // delete (*it)->item_id;
            // delete (*it)->vec;
            break;
        }
        it++;
    }



    for (int i = 0; i < vectors.size(); i++)
    {
        // cout << "Compute distances from vector: " << *vectors.at(i)->item_id << " curr_centroid " << curr_centroid_id << endl;
        vector<double> distances;
        if (!this->represantation.compare("vector")) {
            distance = euclidean_distance(vectors.at(i)->vec, curr_centroid);
        }
        else if (!this->represantation.compare("curve")) {
            vector<T> mean_curve;
            distance = Discrete_Frechet_distance(vectors.at(i)->vec, curr_centroid, mean_curve);
        }
        else {
            cout << "Cann't run k-means++ with represantation: " << this->represantation << " and metric " << this->metric << endl;
            exit(EXIT_FAILURE);
        }

        distances.push_back(distance);
        dist_from_allcentroids.push_back( make_pair(*vectors.at(i)->item_id, distances) );
    }

    // cout << "Size of dist from all centroid: " << dist_from_allcentroids.size() << endl;
    bool first_centroid = true;
    int t = 0;
    dataNode<T>* centroid_node;


    while(centroids.size() < k)     // keep searching centroids until to find all of them
    {
        if (first_centroid == false) {  // we have already computed the distance from first centroid
            for (int i = 0; i < dist_from_allcentroids.size(); i++) 
            {
                //we calculate the distance from the 't' centroid for every point
                if (!this->represantation.compare("vector")) {
                    distance = euclidean_distance(vectors.at(i)->vec, curr_centroid);
                }
                if (!this->represantation.compare("curve") && !this->metric.compare("discrete") ) {
                    vector<T> mean_curve;
                    distance = Discrete_Frechet_distance(vectors.at(i)->vec, curr_centroid, mean_curve);
                }
                else {
                    cout << "Cann't run k-means++ with represantation: " << this->represantation << " and metric " << this->metric << endl;
                    exit(EXIT_FAILURE);
                } 
                dist_from_allcentroids.at(i).second.push_back(distance);
            }
        }


        for (int i = 0; i < dist_from_allcentroids.size(); i++)
        {
            vector_with_distances = dist_from_allcentroids.at(i).second;
            max_Di = *max_element(vector_with_distances.begin(), vector_with_distances.end());
        
            transform(vector_with_distances.begin(), vector_with_distances.end(), vector_with_distances.begin(), [max_Di](double &Di){ return Di/max_Di; });

            Di = *min_element(vector_with_distances.begin(), vector_with_distances.end());    // get the min distance from centroids
            Di = pow(Di, 2);   
            pi += Di;
            Pi.push_back(pi);   // P(r) = Sum_i_r ( D(i) ^ 2 );

            prob.push_back( make_pair(dist_from_allcentroids.at(i).first, pi) );
            sum_of_Di += Di;    // sum of min distances
        }

        int highest_probability = -1;
        string next_centroid_id = "";

        // find the next centroid based on proability
        for (auto item : prob) {
            long double prob = item.second / sum_of_Di;
            if (prob > highest_probability) {
                highest_probability = prob;
                next_centroid_id = item.first;
            }
        }

        centroid_node = this->get_vector(next_centroid_id);
        centroids.push_back(centroid_node);

        // remove from vectors the next centroid
        vector<struct dataNode<double>*>::iterator it = vectors.begin();

        while(it != vectors.end())  
        {
            if (!((*it)->item_id->compare(next_centroid_id)))    
            {   
                vectors.erase(it);
                // delete (*it)->item_id;
                // delete (*it)->vec;
                break;
            }
            it++;
        }


        // remove from dist_from_allcentroids the next centroid
        vector < pair<string, vector<double> > >::iterator it2 = dist_from_allcentroids.begin();
        // cout << "Searching centroid: " << curr_centroid_id << endl;
        while(it2 != dist_from_allcentroids.end())  
        {
            if (!((*it2).first.compare(next_centroid_id)))    
            {   
                dist_from_allcentroids.erase(it2);
                // delete (*it)->item_id;
                // delete (*it)->vec;
                break;
            }
            it2++;
        }


        first_centroid = true;
        t++;
    }

    for (auto c : centroids)
        centroids_vectors.push_back(*c->vec);
    

    // cout << "print dataset\n";
    // for (auto c : this->items) {
    //     cout << *(c->item_id) << endl;
    //     for (auto v : *c->vec)
    //         cout << v << " ";
    //     cout << endl;

    // }

    return centroids_vectors;

}





template<typename T>
void Dataset<T>::print(int num_of_col_print) 
{
    if (num_of_col_print > this->items.size())
        num_of_col_print = this->items.size();
    struct dataNode<T>* curr_node;

    for (int i = 0; i < num_of_col_print; i++) {
        curr_node = this->items.at(i);
        cout << *(curr_node->item_id) << " ) ";     // print item_id
        for (int i = 0; i < (curr_node->vec)->size(); i++) {
            cout << curr_node->vec->at(i) << " ";       // print i coordinate
        }
        cout << endl;
    }
}