#include <iostream>
#include "Cluster.h"

using namespace std;

template <typename T>
Cluster<T>::Cluster(Dataset<T>* data, string assignment, string represantation, string conf_file, string update_method)
{
    int w = 0, D = 0, k_lsh = 0, k_cube = 0, l = 0, tableSize, m = 0, probes = 0;
    this->dataset = data;
    this->represantation = represantation;
    this->update_method = update_method;
    this->assignment = assignment;

    double delta = 1.9;

    //initialize all values
    this->get_arguments(conf_file);

    //initialize the centroids
    this->K_means_plusplus(this->number_of_clusters);

    // compute the number of buckets for hash tables based on number of vectors in dataset
    tableSize = ceil(data->get_num_vectors() / 4.0); 

    w = this->dataset->calculate_W()/2;     // calculate the w based on a data subset

    D = dataset->get_dimensions();
    cout << "Dimensions of curves/vectors: " << D << endl;
    cout << "Represantation: " << this->represantation << endl;
    if (this->assignment != "Classic") {
        cout << "W: " << w << endl;
        cout << "Table size: n/4 = " << tableSize << endl;

        if (!this->represantation.compare("curve")) {
            cout << "Delta: " << delta << endl;
        }
    }

    k_lsh = this->number_of_vector_hash_functions;
    l = this->number_of_vector_hash_tables;
    m = this->max_number_M_hypercube;
    probes = this->number_of_probes;
    k_cube - this->number_of_hypercube_dimensions;
        
    if (!assignment.compare("LSH") || !assignment.compare("LSH_Frechet")) {
        this->lsh = new LSH<T>(tableSize, k_lsh, l, w, D, delta, represantation, "discrete");
    }
    else {
        this->lsh = NULL;
    }
    if (!assignment.compare("Hypercube")) {
        this->hypercube = new Hypercube<T>(k_cube, w, m, probes, D);
    }
    else {
        this->hypercube = NULL;
    }

    dataNode<T>* vec;

    if (this->assignment != "Classic") {
        // insert vectors in Hash Tables (LSH) and hypercube
        for (int i = 0; i < dataset->get_num_vectors(); i++)
        {
            vec = dataset->get_i_vector(i);     //get a pointer from dataset
            if (!this->assignment.compare("LSH") || !this->assignment.compare("LSH_Frechet"))
                this->lsh->insert_vector(vec);            //insert pointer to lsh structs
            if (!this->assignment.compare("Hypercube"))
                this->hypercube->insert_vector(vec);
        }
    }
}



template<typename T>
Cluster<T>::~Cluster()
{
    node_cluster<T>* node;

    for (int i = 0; i < this->clusters.size(); i++)
    {
        node = this->clusters.at(i);
        delete node->vectors;
        free(node);
    }

    if (!this->assignment.compare("LSH") || !this->assignment.compare("LSH_Frechet"))
        delete lsh;
    else if (!this->assignment.compare("Hypercube"))
        delete hypercube;
}



template <typename T>
void Cluster<T>::K_means_plusplus(int number_of_clusters)
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::seconds;
    using std::chrono::milliseconds;

    auto t_start = high_resolution_clock::now();
    this->centroids = dataset->k_means_plusplus(number_of_clusters, this->represantation);   // get the initialized centroids
    auto t_end = high_resolution_clock::now(); 

    duration<double, std::milli> t = t_end - t_start;

    printf("K-Means++ time : %.2f ms\n", t.count());
}






template<typename T>
void Cluster<T>::get_arguments(string file)
{
    string line = "";
    // get the arguments for cluster from the file cluster.conf, initialize all parametres.
    ifstream cluster_file(file);
    if (cluster_file.is_open()) 
    {
        size_t pos_punctuation;
        string argument = "";

        while(getline(cluster_file, line)) {
            pos_punctuation = line.find(":");
            argument = line.substr(0, pos_punctuation);

            if (argument == "number_of_clusters") {
                number_of_clusters = stoi(line.substr(pos_punctuation+1));
            }
            else if (argument == "number_of_vector_hash_tables") {
                number_of_vector_hash_tables = stoi(line.substr(pos_punctuation+1));
            }
            else if (argument == "number_of_vector_hash_functions") {
                number_of_vector_hash_functions = stoi(line.substr(pos_punctuation+1));
            }
            else if (argument == "max_number_M_hypercube") {
                max_number_M_hypercube = stoi(line.substr(pos_punctuation+1));
            }
            else if (argument == "number_of_hypercube_dimensions") {
                number_of_hypercube_dimensions =  stoi(line.substr(pos_punctuation+1)); 
            }
            else if (argument == "number_of_probes") {
                number_of_probes =  stoi(line.substr(pos_punctuation+1));
            }
        }
    }
    else 
        cout << "Uable to open cluster.conf file";


    // check if arguments is  greater that 0 if it doesn't use the default values
    number_of_vector_hash_tables = (number_of_vector_hash_tables <= 0) ? NUMBER_OF_VECTOR_HASH_TABLES : number_of_vector_hash_tables;
    number_of_vector_hash_functions = (number_of_vector_hash_functions <= 0) ? NUMBER_OF_VECTOR_HASH_FUNCTION : number_of_vector_hash_functions;
    max_number_M_hypercube = (max_number_M_hypercube <= 0) ? MAX_NUMBER_M_HYPERCUBE : max_number_M_hypercube;
    number_of_hypercube_dimensions = (number_of_hypercube_dimensions <= 0) ? NUAMBER_OF_HYPERCUBE_DIMENSIONS : number_of_hypercube_dimensions;
    number_of_probes = (number_of_probes <= 0) ? NUMBER_OF_PROBES : number_of_probes; 
}




// With Binary tree
template <typename T>
vector<vector<T> > Cluster<T>::update_centroids_Mean_Frechet()     // mean of all curve in cluster with frechet
{

    vector<vector<T>> new_centroids;
    int num_curves = 0;
    

    // compute the mean of each cluster (immplement with vector)
    for (int index_cluster = 0; index_cluster < this->centroids.size(); index_cluster++)
    {
        vector<T> mean_curve;
        node_cluster<T>* curr_cluster = this->clusters.at(index_cluster); 
        if ((curr_cluster->vectors->size() <= 1)) {
            new_centroids.push_back(centroids.at(curr_cluster->centroid));
            continue;
        }

        num_curves = curr_cluster->vectors->size();
        BinaryTree<T>* binaryTree = new BinaryTree<T>(num_curves);  // create the binary tree for calculate the mean curve

        // get all curves for the current cluster inside to binary tree
        for (int i = 0; i < num_curves; i++) {
            binaryTree->insert_curve(curr_cluster->vectors->at(i)->vec);
        }

        binaryTree->get_mean_curve(mean_curve);

        new_centroids.push_back(mean_curve);

        delete binaryTree;
    }

    return new_centroids;
}





// calculate the means of vectors in each cluster, so we have the new centroids
template <typename T>
vector<vector<T>> Cluster<T>::update_centroids_Mean_Vectors()
{
    vector<vector<T>> new_centroids;
    vector<T> new_centroid;
    node_cluster<T>* cluster = NULL;
    int number_of_vectors = 0;

    // calculate the new centroids
    for (int i = 0; i < this->centroids.size(); i++)
    {
        cluster = this->clusters.at(i); // take the i cluster


        //initialize the new centroid
        for (int d = 0; d < dataset->get_dimensions(); d++){
            new_centroid.push_back(0.0);
        }

        if (cluster->vectors->size() == 0) {
            new_centroids.emplace_back(this->centroids.at(i));
            new_centroid.clear();   //go to the next centroid
            continue;
        }

        // add all coordinates from all vectors in cluster, so we make a new vector
        for (int j = 0; j < cluster->vectors->size(); j++)
        {
            transform(new_centroid.begin(), new_centroid.end(), cluster->vectors->at(j)->vec->begin(), new_centroid.begin(), std::plus<T>());
        }

        number_of_vectors = cluster->vectors->size();

        // divide all coordinates with the number of vectors in cluster, so we have a new centroid
        std::for_each(new_centroid.begin(), new_centroid.end(), [number_of_vectors](T &n){ n /= number_of_vectors; });

        new_centroids.emplace_back(new_centroid);

        new_centroid.clear();   //go to the next centroid
    }

    return new_centroids;
}



// runt the Lloyds algorithm, print the expected output
template <typename T>
void Cluster<T>::Lloyds()
{
    struct dataNode<T> *curr_vector , *new_vector;
    node_cluster<T>* cluster = NULL;
    vector<vector<T>> new_centroids;
    vector<bool> flags;         // flags for checking the centroids ( if all flags is true we have find the optimal centroids)
    int count_flags = 0;        // if we find at least > k/2 false we stop the process
    int iterations = ITERATIONS;

    vector<float> distances;        // for every centroid keep the distance between point
    float distance = 0.0;
    int centroid_index = 0;
    bool already_exist = false;
    int i = 0;  
    while(1)      // if you dont have find the optimal centroid after 100 iterations stop the algorithm
    {
        count_flags = 0;

        // for every point we must calculate the distance from all centroids
        for (int i = 0; i < this->dataset->get_num_vectors(); i++)
        {
            already_exist = false;

            curr_vector = dataset->get_i_vector(i);
            new_vector = curr_vector;

            //calculate the distance between vector from all centroids
            for (int j = 0; j < this->centroids.size(); j++) 
            {
                if (!this->represantation.compare("vector")) {
                    distance = euclidean_distance(curr_vector->vec, &this->centroids.at(j));
                }
                else if (!this->represantation.compare("curve")) {
                    vector<T> mean_curve;
                    distance = Discrete_Frechet_distance(curr_vector->vec, &this->centroids.at(j), mean_curve);
                }

                distances.push_back(distance);
            }   

            // keep only the minimum distance 
            //(we will find the index and we will take the centroid with this index)
            centroid_index = min_element(distances.begin(), distances.end()) - distances.begin();

            // Assign every point to centroid
            for (int j = 0; j < this->clusters.size(); j++)
            {
                cluster = this->clusters.at(j);

                if (cluster->centroid == centroid_index) {  // check if cluster already exist
                    cluster->vectors->push_back(new_vector);
                    already_exist = true;
                }
            }

            if (already_exist == false) { //match the vector with centroid if we dont have already find this centroid
                cluster = (node_cluster<T>*) malloc(sizeof(node_cluster<T>));
                cluster->centroid = centroid_index;
                cluster->vectors = new vector<struct dataNode<T>*>();
                cluster->vectors->push_back(new_vector);   
                clusters.push_back(cluster);
            }

            distances.clear();
        }


        // cout << "Update centroids via " << update_method << endl;
        if (!this->update_method.compare("Mean Vector")) {
            new_centroids = this->update_centroids_Mean_Vectors();
        }
        else if (!this->update_method.compare("Mean Frechet")) {
            new_centroids = this->update_centroids_Mean_Frechet();
        }

        // check the distance between new centroids and the old centroids
        // if distance between old_centroids and new_centroids
        // is smaller that 1, keep going until you find a better centroid
        for (int i = 0; i < this->centroids.size(); i++)
        {
            float distance_bw_centroids = 0.0;
            if (!this->represantation.compare("vector")) {
                distance_bw_centroids = euclidean_distance(&this->centroids.at(i), &new_centroids.at(i));
            }
            else if (!this->represantation.compare("curve")) {
                vector<T> mean_curve;
                distance_bw_centroids = Discrete_Frechet_distance(&this->centroids.at(i), &new_centroids.at(i), mean_curve);
            }

            // cout << "Distance between old and new centroids: " << distance_bw_centroids << endl;
            if ((distance_bw_centroids) <= 0.5)    // we have the optimal centroid
            {
                flags.push_back(true);
            }
            else {      // we need to find another centroid
                flags.push_back(false);
            }
        }

        for(auto b : flags)
        {
            if ( b == true)
                count_flags++;
        }


        // if at least the half centroids are optimal stop the algoritm
        if ((count_flags >= this->number_of_clusters/2) || iterations <= 0)
        {
            break;
        }
        else {  //clear current cluster and go to the next
            // this->print();
            for (int i = 0; i < this->clusters.size(); i++)
            {
                cluster = this->clusters.at(i);
                delete cluster->vectors;
                free(cluster);
            }
            this->clusters.clear();
            count_flags = 0;
        }
        
        // swap the old centroids with new one
        this->centroids.swap(new_centroids);

        new_centroids.clear();  // here is the olds centroids after swap so we have to clear
        flags.clear();
        iterations--;
    }
}





template <typename T>
void Cluster<T>::Range_Search(string method)    // method can be LSH or Hypercube LSH_Frechet
{
    if (this->number_of_clusters < 2) {
        cout << "Too less number of clusters";
        exit(EXIT_FAILURE);
    }

    std::chrono::steady_clock::time_point begin, end;
    double time_for_brute_force = 0.0;

    unordered_map<string, bool> flags_vectors;    // when a vector belong in a cluster the flag is true, otherwise is false
    unordered_map<string, int> assigned_vectors;        // keep all teh assigments
    vector<bool> flags_centroids;                // flags for checking the centroids ( if all flags is true we have find the optimal centroids)
    string item_id = "";
    long int number_of_vectors = dataset->get_num_vectors();
    node_cluster<T> *curr_cluster = NULL, *curr_assigned_cluster = NULL;
    struct dataNode<T>* curr_vector;       // current vector  with distance between query and vector
    struct dataNode<T>* curr_node;      // current vector from struct of dataset
    vector<double> distances, distances_radius;        // for every centroid keep the distance between point
    vector<vector<T>> new_centroids;
    vector<T> *curr_centroid, *next_centroid;
    vector<T> new_centroid;
    vector<vector<struct dataNode<T>>*> vectors_in_buckets_for_all_centroids;     // for every centroid keep the vectors from L buckets
    vector<struct dataNode<T>>* vectors_in_buckets;     // for every centroid keep the vectors from L buckets
    vector<float> radius;
    double curr_radius = 0.0;
    struct dataNode<T>* new_vector;
    vector<T> *curr_assigned_centroid;
    int curr_assigned_centroid_index = 0;
    double distance = 0.0, min_distance;
    int min_distance_index;
    int unassigned_vectors = number_of_vectors;
    int count_flags = 0;
    double min_radius = 0.0 ,min_radius_temp = 0.0;
    int tableSize;
    bool overflow_radius = false;
    int end_point = ceil((float)number_of_vectors * 5/100.0);     // when 10% of n are unassigned start the brute force (Lloyd's)
    int centroid_index, centroid_index2;
    int all_iterations = 1;
    bool centroids_changed = true;
    int vectors_in_cluster = 0;
    int inserts_per_iterations = 0;
    vector<int> insertions_per_centroids;       // keep the number of insertions for every centroid
    int r = 1;
    int loop = 0;
    bool initialized_cluster = false;
    int points_which_assignned_with_brute_force = 0;
    int points_which_assignned_with_method = 0;
    bool collision_centroids = false;
    int stable_clusters = 0;
    int iterations = ITERATIONS;

    //each centroids insert in lsh like query and executed range search 
    while(1)
    {   

        /* for each centroid initialize the radius , save them in a vector*/
        // also get all vectors from L/probes buckets
        if (centroids_changed == true)  
        {

            for (int i = 0; i < this->centroids.size(); i++)
            {
                curr_centroid = &(this->centroids.at(i));  // get a centroid

                // for each neighbor centroid
                for (int j = i; j < this->centroids.size(); j++)    
                {
                    next_centroid = &(this->centroids.at(j));

                    if (i == j)     // dont compute the distance between same centers (it is 0)
                        continue;

                    // calculate the distance between 2 centroids
                    if (!this->represantation.compare("vector")) {
                        distance = euclidean_distance(curr_centroid, next_centroid);
                    }
                    else if (!this->represantation.compare("curve")) {
                        vector<T> mean_curve;
                        distance = Discrete_Frechet_distance(curr_centroid, next_centroid, mean_curve);
                    }

                    distances_radius.push_back(distance);
                }


                // get all vectors from L/probes buckets for each centroid
                if (!method.compare("LSH") || !method.compare("LSH_Frechet")) {
                    vectors_in_buckets = this->lsh->get_all_vectors_from_L_buckets(curr_centroid);
                }
                else if (!method.compare("Hypercube")) {
                    // each centroid has M vectors to compare
                    vectors_in_buckets = this->hypercube->get_M_vectors_from_buckets(curr_centroid); 
                }   

                vectors_in_buckets_for_all_centroids.push_back(vectors_in_buckets); 
            }
        }

        // get the (minimum distance / 2)
        min_radius = *min_element(distances_radius.begin(), distances_radius.end()) / 2; 
        distances.clear();

        // for each centroid run the range search with lsh/hypercube
        for (centroid_index = 0; centroid_index < this->centroids.size(); centroid_index++)
        {
            curr_centroid = &(this->centroids.at(centroid_index));

            // initialize each cluster 
            if (initialized_cluster == false) 
            {  
                curr_cluster = (node_cluster<T>*) malloc(sizeof(node_cluster<T>));
                curr_cluster->centroid = centroid_index;
                curr_cluster->vectors = new vector<struct dataNode<T>*>();
                this->clusters.push_back(curr_cluster); 
            }

            //  Double the radius
            curr_radius = r * min_radius;
            // curr_radius = pow(min_radius, 2);
            
            curr_cluster = this->clusters.at(centroid_index);   // get the current cluster
           
            // get the vectors from L buckets for the centroid (doesn't matter the radius)
            vectors_in_buckets = vectors_in_buckets_for_all_centroids.at(centroid_index);

            // Range search for the centroid
            // for each vector check if it is inside radius
            for (int i = 0; i < vectors_in_buckets->size(); i++)
            {
                curr_vector = &vectors_in_buckets->at(i);

                if (!this->represantation.compare("vector")) {
                    distance = euclidean_distance(curr_vector->vec, curr_centroid);
                }
                else if (!this->represantation.compare("curve")) {
                    vector<T> mean_curve;
                    distance = Discrete_Frechet_distance(curr_vector->vec, curr_centroid, mean_curve);
                }


                if (distance > curr_radius) { // if vector isn't inside radius pass
                    continue;
                } 

                item_id = *curr_vector->item_id;    // get the id from the vector

                //  assign the vector if we already didn't
                if ( (flags_vectors.find(item_id) == flags_vectors.end()) ) {
                    // insert the vector in current cluster
                    new_vector = dataset->get_vector(item_id);
                    curr_cluster->vectors->emplace_back(new_vector);
                    inserts_per_iterations++;
                    // points_which_assignned_with_lsh++;

                    // mark the vector it is inside a cluster
                    flags_vectors[item_id] = true;      
                    assigned_vectors[item_id] = centroid_index;  
                }
                else if ((flags_vectors.find(item_id) != flags_vectors.end()) && ((flags_vectors[item_id] == true) && (assigned_vectors.at(item_id) != centroid_index))) {
                    // if a vector belong in a cluster but an another cluster want this vector
                    // we calculate the distance between centroids and we take the nearest centroid 

                    // get the current centroid for the vector
                    curr_assigned_centroid_index = assigned_vectors.at(item_id);  
                    curr_assigned_centroid = &this->centroids.at(curr_assigned_centroid_index);
                    curr_assigned_cluster = this->clusters.at(curr_assigned_centroid_index);


                    /* find the nearest centroid */
                    // find the distance between vector and assigned centroid
                    if (!this->represantation.compare("vector")) {
                        distance = euclidean_distance(curr_vector->vec, curr_assigned_centroid);
                    }
                    else if (!this->represantation.compare("curve")) {
                        vector<T> mean_curve;
                        distance = Discrete_Frechet_distance(curr_vector->vec, curr_assigned_centroid, mean_curve);
                    }
                    distances.push_back(distance);

                    // get the distance between current centroid and vector
                    distance = distance;
                    distances.push_back(distance);

                    // get the index from minimum distance in the vector of distances
                    min_distance_index = min_element(distances.begin(), distances.end()) - distances.begin();
                    distances.clear();

                    if (min_distance_index == 1) { // if  min_distance_index == 1, the current centroid is nearest to the vector 

                        // assigned the vector to the new centroid
                        new_vector = dataset->get_vector(*curr_vector->item_id);

                        assigned_vectors[item_id] = centroid_index; // mark the vector with the new centroid
                        curr_cluster->vectors->push_back(new_vector);
                        inserts_per_iterations++;


                        // remove the vector from the assigned cluster
                        vector<struct dataNode<Type>*>::iterator it =  curr_assigned_cluster->vectors->begin();
                        while(it != curr_assigned_cluster->vectors->end())  
                        {
                            if (!((*it)->item_id->compare(item_id)))    
                            {   // we have a collision between 2 radius
                                curr_assigned_cluster->vectors->erase(it);
                                collision_centroids = true;
                            }
                            else 
                                ++it;
                        }
                    }
                }

            }

            // save the number of insertions for this cluster
            insertions_per_centroids.push_back(inserts_per_iterations);
            unassigned_vectors -= inserts_per_iterations;  
            inserts_per_iterations = 0;
        }


        // count the number of clusters that did not initialize any vector
        stable_clusters = count(insertions_per_centroids.begin(), insertions_per_centroids.end(), 0);
        insertions_per_centroids.clear();

        // cout << "Check for stable clusters\n";
        // if mostly centroids dont isnert vectors  {or we have at least one collision with radiuses, (optional)} start the brute force
        if ((stable_clusters >= this->number_of_clusters/2 ) /*|| ( collision_centroids == true)*/)   
        {   
            // if there are unassigned vectors, assign with brute force
            if (unassigned_vectors > 0) {  

                begin = std::chrono::steady_clock::now();   // keep the time for assign vectors with brute force

                // cout << "Start brute force\n";
                for (int i = 0; i < dataset->get_num_vectors(); i++)
                {
                    curr_node = dataset->get_i_vector(i);
                    item_id = *curr_node->item_id;
                    
                    
                    if ((flags_vectors.find(item_id) != flags_vectors.end()) & (flags_vectors[item_id] == true)) { // pass the assigned vectors
                        points_which_assignned_with_method++;
                        continue;
                    } 
  
                    points_which_assignned_with_brute_force++;

                    // get the unassigned vector
                    new_vector = curr_node;  


                    // calculate the distance between all centroids
                    for (int i = 0; i < this->centroids.size(); i++) {
                        if (!this->represantation.compare("vector")) {
                            distance = euclidean_distance(curr_node->vec, &this->centroids.at(i));
                        }
                        else if (!this->represantation.compare("curve")) {
                            vector<T> mean_curve;
                            distance = Discrete_Frechet_distance(curr_node->vec, &this->centroids.at(i), mean_curve);
                        }
                        distances.push_back(distance);
                    }

                    //find the index for the nearest centroid
                    centroid_index = min_element(distances.begin(), distances.end()) - distances.begin();
                    distances.clear();

                    // assigned the vector to the nearest centroid
                    curr_cluster = this->clusters.at(centroid_index);
                    curr_cluster->vectors->push_back(new_vector);

                }

                end = std::chrono::steady_clock::now();
                time_for_brute_force += chrono::duration_cast<chrono::microseconds>(end - begin).count()/1000000.0;          
            }
            
            // calculate the new centroids (means of vectors in each cluster / mean frechet, based on represantation)
            if (!this->update_method.compare("Mean Vector")) {
                new_centroids = this->update_centroids_Mean_Vectors();
            }
            else if (!this->update_method.compare("Mean Frechet")) {
                new_centroids = this->update_centroids_Mean_Frechet();
            }


            // check the distance between new centroids and the old centroids
            // if distance between old_centroids and new_centroids
            // is smaller than 0.1, mark the centroid, has the optimal position
            for (int i = 0; i < this->centroids.size(); i++)
            {
                float distances_bw_oldcentroids = 0.0;
                if (!this->represantation.compare("vector")) {
                    distances_bw_oldcentroids = euclidean_distance(this->centroids.at(i), new_centroids.at(i));
                }
                else if (!this->represantation.compare("curve")) {
                    vector<T> mean_curve;
                    distances_bw_oldcentroids = Discrete_Frechet_distance(&this->centroids.at(i), &new_centroids.at(i), mean_curve);
                }

                cout << "Distance between previous centroids: " << distances_bw_oldcentroids << endl;
                if (distances_bw_oldcentroids < 1.0)    
                {   // we have the optimal centroid
                    flags_centroids.push_back(true);
                }
                else {      
                    // we need to find another centroid
                    flags_centroids.push_back(false);
                }
            }

            // count the non-movable centroids
            for(auto b : flags_centroids)
            {
                if ( b == true)
                    count_flags++;
            }  

            cout << "unmovable centroids = " << count_flags << endl;
            if ((count_flags >= this->number_of_clusters/2) || (iterations <= 0)) // if the half or more centroids don't move or we have consume all iterations
            {  
                for (int i = 0; i < vectors_in_buckets_for_all_centroids.size(); i++) // free the memory
                {
                    delete vectors_in_buckets_for_all_centroids.at(i);
                }
                break;  // stop the proces, we have the optimal clusters
            }


            // clear all clusters and start the process from begin with the new centroids
            for (int i = 0; i < this->clusters.size(); i++) 
            {
                curr_cluster = this->clusters.at(i);
                delete curr_cluster->vectors;
                free(curr_cluster);
            }
            
            this->clusters.clear();
            centroids_changed = true;
            number_of_vectors = dataset->get_num_vectors();

            // insert the new centroidsin clusters
            this->centroids.swap(new_centroids);

            /* clear all structures and counters */
            new_centroids.clear();     
            count_flags = 0;
            all_iterations++;
            flags_centroids.clear();   

            for (int i = 0; i < vectors_in_buckets_for_all_centroids.size(); i++)
            {
                delete vectors_in_buckets_for_all_centroids.at(i);
            }
            vectors_in_buckets_for_all_centroids.clear(); 
            flags_vectors.clear();
            assigned_vectors.clear();
            unassigned_vectors = number_of_vectors;
            initialized_cluster = false;
            curr_radius = 0;
            min_radius_temp = 0;
            distances_radius.clear();
            collision_centroids = false;
            stable_clusters = 0;
            r = 1;
        }
        else {  // double the radius to each cluster

            centroids_changed = false;
            initialized_cluster = true;
            r = r * 2;
        }

        iterations--;
    }

    // cout << "Points which assigned with brute force: " << points_which_assignned_with_brute_force  << endl;
    // cout << "Points which assigned with " << method << " : "  << points_which_assignned_with_method  << endl;
    // cout << "Time spent on assigning points with brute force: " << time_for_brute_force << endl;
}





//s(i) = [ b(i) - a(i) ] / max { a(i), b(i)};
//For 1 ≤ i ≤ n, a ( i ) = average distance of i to objects in same cluster
//Let b ( i ) = average distance of i to objects in next best (neighbor) cluster, i.e. cluster of 2nd closest centroid.
template<typename T>
void Cluster<T>::Silhouette(std::ofstream& output_file)
{
    if (this->number_of_clusters <= 1 )
    {
        output_file << "No silhouette for 1 cluster\n";
        return;
    }
    
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end; 

    std::unordered_map<string, float> map_distances;    // avoid  calculate the distance between 2 vectors several times (0 - 1 == 1 - 0)
    string curr_vector_id = "", hash_key = "", reverse_key = "";

    struct node_cluster<T>* curr_cluster = NULL;
    struct node_cluster<T>* nearest_cluster = NULL;
    struct dataNode<T>* curr_vector = NULL;
    struct dataNode<T>* next_vector = NULL;
    float a_i = 0.0, b_i = 0.0, s_i = 0.0;
    vector<float> distances_bw_centroids;
    float distance = 0.0;
    int nearest_centroid_index;
    vector<float> silhouette;
    float s_average = 0.0;
    
    begin = std::chrono::steady_clock::now();

    // for every cluster
    for (int i = 0; i < this->clusters.size(); i++)
    {
        s_i = 0.0;

        // if we have only one element in cluster mark the cluster (no silhouette)
        if ((this->clusters.at(i)->vectors->size() <= 1 ))
        {
            s_i = -2;
            silhouette.push_back(s_i);
            continue;
        }

        curr_cluster = this->clusters.at(i); // take a cluster

        // find the nearest cluster from this cluster
        for (int c = 0; c < this->centroids.size(); c++)
        {

            if (!this->represantation.compare("vector")) {
                distance = euclidean_distance(&centroids.at(curr_cluster->centroid), &centroids.at(c));
            }
            else if (!this->represantation.compare("curve")) {
                vector<T> mean_curve;
                distance = Discrete_Frechet_distance(&centroids.at(curr_cluster->centroid), &centroids.at(c), mean_curve);
            }
            
            distances_bw_centroids.push_back(distance); // included the distance 0 between the centroid itself
        }

        // find the index from second smaller distance (reminder: the first smaller distance is the 0)
        nearest_centroid_index = index_of_second_smallest_number(distances_bw_centroids);

        // get the nearest cluster
        nearest_cluster = this->clusters.at(nearest_centroid_index);

        // take a vector and caclucate the average distance between all vectors in cluster (a(i))
        for (int j = 0; j < curr_cluster->vectors->size(); j++)
        {
            curr_vector = curr_cluster->vectors->at(j);  
            curr_vector_id = *curr_vector->item_id;       // create the hash key for unorder map


            // calculate the distance between current vector and all vectors in cluster
            for (int v = 0; v < curr_cluster->vectors->size(); v++)
            {
                next_vector = curr_cluster->vectors->at(v);  // get the next vector in cluster
                hash_key = curr_vector_id + " - ";
                hash_key += *next_vector->item_id;
                reverse_key = reverseString(hash_key);

                if (*curr_vector->item_id == *next_vector->item_id)
                {
                    distance = 0;
                }   
                else if (map_distances.find(hash_key) != map_distances.end())   // search in map, maybe we have already calculate the distance
                {
                    distance = map_distances.at(hash_key);
                }
                else if (map_distances.find(reverse_key) != map_distances.end())
                {
                    distance = map_distances.at(reverse_key);
                }
                else {
                    // if we dont find the distance between 2 vectors in map calculate it
                    if (!this->represantation.compare("vector")) {
                        distance = euclidean_distance(curr_vector->vec, next_vector->vec);  
                    }
                    else if (!this->represantation.compare("curve")) {
                        vector<T> mean_curve;
                        distance = Discrete_Frechet_distance(curr_vector->vec, next_vector->vec, mean_curve);
                    }

                    map_distances[reverse_key] = distance; 
                    map_distances[hash_key] = distance;     //(this is dumb)
                }

                a_i += distance;
                hash_key.clear();
            }

            hash_key.clear();
            reverse_key.clear();


            a_i /= (curr_cluster->vectors->size()-1);    // take the average distance 
            
            // take the average distance from points in nearest cluster
            for (int v = 0; v < nearest_cluster->vectors->size(); v++)
            {
                next_vector = nearest_cluster->vectors->at(v);
                hash_key  = curr_vector_id + " - ";
                hash_key += *next_vector->item_id;
                reverse_key = reverseString(hash_key);
                
                // search in map the distance between current vector and next vector 
                if (map_distances.find(hash_key) != map_distances.end())
                {
                    distance = map_distances.at(hash_key);
                }
                else if (map_distances.find(reverse_key) != map_distances.end())
                {
                    distance = map_distances.at(reverse_key);
                }
                else {  // if didnt foud the distance in map calculate and insert it in map 
                    
                    if (!this->represantation.compare("vector")) {
                        distance = euclidean_distance(curr_vector->vec, next_vector->vec);  
                    }
                    else if (!this->represantation.compare("curve")) {
                        vector<T> mean_curve;
                        distance = Discrete_Frechet_distance(curr_vector->vec, next_vector->vec, mean_curve);
                    }
                    map_distances[reverse_key] = distance;  // e.g distance between 0 - 1 and 1 - 0 is equal
                    map_distances[hash_key] = distance;
                }

                b_i += distance;
                hash_key.clear();
            }

            hash_key.clear();
            reverse_key.clear();

            if (nearest_cluster->vectors->size() == 0 ) {
                b_i /= 1;
            }
            else {
                b_i /= nearest_cluster->vectors->size();
            }

            s_i += (b_i - a_i) / max_number(a_i, b_i);      // get the si

            b_i = 0.0;
            a_i = 0.0;
        }

        // get the average silhouette
        if (curr_cluster->vectors->size() <= 1) {  // no silhouette for clusters with <=1 vectors 
            s_i = -2;
        }
        else {
            s_i /= curr_cluster->vectors->size(); 
        }

        silhouette.push_back(s_i);
        distances_bw_centroids.clear();
    }

    end = std::chrono::steady_clock::now();


    //print the silhouette
    output_file << "Silhouette: [";
    for (auto s : silhouette)
    {
        if (s == -2 )
            output_file << "No silhouette, ";   
        else {
            output_file << s << ", ";
            s_average += s;
        }
    }

    output_file << s_average / silhouette.size() << "]\n";
    cout << "silhouette time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000000.0 << "[seconds]" << std::endl;
}








template <typename T>
void Cluster<T>::print(ofstream& output_file, bool complete)
{
    node_cluster<T>* curr_node;
    struct dataNode<T>* curr_dataNode;
    int index;

    for (int i = 0; i < this->clusters.size(); i++)
    {
        curr_node = this->clusters.at(i);
        index = curr_node->centroid;

        output_file << "Cluster - " <<  index+1 << ": { ";
        output_file << "size: " <<  curr_node->vectors->size() << ", (";
        for (auto c : centroids.at(i))
            output_file << c << ", ";
        output_file << ") ";

        if (complete) {  // print the items in cluster
            output_file << ",  ";
            for (int j = 0; j < curr_node->vectors->size(); j++)
            {
                //dont print the centroid as a vector in cluster-i  
                // if (*curr_node->vectors->at(j)->item_id != *this->centroids.at(index)->item_id) 
                output_file << *curr_node->vectors->at(j)->item_id << ", ";
            }

        }


        output_file << "}" << endl;
    }




}