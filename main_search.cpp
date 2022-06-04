#include <iostream>
#include "Utilities.h"
#include "Dataset.h"
#include "Dataset.cpp"
#include "LSH.h"
#include "LSH.cpp"
#include "Hypercube.h"
#include "Hypercube.cpp"

#define input_file_path "/home/nick/Desktop/project_Emiris/Datasets/nasdaq2017_LQ.csv"
#define query_file_path ""

int main(int argc, char* argv[])
{
    if ((argc < 9) || (argc > 20) ) {
        cout << "Failure usage: ./search -i <input file> -q <query file> -k <int> -L <int> -M <int> -probes <int> -o <output file> -algorithm <LSH or Hypercube or Frechet> -metric <discrete or continuous | only for -algorithm Frechet> -delta <double>" << endl;
        exit(EXIT_FAILURE);
    }

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    // for the time
    std::chrono::_V2::system_clock::time_point tAlgorithm_start, tAlgorithm_end, tTrue_start, tTrue_end;    
    std::chrono::steady_clock::time_point begin, end, begin_filling, end_filling;
    double time_for_filling_structs = 0.0;
    double sum_tTrue = 0.0;
    double sum_tAlgorithm = 0.0;

    string line;

    string input_file_str = "";
    string query_file_str = "";
    string output_file_str = "output_file.txt";
    int k = 0, probes = 0, m = 0, l = 0;
    int n = 0, D, w;
    double delta = 0.0;
    string algorithm_str = "", metric = "";
    string represantation = "";
    int not_optionals = 0;
    int check_delta = 0;
    int number_of_queries = 0;
    int number_of_nearest = 1;
    vector<double> approximation_factors;
    double approximation_factor, Max_approximation_factor;
    double sum_of_lshDistances = 0.0, sum_of_trueDistances = 0.0;

    
    

    // get the arguments from command line
    for (int i = 0; i < argc; i++)
    {
        if (!strcmp(argv[i], "-i")) {
            input_file_str = argv[++i];
            not_optionals++;
        }
        else if (!strcmp(argv[i], "-q")) {
            query_file_str = argv[++i];
            not_optionals++;
        }
        else if (!strcmp(argv[i], "-o")) {
            output_file_str = argv[++i];
            not_optionals++;
        }
        else if (!strcmp(argv[i], "-k")) {   //optional
            k = stoi(argv[++i]);
        }
        else if (!strcmp(argv[i], "-M")) {   //optional
            m = stoi(argv[++i]);
        }
        else if (!strcmp(argv[i], "-probes")) { // optional
            probes = stoi(argv[++i]);
        }
        else if (!strcmp(argv[i], "-L")) {   //optional
            l = stoi(argv[++i]);
        }
        else if (!strcmp(argv[i], "-algorithm")) {   
            algorithm_str = argv[++i];
            if (!algorithm_str.compare("LSH") || !algorithm_str.compare("Hypercube")) {
                represantation = "vector";
            }
            else if (!algorithm_str.compare("Frechet")) {
                represantation = "curve";
                check_delta++;
            }
            not_optionals++;
        }
        else if (!strcmp(argv[i], "-metric")) {   //only for Frechet
            metric = argv[++i];
            if (!algorithm_str.compare("LSH") || !algorithm_str.compare("Hypercube")) {
                cout << "Failure usage: ./search -i <input file> -q <query file> -k <int> -L <int> -M <int> -probes <int> -o <output file> -algorithm <LSH or Hypercube or Frechet> -metric <discrete or continuous | only for -algorithm Frechet> -delta <double>" << endl;
                exit(EXIT_FAILURE);
            }
            if (!metric.compare("discrete")) {
                algorithm_str = "LSH_Frechet_Discrete";
                check_delta++;
            }
            else if (!metric.compare("continuous")) {
                algorithm_str = "LSH_Frechet_Continuous";
                check_delta++;
            }
        }
        else if (!strcmp(argv[i], "-delta")) {   //optional
            delta = atof(argv[++i]);
            check_delta++;
        }
    }


    if ((not_optionals != 4) || ((check_delta != 0) && (check_delta < 3))) {
        cout << "Failure usage: ./search -i <input file> -q <query file> -k <int> -L <int> -M <int> -probes <int> -o <output file> -algorithm <LSH or Hypercube or Frechet> -metric <discrete or continuous | only for -algorithm Frechet> -delta <double>" << endl;
        exit(EXIT_FAILURE);
    }


    Dataset<double>* dataset = new Dataset<double>(represantation, metric); // create the dataset
    ifstream input_file(input_file_str);
    ifstream query_file(query_file_str);
    ofstream output_file;
    output_file.open(output_file_str, std::ofstream::out | std::ofstream::trunc);  
    LSH<double>* lsh;
    Hypercube<double>* hypercube;
    vector<vector_distances*> DistancesTrue;    // an array wich will keep the true nearest neighbors
    vector<vector_distances*> Distances;
    _vector* vec_query = NULL;
    dataNode<double>* current_vector;
    string item_id;
    static double M = 0.0;      // for padding




    // check if input/query file exist, if it doesn't request the correct files
    while (!exists_file(input_file_str)) {
        cout << "Please enter a valid input file name: " << flush;
        input_file_str.clear();
        getline(cin, input_file_str);
        input_file.open(input_file_str.c_str());
        if (input_file) break;
    }

    while (!exists_file(query_file_str)) {
        cout << "Please enter a valid query file name: " << flush;
        query_file.clear();
        getline(cin, query_file_str);
        query_file.open(query_file_str.c_str());
        if (query_file) break;
    }


    begin = std::chrono::steady_clock::now();       // begin program after got the right input files
    begin_filling = std::chrono::steady_clock::now();   // begin filling the dataset


    // insert vectors in dataset
    if (input_file.is_open())        
    {
        while(getline(input_file, line)) {
            dataset->insert(line);
        }
    }
    else {
        cout << "Uable to open Data file";
        exit(EXIT_FAILURE);
    }

    // dataset->print(1);      // print the n first curves


    n = dataset->get_num_vectors();
    cout << "\n\nNumber of vectors in dataset: " << n << endl;

    D = dataset->get_dimensions();
    w = ceil(dataset->calculate_W());     // calculate w based on a data subset

    if (!algorithm_str.compare("LSH"))  // if algorithm is LSH initialize structures for lsh
    {
        cout << "W: " << w/3 << endl;

        // compute the number of buckets for hash tables based of number of vectors in dataset
        n = ceil(n / 16.0);

        lsh = new LSH<double>(n, k, l, w/3, D, delta, represantation, metric);

        cout << "Number of buckets for each hash table: " << n << endl;
        cout << "Number of hash functions h: " << lsh->get_num_of_h() << endl;
        cout << "Number of hash tables: " << lsh->get_num_hashTables() << endl;

        // insert vectors in Hash Tables (LSH)
        for (int i = 0; i < dataset->get_num_vectors(); i++) {
            current_vector = dataset->get_i_vector(i);     //get a pointers from dataset
            lsh->insert_vector(current_vector);        //insert pointers to lsh structs
        }

    }
    else if (!algorithm_str.compare("Hypercube"))    // if algorithm is Hypercube initialize the structure for hypercube
    {
        cout << "W: " << w/3 << endl;

        hypercube = new Hypercube<double>(k, w/3, m, probes, D);
        cout << "Number of hash functions h: " << hypercube->get_num_of_h() << endl;
        cout << "Number of M vectors for check: " << hypercube->max_num_of_vectors() << endl;
        cout << "Number of vertices for check: " << hypercube->max_num_of_vertices() << endl;

        // insert vectors in Hash Table (Hypercube)
        for (int i = 0; i < dataset->get_num_vectors(); i++) {
            current_vector = dataset->get_i_vector(i);           //get a pointers from dataset
            int k = hypercube->insert_vector(current_vector);    //insert pointers to lsh structs

        }
    }
    else if (!algorithm_str.compare("LSH_Frechet_Discrete")) 
    {
        cout << "W: " << w/2 << endl;

        n = ceil(n / 16.0);     // compute the number of buckets for hash tables based of number of vectors in dataset
        
        lsh = new LSH<double>(n, k, l, w/2, D, delta, represantation, metric);

        cout << "Number of buckets for each hash table: " << n << endl;
        cout << "Number of hash functions h: " << lsh->get_num_of_h() << endl;
        cout << "Number of hash tables: " << lsh->get_num_hashTables() << endl;
        cout << "delta for Grid: " << delta << endl;
        cout << "Value for padding: " << value_for_padding(0.0) + 1000 << " (maximum value in dataset + 1000)" << endl;

        vector<int> buckets; 
        vector<int> all_buckets;
        // insert vectors in Hash Tables (LSH)
        for (int i = 0; i < dataset->get_num_vectors(); i++) {
            current_vector = dataset->get_i_vector(i);     //get a pointers from dataset
            lsh->insert_vector(current_vector);        //insert pointers to lsh structs
        }
    }
    else if (!algorithm_str.compare("LSH_Frechet_Continuous")) 
    {
        cout << "W: " << w/2 << endl;
        w = 50;
        n = ceil(n / 16.0);     // compute the number of buckets for hash tables based of number of vectors in dataset
        l = 1;          
        
        lsh = new LSH<double>(n, k, l, w/2, D, delta, represantation, metric);

        cout << "Number of buckets for each hash table: " << n << endl;
        cout << "Number of hash functions h: " << lsh->get_num_of_h() << endl;
        cout << "Number of hash tables: " << lsh->get_num_hashTables() << endl;
        cout << "delta for Grid: " << delta << endl;
        cout << "Value for padding: " << value_for_padding(0.0) + 1000 << " - maximum value in dataset + 1000 -" << endl;

        vector<int> all_buckets;
        vector<int> buckets;


        // insert vectors in Hash Tables (LSH)
        for (int i = 0; i < dataset->get_num_vectors(); i++) {
            current_vector = dataset->get_i_vector(i);     //get a pointers from dataset
            buckets = lsh->insert_vector(current_vector);        //insert pointers to lsh structs
            for (auto b : buckets)
                all_buckets.push_back(b);
        }

        all_buckets.erase(std::unique(all_buckets.begin(), all_buckets.end(), [](int a, int b)
        {
            return (a == b);
        }
        ), all_buckets.end());

    }

    // lsh->print();

    end_filling = std::chrono::steady_clock::now();
    time_for_filling_structs = chrono::duration_cast<chrono::microseconds>(end_filling - begin_filling).count()/1000000.0;          
    cout << "\nTime for filling structures: " << time_for_filling_structs << " seconds" << endl;
    int num_query = 1;
    string answer = "";

    while(1) {
        approximation_factors.clear();
        sum_tAlgorithm = 0.0;
        sum_tTrue = 0.0;
        number_of_queries = 0;
        num_query = 1;

        if (query_file.is_open()) {


            while(getline(query_file, line)) {

                number_of_queries++;

                item_id =  line.substr(0, line.find_first_of(" \t")).c_str();
                line = line.substr(line.find_first_of(" \t")+1);    // remove the item_id from vector
                vec_query = convert_string_to_vector(line); 

                if (!represantation.compare("curve") && !metric.compare("continuous")) {    
                    // filtering the queries
                    filtering(vec_query->vec);
                }


                output_file << "Query: " << item_id << endl;
                output_file << "Algorithm: " << algorithm_str << endl;  //{LSH_Vector, Hypercube, LSH_Frechet_Continuous, LSH_Frechet_Discrete}

                // run algorithm for nearest (LSH/ Hypercube/ Frechet with LSH) and keep the run time
                tAlgorithm_start = high_resolution_clock::now();
                if (!algorithm_str.compare("LSH") || !algorithm_str.compare("LSH_Frechet_Discrete") || !algorithm_str.compare("LSH_Frechet_Continuous")) {

                    Distances = lsh->findNN_LSH(vec_query->vec, number_of_nearest); 
                    if (Distances.size() == 0)
                        output_file << "Approximate Nearest neighbor: Not found" << endl;

                }
                else if (!algorithm_str.compare("Hypercube")) {
                    Distances = hypercube->findNN_Hypercube(vec_query->vec, number_of_nearest);
                    if (Distances.size() == 0)
                        output_file << "Approximate Nearest neighbor: Not found" << endl;
                }
            
                tAlgorithm_end = high_resolution_clock::now();   

                // run brute force and keep the run time
                auto tTRUE_start = high_resolution_clock::now();
                DistancesTrue = dataset->findNN_true(vec_query->vec, number_of_nearest);
                auto tTRUE_end = high_resolution_clock::now();    


                for (int i = 0; i < number_of_nearest; i++)
                {
                    if (Distances.size() > 0) {
                        output_file << "Approximate Nearest neighbor: " << *Distances.at(i)->item_id << endl;
                    }
                    output_file << "True Nearest neighbor: " <<  *DistancesTrue.at(i)->item_id << endl;

                    if (Distances.size() > 0) {
                        output_file << "distanceApproximate: " <<  std::fixed 
                            << std::setprecision(3) << Distances.at(i)->distance_with_query << endl;
                    } 

                    output_file << "distanceTrue: " << std::fixed
                        << std::setprecision(3) << DistancesTrue.at(i)->distance_with_query << endl;
                    
                    if (Distances.size() > 0) {
                        approximation_factor = Distances.at(0)->distance_with_query / DistancesTrue.at(0)->distance_with_query; // true distance / approximate distance 
                        approximation_factors.push_back(approximation_factor);
                    }
                }

                if (Distances.size() > 0) { // for accuracy
                    sum_of_lshDistances += Distances.at(0)->distance_with_query;
                    sum_of_trueDistances += DistancesTrue.at(0)->distance_with_query;
                }


                // print run time for LSH and brute force
                duration<double, std::milli> tTRUE = tTRUE_end - tTRUE_start;
                duration<double, std::milli> tAlgorithm = tAlgorithm_end - tAlgorithm_start;
                sum_tTrue += tTRUE.count();
                sum_tAlgorithm += tAlgorithm.count();

                // output_file << "tApproximate: " << std::fixed << std::setprecision(3) << tAlgorithm.count() << endl;
                // output_file << "tTRUE: " << std::fixed << std::setprecision(3) << tTRUE.count() << endl;

                // free the keepers for this query go to the next query
                for (int i = 0; i < DistancesTrue.size(); i++) {
                    free(DistancesTrue.at(i));
                }
                for (int i = 0; i < Distances.size(); i++) {
                    free(Distances.at(i));
                }

                output_file << "\n\n\n";

                delete vec_query->vec;
                free(vec_query);

            }

        }
        else {
            cout << "Uable to open query file";
            exit(EXIT_FAILURE);
        }

        output_file << "tApproximateAverage: " << std::fixed << std::setprecision(3) << sum_tAlgorithm / number_of_queries << endl;
        output_file << "tTrueAverage: " << std::fixed << std::setprecision(3) << sum_tTrue / number_of_queries << endl;
        if (approximation_factors.size() > 0) {
            output_file << "MAF: " << *max_element(approximation_factors.begin(), approximation_factors.end()) << endl;
        }
        else  {
            output_file << "MAF: " << "Not Found\n";
        }
        if (sum_of_lshDistances > 0) {
            cout << "accuracy: " << std::fixed
                << std::setprecision(2)<< (sum_of_trueDistances / sum_of_lshDistances)*100 << " %" << endl;
        }
        // cout << algorithm_str << " is " << sum_tTrue / sum_tAlgorithm << " ms faster than brute force"<< endl;



        cout << "Would you like to continue with another query file (y or n): ";
        fflush(stdout);
        getline(cin, answer);
        fflush(stdin);

        if (!answer.compare("y")) {
            query_file.close();
            query_file_str = "";
            query_file.clear();
            while (!exists_file(query_file_str)) {
                cout << "Please enter a valid query file name: " << flush;
                query_file.clear();
                getline(cin, query_file_str);
                query_file.open(query_file_str.c_str());
                if (query_file) break;
            }

            cout << "Enter a output file name: " << flush;

            output_file_str.clear();
            output_file.close();

            getline(cin, output_file_str);

            output_file.open(output_file_str, std::ofstream::out | std::ofstream::trunc);  

        }
        else 
            break;

    }


    delete dataset;
    if (!algorithm_str.compare("LSH") || !algorithm_str.compare("LSH_Frechet_Discrete") || !algorithm_str.compare("LSH_Frechet_Continuous") ) {
        delete lsh;
    }
    else if (!algorithm_str.compare("Hypercube")) {
        delete hypercube;
    }

    return 0;
}