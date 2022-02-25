#include <iostream>
#include "Cluster.h"
#include "Cluster.cpp"
#include <unistd.h>

using namespace std;

int main(int argc, char *argv[])
{
    if ((argc < 9) || (argc > 14)) {
        cout << "Failure usage: ./cluster -i <input file> -c <configuration file> -o <output file> -update <Mean Frechet or Mean Vector> -assignment <Classic or LSH or Hypercube or LSH_Frechet> -complete <optional> -silhouette <optional>" << endl;
        exit(EXIT_FAILURE);
    }

    string metric = "";     
    string represantation = "";
    string input_file_str = "";
    string query_file_str = "";
    string output_file_str;
    string conf_file_str = "";
    string update_method = "";         //can be Mean Frechet or Mean Vector
    string assignment = "";     // can be Classic or LSH or Hypercube or LSH_Frechet
    bool complete = false;
    bool shilhouette = false;


    // get the arguments from command line
    for (int i = 0; i < argc; i++)
    {
        if (!strcmp(argv[i], "-i")) {
            input_file_str = argv[++i];
        }
        else if (!strcmp(argv[i], "-q")) {
            query_file_str = argv[++i];
        }
        else if (!strcmp(argv[i], "-o")) {
            output_file_str = argv[++i];
        }
        else if (!strcmp(argv[i], "-c")) {   
           conf_file_str = argv[++i];
        }
        else if (!strcmp(argv[i], "-update")) {
            update_method = argv[++i];
            update_method += " ";
            update_method += argv[++i]; 

            if (!update_method.compare("Mean Vector")) {
                represantation = "vector";
            }
            else if (!update_method.compare("Mean Frechet")) {
                represantation = "curve";
            }
        }
        else if (!strcmp(argv[i], "-assignment")) {   
            assignment = argv[++i];
        }
        else if (!strcmp(argv[i], "-complete")) {   //optional
            complete = true;
        }
        else if (!strcmp(argv[i], "-silhouette")) {
            shilhouette = true;
        }
    }

    if ((update_method == "Mean Frechet") && (represantation != "curve")) {
        cout << "Can't run clustering with update: " << update_method << " and assignment " << assignment << endl; 
        exit(EXIT_FAILURE);
    }


    // files
    ifstream input_file(input_file_str);
    ifstream conf_file(conf_file_str);
    ofstream output_file(output_file_str);


    // check if input/query file exist, if it doesn't request the correct files
    while (!exists_file(input_file_str)) {
        cout << "Please enter a valid input file name: " << flush;
        input_file_str.clear();
        getline(cin, input_file_str);
        input_file.open(input_file_str.c_str());
        if (input_file) break;
    }

    // check if file.conf file exist, if it doesn't, request the correct files
    while (!exists_file(conf_file_str)) {
        cout << "Please enter a valid file_name.conf: " << flush;
        conf_file_str.clear();
        getline(cin, conf_file_str);
        conf_file.open(conf_file_str.c_str());
        if (conf_file) break;
    }


    std::chrono::steady_clock::time_point begin_Lloyds, begin_all;
    std::chrono::steady_clock::time_point end_Lloyds, end_all; 

    begin_all = std::chrono::steady_clock::now();  

    cout << "\n\nUpdate method: " << update_method << endl;
    cout << "Assignment: " << assignment << endl;


    Dataset<double>* dataset = new Dataset<double>(represantation, metric); // create the dataset
    string line;
    int D = 0;
    int hypercube = 0, lsh = 1;
    if (input_file.is_open())        // insert the vectors in dataset
    {
        while(getline(input_file, line)) {
            dataset->insert(line);
        }
    }
    else 
        cout << "Uable to open Data file";  

    Cluster<double>* cluster = new Cluster<double>(dataset, assignment, represantation, conf_file_str, update_method);        // create the cluster
    output_file << "Algorithm: " << assignment << ", " << update_method << endl;
    
    begin_Lloyds = std::chrono::steady_clock::now();    
    if (!assignment.compare("Classic")) {
        cluster->Lloyds();           // run algorithm Lloyd's 
    }  
    else if (!assignment.compare("LSH") || !assignment.compare("Hypercube") || !assignment.compare("LSH_Frechet")) {
        cluster->Range_Search(assignment);                            //run algorithm range search lsh or hypercube
    }


    end_Lloyds = std::chrono::steady_clock::now();

    cluster->print(output_file, complete);   // print all {size,centroid}
    

    output_file << "\nclustering time: " << chrono::duration_cast<std::chrono::microseconds>(end_Lloyds - begin_Lloyds).count()/1000000.0 << " seconds" << endl; // print the run time for clustering

    if (shilhouette == true) {
        cluster->Silhouette(output_file);  // run the silhouette and print the results
    }

    output_file.close();
    
    delete cluster;
    delete dataset;

    end_all = std::chrono::steady_clock::now();
    cout << "\nRun time: " << chrono::duration_cast<chrono::microseconds>(end_all - begin_all).count()/1000000.0 << " seconds" << endl;
}   