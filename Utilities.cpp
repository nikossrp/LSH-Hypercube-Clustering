#include <iostream>
#include "Utilities.h"

using namespace std;

vector<vector<float>*> vectors_v(int k, int D) 
{
    vector<vector<float>*> vector_v;
    vector<float>* vec;

    // create the random vectors v with normal distribution v ~ N(0, 1)
    unsigned seed = chrono::steady_clock::now().time_since_epoch().count();
    default_random_engine e(seed);
    normal_distribution<double> distrN{0.0, 1.0};
    for (int i = 0; i < k; i++) 
    {
        vec = new vector<float>;

        for (int i = 0; i < D; i++) {
            vec->push_back(distrN(e));
        }
        vector_v.push_back(vec);
    }

    return vector_v;
}

vector<float> vector_t(int k) 
{
    vector<float> t;
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0, W-1);
    for (int n = 0; n < k; n++) {
        // Use dis to transform the random unsigned int generated by gen into a 
        // double in [0, W). Each call to dis(gen) generates a new random double
        t.push_back(dis(gen));
    }

    return t;
}





double euclidean_distance(vector<Type>* vec1, vector<Type>* vec2) {
    Type sum_all = 0.0;
    int vectorsize1 = vec1->size();
    int vectorsize2 = vec2->size();

    // if (vectorsize1 != vectorsize2) {  
    //     cout << "Error: Vectors haven't the same size\n";
    //     cout << "Vector size 1: " << vectorsize1 << endl;
    //     cout << "Vector size 2: " << vectorsize2 << endl;

    //     exit(EXIT_FAILURE);
    // }

    int size = vectorsize1;
    float Nsum = 0.0;

    vector<Type>::iterator vec1_it = vec1->begin();
    vector<Type>::iterator vec2_it = vec2->begin();

    for (int i = 0; i < size; i++) { //(remainder: first item in vector is id)
        Nsum = vec1_it[i] - vec2_it[i];
        sum_all += pow(Nsum, 2.0);
    }
    return sqrt(sum_all); 
}


double euclidean_distance(vector<Type>& vec1, vector<Type>& vec2) {
    Type sum_all = 0.0;
    int vectorsize1 = vec1.size();
    int vectorsize2 = vec2.size();

    // if (vectorsize1 != vectorsize2) {  
    //     cout << "Error: Vectors haven't the same size\n";
    //     cout << "Vector1 size: " << vectorsize1 << endl;
    //     cout << "Vector2 size: " << vectorsize2 << endl;

    //     exit(EXIT_FAILURE);
    // }

    int size = vectorsize1;
    float Nsum = 0.0;

    vector<Type>::iterator vec1_it = vec1.begin();
    vector<Type>::iterator vec2_it = vec2.begin();

    for (int i = 0; i < size; i++) { //(remainder: first item in vector is id)
        Nsum = vec1_it[i] - vec2_it[i];
        sum_all += pow(Nsum, 2.0);
    }
    if (sum_all > 0)
        return sqrt(sum_all);
    else if (sum_all == 0)
        return 0; 
        
    exit(EXIT_FAILURE);
}



long int mod(long int x, long int y)
{
    // long long int res = x - ((floor( (long double)x / (long double)y) ) * y);
    long int res = ((x % y) + y) % y;
    return res;
}


_vector* convert_string_to_vector(string vec_string)
{
    double curr_coordinate = 0.0;
    _vector* vec = (_vector*)malloc(sizeof(_vector));
    int t = 1;
    std::vector<Type>* result = new vector<Type>(); 
    std::istringstream iss(vec_string); 
    int d = 0;


    for(std::string vec_string; iss >> vec_string; )  {
        curr_coordinate = atof(vec_string.c_str());
        result->push_back (curr_coordinate); 
        value_for_padding(curr_coordinate); // keep the greater value for padding
        d++;
    }

    vec->d = d;             //keep and check the dimensions (avoid trash files, all dimensions should be the same)
    vec->vec = result;
    return vec;
}


std::string toBinary(int n)
{
    std::string r;
    while(n!=0) {r=(n%2==0 ?"0":"1")+r; n/=2;}
    return r;
}


size_t count (const string & src, const string & str) {
    size_t cnt = 0, fnd = 0;
    while ((fnd = ((src.find(str, fnd)) != string::npos))) {
        cnt++; fnd++;
    }
    return cnt;
}


void get_vertices(char* str, int i, int changesLeft, std::vector<std::string>& all_nearest_vertices) {
        if (changesLeft == 0) {
            std::string str_str(str);
            all_nearest_vertices.push_back(str_str);
            return;
        }
        if (i < 0) return;
        // flip current bit
        str[i] = str[i] == '0' ? '1' : '0';
        get_vertices(str, i-1, changesLeft-1, all_nearest_vertices);
        // or don't flip it (flip it again to undo)
        str[i] = str[i] == '0' ? '1' : '0';
        get_vertices(str, i-1, changesLeft, all_nearest_vertices);
}

int binary(vector<float> arr, float to_search,int arr_size)
{ 
    int i=arr_size/2,index = -1; 
    bool next_small_bool = false, next_large_bool = false; 
    while((i>=0)&&(i<arr_size)){ 
     
        if(arr[i]==to_search){ 
     
            return i; 
        } 
        else{ 
            if(arr[i]>to_search){ 
                if(next_small_bool)break; 
                //dif = arr[i] - to_search; 
                if(index!=-1){ 
                    if(arr[index]>arr[i]){ 
                     index = i; 
                    } 
                } 
                else{ 
                    index = i; 
                } 
                i--; 
                next_large_bool = true; 
                 
            } 
     
            else if(arr[i]<to_search){ 
                if(next_large_bool)break; 
                //dif = arr[i] - to_search; 
                if(index!=-1) 
                { 
                        if(arr[index]<arr[i]){ 
                            index = i; 
                        } 
                } 
                else 
                { 
                    index = i; 
                } 
                i++; 
                next_small_bool = true; 
                } 
        } 
        } 
    return index; 
} 

int index_of_second_smallest_number(vector<float> v)
{
    //initialize values
    float smallest = *max_element(v.begin(), v.end());
    float second = smallest;
    
    for (int i = 0; i < v.size(); i++)
    {
        if (v.at(i) <= smallest)
        {
            second = smallest;
            smallest = v.at(i);
        }
        else if (v.at(i) < second)
        {
            second = v.at(i);
        }
    }
    
    //return the index of the second smaler number
    return find(v.begin(), v.end(), second) - v.begin();
    
}

float max_number(float x, float y)
{
    return (x < y) ? y : x;
}


// Function to reverse the given string
string reverseString(string str)
{
 
    // Reverse str using inbuilt function
    reverse(str.begin(), str.end());
 
    // Add space at the end so that the
    // last word is also reversed
    str.insert(str.end(), ' ');
 
    int n = str.length();
 
    int j = 0;
 
    // Find spaces and reverse all words
    // before that
    for (int i = 0; i < n; i++) {
 
        // If a space is encountered
        if (str[i] == ' ') {
            reverse(str.begin() + j,
                    str.begin() + i);
 
            // Update the starting index
            // for next word to reverse
            j = i + 1;
        }
    }
 
    // Remove spaces from the end of the
    // word that we appended
    str.pop_back();
 
    // Return the reversed string
    return str;
}


bool exists_file (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}


double max(double x, double y)
{
    return (x > y ) ? x : y;
}


double Discrete_Frechet_distance(std::vector<Type>* curve1, std::vector<Type>* curve2, std::vector<Type>& mean_curve)
{
    // if (curve1->size() != curve2->size()) {
    //     cout << "Error: Trying to compute discrete frechet distance with differend complexity\n";
    //     exit(EXIT_FAILURE);
    // }

    int sampling_frequency = 1;     // 1 day
    int time = 0;


    vector<pair<double, double>> curve1_R2;
    vector<pair<double, double>> curve2_R2;

    // Convert the curves from R^1 to R^2
    for (int i = 0; i < curve1->size(); i++) {
        time += sampling_frequency;      // get the x coordinate for each step

        curve1_R2.push_back(make_pair(time, curve1->at(i)));
        curve2_R2.push_back(make_pair(time, curve2->at(i)));
    }

    double distance_curves = 0.0;
    double distance_points;
    double distance_time;
    double x1, x2;
    double y1, y2;
    double curr_cord_curve1_x, curr_cord_curve2_x;
    int m1 = curve1->size();
    int m2 = curve2->size();

    // create the 2D array
    double** array = new double*[m1];
    for(int i = 0; i < m1; ++i)
        array[i] = new double[m2];


    // filling the array (Dynamic programming)
    // the distance between 2 curves in the axis x is 1
    for (int i = 0; i < m1; i ++)
    {
        for (int j = 0; j < m2; j++)    
        {
            // for curve 1
            x1 = curve1_R2[i].first;
            y1 = curve1_R2[i].second;

            // for curve 2
            x2 = curve2_R2[i].first;
            y2 = curve2_R2[i].second;

            // norm ||p1 - p2|| in R^2 
            distance_points = sqrt( pow(x1-x2, 2) + pow(y1-y2, 2) );   // sqrt ( (x1-x2)^2  +  (y1-y2)^2 )


            if ((i == 0) && (j == 0)) {   // initialize array    || p1 - p2 ||
                
                array[0][0] = distance_points;            
            
            }
            else if ((i == 0) && (j > 0)) {   // If i = 1; j > 1, then c(1, j) = max{ c(1, j-1), ||p1 - qj|| };

                array[0][j] = std::max ({ array[0][j-1], distance_points });
            
            }
            else if ((i > 0) && (j == 0)) {   // If i > 1; j = 1, then c(i, 1) = max{ c(i-1, 1), ||pi - q1|| };
                
                array[i][0] = std::max ({ array[i-1][0], distance_points });
            
            }
            else {  // if i > 1, j > 1 then : c(i, j) = max { min{ c(i-1, j), c(i-1, j-1), c(i, j-1) }, |pi - pj| }.
                
                array[i][j] = std::max ({ std::min({ array[i-1][j], array[i-1][j-1], array[i][j-1] }), distance_points });
            
            }
        }   
    }

    // the distance will be at the end of the array
    // with back tracking we can find the mean curve (?)
    distance_curves = array[m1-1][m2-1];   

    
    // if mean_curve isn't empty we have to calculate it, 
    // otherwise dont compute it and save time (we need it at  question B of the project_2)
    if (!mean_curve.empty()) {
        int i = m1-1, j = m2-1;
        int index_of_min;

        // find the optimal traversal
        vector<std::pair<int, int> > optimal_traversal;
        optimal_traversal.push_back(make_pair(i, j));

        // optimal traversal starting from c[m1][m2] -- and we go to --> c[0][0], we getting the minimum value each time, 
        // after that we will need to reverse this to have the right curve (see page. 28 from curves.pdf)
        while(1)   
        {
            if ((j == 0) && (i == 0))
                break;

            vector<double> steps;

            if (i > 0) {
                steps.push_back(array[i-1][j]);
            }
            if (j > 0) {
                steps.push_back(array[i][j-1]);
            }
            if ((i > 0) && (j > 0)) {
                steps.push_back(array[i-1][j-1]);
            }

            if (steps.size() == 3) {
                index_of_min = min_element(steps.begin(), steps.end()) - steps.begin();
            }
            else if ((i == 0) && (j > 0)) {
                index_of_min = 1;
            }
            else if ((i > 0) && j == 0) {
                index_of_min = 0;
            }

            if (index_of_min == 0) {
                i--;
                optimal_traversal.push_back(make_pair(i, j));
            }
            else if (index_of_min == 1) {
                j--;
                optimal_traversal.push_back(make_pair(i, j));
            }
            else {      // index_of_min == 2
                i--;
                j--;
                optimal_traversal.push_back(make_pair(i, j));
            }

            steps.clear();
        }
        
        // reverse the optimal traversal
        reverse(optimal_traversal.begin(), optimal_traversal.end());

        int max_size = curve1->size();

        int mean_curve_size = optimal_traversal.size();


        mean_curve.clear();

        pair<int, int> curr_pair;

        int num_unwanted_points = floor((float)(mean_curve_size - max_size) / 2.0);
        int start = num_unwanted_points;
        int end =   max_size + num_unwanted_points;

        // find the mean curve (mean curve can grow up too much, so we need a max size)
        // we cut the half unwanted point from the start of the curve
        // and the other half from the end of the curve
        for (int i = start; ; i++) {
            if (max_size <= 0) {
                break;
            }
            curr_pair = optimal_traversal.at(i);

            int index_curve1 = curr_pair.first;
            int index_curve2 = curr_pair.second;

            double mean_point;

            mean_point = (curve1->at(index_curve1) + curve2->at(index_curve2)) / 2;

            mean_curve.push_back(mean_point);

            max_size--;
        }

    }
    

    // free the 2D array
    for(int i = 0; i < m1; ++i) {
        delete [] array[i];
    }
    delete [] array;

    return distance_curves;
}



double Continuous_Frechet_distance(vector<Type>* my_curve1, vector<Type>* my_curve2)
{
    double distance;
    Curve curve1(my_curve1->size());
    Curve curve2(my_curve1->size());
    
    for (auto coord : *my_curve1) {
        Point point(1);
        point.set(0, coord);
        curve1.push_back(point);
    }

    for (auto coord : *my_curve2) {
        Point point(1);
        point.set(0, coord);
        curve2.push_back(point);
    }

    distance = Frechet::Continuous::distance(curve1, curve2).value;


    return distance;
}


//just keep the greater value from dataset on stack 
double value_for_padding(double m)
{
    static double M = 0.0;
    if (M < m)
        M = m;
    return M;
}


// padding the snapped_curve with the greater element M
void padding(std::vector<Type>* snapped_curve, int D)
{
    //otherwise get the max number from dataset;
    double M = value_for_padding(0.0) + 1000;  // get the greater value from coordinates
    vector<Type>* final_vector;

    int iterations = D - snapped_curve->size();   // we must push iterations M elements in the vector 

    for (int i = 0; i < iterations; i++)
    {
        snapped_curve->push_back(M);
    }
}


void filtering(vector<Type>* curve)
{
    double e;
    double average_change = 0.0;
    double curr_change = 0.0;
    int index = 0;

    for (int i = 0; i < curve->size()-1; i++) {
        curr_change = abs(curve->at(i) - curve->at(i+1));
        average_change += curr_change;
        index++;
    }
    e = average_change / curve->size();
    e = 2*e;

    vector<double>::iterator it =  curve->begin();
    
    while(((it+1) != curve->end()) && (it != curve->end()) && ((it+2) != curve->end()))  
    {
        // if a, b, c is consecutive points
        // if |a - b| < e and |b - c| < e, then remove b

        if ((abs((*it) - (*(it+1))) < e) &&  (abs((*(it+1)) - (*(it+2))) < e))
            curve->erase(it+1);
        
        it ++;
    }

    it = curve->end()-1;        // go to the last element
    
    // if |a - b| < e and |b - c| < e, then remove b
    if ((abs((*it-2) - (*(it-1))) < e) &&  (abs((*(it-1)) - (*(it))) < e))
        curve->erase(it-1);
    it++;
}


