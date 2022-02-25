#include "Grid.h"
#include <iostream>

using namespace std;


template <typename T>
Grid<T>::Grid(double delta)
{
    this->delta = delta;

    // initialize the t with uniform distribution
    std::random_device rd;      // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd());     // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0, this->delta);
    this->t = dis(gen);
}


template <typename T>
vector<T>* Grid<T>::get_snapped_curve_1D(vector<T>* curve)
{
    vector<T>* snapped_curve = new vector<T>();     // initialize a vector

    vector<T>* curve_vector = curve;
    T current_coordinate;
    double grid_coordinate;

    // snapping
    for (int i = 0; i < curve_vector->size(); i++) 
    {   
        current_coordinate = curve_vector->at(i);
        grid_coordinate = floor(((current_coordinate / this->delta) + 0.5)) * this->delta;  // see  lecture 13/12/2021
        snapped_curve->push_back(grid_coordinate);
    }

    this->MinimaMaxima(snapped_curve); // remove duplicats coordinates

    return snapped_curve;
}


template <typename T>
vector<T>* Grid<T>::get_snapped_curve_2D(vector<T>* curve)
{
    vector< pair<T, T> >* snapped_curve_2D = new vector< pair<T, T> >();     // initialize a vector
    vector<T>* snapped_curve_1D = new vector<T> ();
    pair<T, T> curr_coordinate;

    vector<T>* curve_vector = curve;
    T current_coordinate_x, current_coordinate_y;
    
    double grid_coordinate_x, grid_coordinate_y;
    
    int time = 1;
    int sampling_frequency = 1;   // is posible this to change with different dataset

    // snapping
    for (int i = 0; i < curve_vector->size(); i++) 
    {   
        current_coordinate_x = time;
        current_coordinate_y = curve_vector->at(i);

        // see page 17 from 2.curves.pdf      - Pij -> floor(x/δ + 1/2) * δ -
        grid_coordinate_x = floor(((abs(current_coordinate_x-t)) / this->delta) + 0.5) * this->delta;    
        grid_coordinate_y = floor(((abs(current_coordinate_y-t)) / this->delta) + 0.5) * this->delta;  
        
        // shift the curve by t
        grid_coordinate_x = grid_coordinate_x + this->t;        
        grid_coordinate_y = grid_coordinate_y + this->t; 


        snapped_curve_2D->push_back( make_pair(grid_coordinate_x, grid_coordinate_y) );

        time = time + sampling_frequency; // we go to the next coordinate so next day in dataset
    }

    
    // remove consecutive duplicates
    this->remove_consecutive_duplicates(snapped_curve_2D);


    // unpack the snapped curve to the vector, 
    // e.g: snapped_curve = [{1, 2}, {2, 5}] - after unpack -> [1, 2, 2, 5]
    for (int i = 0; i < snapped_curve_2D->size(); i++)
    {
        curr_coordinate = snapped_curve_2D->at(i);  // get the 2D coordinate
        snapped_curve_1D->push_back(curr_coordinate.first);     // insert the x
        snapped_curve_1D->push_back(curr_coordinate.second);    // insert the y
    }

    delete snapped_curve_2D;

    return snapped_curve_1D;    
}


// for 2D grid, e.g { (1, 1), (1, 1), (2, 2), (1, 1) } ->  {(1, 1), (2, 2), (1, 1)}
template <typename T>
void Grid<T>::remove_consecutive_duplicates(std::vector< std::pair<T, T> >* snapped_curve)
{
    vector< pair<Type, Type> >::iterator it =  snapped_curve->begin();



    while((it+1) != snapped_curve->end() && (it != snapped_curve->end()))  
    {
        if ((*it) == (*(it+1)))     // compare the c(i) with c(i+1)   
        {   
           snapped_curve->erase(it);
        }
        it += 1;
    }

    it = snapped_curve->end()-1;
    if ((*(it-1)) == (*it))     // for the 2 last pairs
    {
        snapped_curve->erase(it++);
    }
}


// for 1D grid, remove middle point between 3 consecutive points in curve
template <typename T>
void Grid<T>::MinimaMaxima(std::vector< T >* snapped_curve)
{
    vector<Type>::iterator it =  snapped_curve->begin();
    double minima, maxima;


    while((it != snapped_curve->end()) && ((it+1) != snapped_curve->end())  &&  ((it+2) != snapped_curve->end()))  
    {
        minima = std::min({ (*it), *(it+2) });
        maxima = std::max({ (*it), *(it+2) });

        if ((minima <= *(it+1)) && (*(it+1) <= maxima))     // compare the c(i) with c(i+1)   
        {   
           snapped_curve->erase(it+1);
        }
        
        it++;
    }

    it = snapped_curve->end()-1;    // go to the last element
    minima = std::min({ (*it-2), *(it) });
    maxima = std::max({ (*it-2), *(it) });

    if ((minima <= *(it-1)) && (*(it-1) <= maxima))     // for the 2 last pairs
    {
        snapped_curve->erase(it-1);
    }
    it++;
}

