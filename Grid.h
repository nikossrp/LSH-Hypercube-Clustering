#ifndef __GRID__
#define __GRID__

#include "Utilities.h"

template<typename T>
class Grid {

    private:
        double delta;
        double t;       // uniform coordinate in range [0, delta)

    public:
        Grid(double delta);   // initialize grid
        ~Grid() { }
        std::vector<T>* get_snapped_curve_1D(std::vector<T>* curve);  
        std::vector<T>* get_snapped_curve_2D(std::vector<T>* curve);                  // the second dimension is time [1, 2, 3,.... ] return an 1D curve
        void remove_consecutive_duplicates(std::vector< std::pair<T, T> >* snapped_curve);     // remove consecutive duplicates for 2D grid - discrete curve - 
        void MinimaMaxima(std::vector< T >* snapped_curve);          // remove middle point find minma and maxima for 1D grid  - continuous curve -

};

#endif  // __GRID__
