// En grid.hpp
#ifndef GRID_HPP
#define GRID_HPP

#include <vector>
#include "progargs.hpp"




class Grid {
public:
    Grid(); // Constructor

    // Método para obtener el vector n
    std::vector<int> getN();
    std::vector<double> calculateBlockSize();

private:
    std::vector<int> n;
};
int nx;
int ny;
int nz;
double sx;
double sy;
double sz;


#endif
