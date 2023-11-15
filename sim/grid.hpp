// En grid.hpp
#ifndef GRID_HPP
#define GRID_HPP

#include <vector>
#include "progargs.hpp"




class Grid {
public:
    Grid(); // Constructor

    // Método para obtener el vector n
    std::vector<double> getN();
    std::vector<double> calculateBlockSize();

private:
    std::vector<double> n;
};
double nx;
double ny;
double nz;
double sx;
double sy;
double sz;


#endif
