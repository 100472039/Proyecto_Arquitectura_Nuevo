// En grid.cpp
#include "grid.hpp"
#include "../constantes.hpp"
#include "progargs.hpp"
#include <cmath>
#include <iostream>

Grid::Grid() {
    // Calcular los tamaños de malla nx, ny y nz
    const int n_x = std::floor((xmax - xmin) / smooth);
    const int n_y = std::floor((ymax - ymin) / smooth);
    const int n_z = std::floor((zmax - zmin) / smooth);
    ::nx = n_x;
    ::ny = n_y;
    ::nz = n_z;
    const double s_x = (xmax - xmin) / nx;
    const double s_y = (ymax - ymin) / ny;
    const double s_z = (zmax - zmin) / nz;
    ::sx = s_x;
    ::sy = s_y;
    ::sz = s_z;
    // Inicializar el vector n con los valores calculados
    n = {nx, ny, nz};
}

std::vector<int> Grid::getN() {
    return n;
}




