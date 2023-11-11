// En grid.cpp
#include "grid.hpp"
#include "../constantes.hpp"
#include "progargs.hpp"
#include <cmath>
#include <iostream>

Grid::Grid() {
    // Calcular los tamaños de malla nx, ny y nz
    const double nx = std::floor((xmax - xmin) / smooth);
    const double ny = std::floor((ymax - ymin) / smooth);
    const double nz = std::floor((zmax - zmin) / smooth);
    std::cout << "nxjjj: " << smooth << "\n";

    // Inicializar el vector n con los valores calculados
    n = {nx, ny, nz};
}

std::vector<double> Grid::getN() {
    return n;
}

std::vector<double> Grid::calculateBlockSize() {
    // Calcular el tamaño del bloque (sx, sy, sz) basado en n
    const double nx = std::floor(n[0]);
    const double ny = std::floor(n[1]);
    const double nz = std::floor(n[2]);

    std::cout << "nx: " << nx << "\n";
    std::cout << "ny: " << ny << "\n";
    std::cout << "nz: " << nz << "\n";

    const double sx = (xmax - xmin) / nx;
    const double sy = (ymax - ymin) / ny;
    const double sz = (zmax - zmin) / nz;

    // Aquí puedes hacer algo con las partículas y aceleración, si es necesario

    return {sx, sy, sz};
}


