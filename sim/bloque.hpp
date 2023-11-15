//
// Created by anamaria on 11/11/23.
//
#include <cmath>
#ifndef FLUID_BLOQUE_HPP
#define FLUID_BLOQUE_HPP
#include "progargs.hpp"
#include "grid.hpp"

#pragma once



struct Block {
    double id;  // Identificador único del bloque
    int i, j, k;  // Índices del bloque
    std::vector<Particle> particles;
    std::vector<double> adyacentes;  // Almacena los identificadores de los bloques adyacentes
};


#endif //FLUID_BLOQUE_HPP