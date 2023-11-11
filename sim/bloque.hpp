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
    int i, j, k;  // índices del bloque
    std::vector<Particle> particles;
    std::vector<Block> adyacentes;
};


#endif //FLUID_BLOQUE_HPP