//
// Created by natha on 30/10/2023.
//

#ifndef PROYECTO_ARQUITECTURA_BLOCK_HPP
#define PROYECTO_ARQUITECTURA_BLOCK_HPP
#include "grid.hpp"
#include <vector>
#include "progargs.hpp"



//extern std::vector<std::vector<std::vector<std::vector<Particle>>>> particleGrid;

int ajustarLimite(int value, int lowerLimit, int upperLimit);

int main_block(std::vector<Particle>& particles);
#endif //PROYECTO_ARQUITECTURA_BLOCK_HPP

