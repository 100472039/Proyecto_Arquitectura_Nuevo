//
// Created by anamaria on 11/11/23.
//

#include "bloque.hpp"
#include "grid.hpp"
#include "../constantes.hpp"
#include <iostream>
#include <cmath>
#include <vector>
#include <set>

int ajustarLim(int value, int lowerLimit, int  upperLimit) {
    if (value < lowerLimit) {
        return lowerLimit;
    } if (value >= upperLimit) {
        return upperLimit-1;
    }
    return value;
}

void anadir_particulas(std::vector<Block>& bloques , std::vector<Particle>& particles){
    for (const Particle &particle : particles) {
        // Calcular los índices de bloque para la partícula
        int pos_i = ajustarLim(std::floor((particle.px - xmin) / sx), 0, nx);
        int pos_j = ajustarLim(std::floor((particle.py - ymin) / sy), 0, ny);
        int pos_k = ajustarLim(std::floor((particle.pz - zmin) / sz), 0, nz);

        bloques[(pos_k * ny * nx) + (pos_j * nx) + pos_i].particles.push_back(particle);
        /*std::cout << "Particle at (" << particle.px << ", " << particle.py << ", " << particle.pz
                  << ") added to block (" << pos_i<< ", " << pos_j << ", " << pos_k << ")\n";*/
    }
}




bool isValidIndices(int i, int j, int k) {
    return i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && k < nz;
}

void anotar_adyacentes(std::vector<Block>& bloques) {
    // Bucle de todos los bloques
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                Block& currentBlock = bloques[(k * ny * nx) + (j * nx) + i];

                // Todas las posibles posiciones de bloques adyacentes
                for (int di = -1; di <= 1; ++di) {
                    for (int dj = -1; dj <= 1; ++dj) {
                        for (int dk = -1; dk <= 1; ++dk) {
                            // Comprobar que es un bloque adyacente y no el mismo bloque
                            if (isValidIndices(i + di, j + dj, k + dk)&&!(di == 0 && dj == 0 && dk == 0)){
                                currentBlock.adyacentes.push_back(bloques[((k + dk) * ny * nx) + ((j + dj) * nx) + (i + di)].id);
                            }
                        }
                    }
                }
            }
        }
    }
}
void recorrer_adyacentes(std::vector<Block>& bloques) {
    for (const Block &bloque : bloques) {
        std::cout << "\nBloque actual: (" << bloque.id << ")     ";
        std::cout << "Bloques adyacentes:";

        for (int adyacente: bloque.adyacentes) {
            std::cout << "  (" << adyacente<< ") ";
            // Realiza las operaciones que desees con el bloque adyacente
        }
    }
}

void crearBloques(std::vector<Block>& bloques,std::vector<Particle>& particles){
    // Primero, creamos los bloques
    int id = 0;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                Block block;
                block.id = id++;
                bloques.push_back(block);
            }
        }
    }
    anotar_adyacentes(bloques);
    recorrer_adyacentes(bloques);

}



void generar_malla(std::vector<Particle>& particles) {
    Grid grid;
    std::vector<Block> bloques;
    crearBloques(bloques,particles);
}