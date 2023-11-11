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





void anotar_adyacentes(std::vector<Block>& bloques) {
    //Bucle de todos los bloques
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                Block& currentBlock = bloques[(k * ny * nx) + (j * nx) + i];

                //Todas las posibles posiciones de bloques adyacentes
                for (int di = -1; di <= 1; ++di) {
                    if (di+i> 0 && di+i<nx){
                        for (int dj = -1; dj <= 1; ++dj) {
                            if (dj+j> 0 && dj+j<ny) {
                                for (int dk = -1; dk <= 1; ++dk) {
                                    if (dk+k> 0 && dk+k<nz) {
                                        //Comprobar que es un bloque adyacente
                                        Block &neighborBlock = bloques[((k + dk) * ny * nx) + ((j + dj) * nx) +
                                                                       (i + di)];
                                        currentBlock.adyacentes.push_back(neighborBlock);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void crearBloques(std::vector<Block>& bloques,std::vector<Particle>& particles){
    // Primero, creamos los bloques
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                Block block;
                bloques.push_back(block);
            }
        }
    }
    anotar_adyacentes(bloques);
    for (int i=0; i<1; i++) {
    //anadir_particulas(bloques, particles);
    //incrementar_dens(bloques, particles);
    }
}



void generar_malla(std::vector<Particle>& particles) {
    Grid grid;
    std::vector<Block> bloques;
    crearBloques(bloques,particles);
}