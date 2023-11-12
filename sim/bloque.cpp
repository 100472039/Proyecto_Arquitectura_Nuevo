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

    for (int i=0;i<4800;i++) {

        // Calcular los índices de bloque para la partícula
        int pos_i = ajustarLim(std::floor((particles[i].px - xmin) / sx), 0, nx);
        int pos_j = ajustarLim(std::floor((particles[i].py - ymin) / sy), 0, ny);
        int pos_k = ajustarLim(std::floor((particles[i].pz - zmin) / sz), 0, nz);


        bloques[(pos_k * ny * nx) + (pos_j * nx) + pos_i].particles.push_back(particles[i]);
    }
}
double diferenciaDistancias(const Particle& particle_i, const Particle& particle_j) {
    double distx = particle_i.px - particle_j.px;
    double disty = particle_i.py - particle_j.py;
    double distz = particle_i.pz - particle_j.pz;

    // Calcula la distancia euclidiana entre las partículas
    double distancia = (distx * distx )+ (disty * disty )+ (distz * distz);

    return distancia;
}

bool distanciaMenosH2(const Particle& particle_i, const Particle& particle_j) {
    double distx = particle_i.px - particle_j.px;
    double disty = particle_i.py - particle_j.py;
    double distz = particle_i.pz - particle_j.pz;
    double distanciaCuadrada = (distx * distx )+ (disty * disty )+ (distz * distz);

    return distanciaCuadrada < smooth * smooth;
}


bool isValidIndices(int i, int j, int k) {
    return i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && k < nz;
}
/*
void anotar_adyacentes(std::vector<Block>& bloques) {
    // Bucle de todos los bloques
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                Block& currentBlock = bloques[(k * ny * nx) + (j * nx) + i];

                // Determinar los límites de las iteraciones para los bloques adyacentes
                int min_di = (i == 0) ? 0 : -1;
                int max_di = (i == nx - 1) ? 0 : 1;
                int min_dj = (j == 0) ? 0 : -1;
                int max_dj = (j == ny - 1) ? 0 : 1;
                int min_dk = (k == 0) ? 0 : -1;
                int max_dk = (k == nz - 1) ? 0 : 1;

                // Todas las posibles posiciones de bloques adyacentes
                for (int di = min_di; di <= max_di; ++di) {
                    for (int dj = min_dj; dj <= max_dj; ++dj) {
                        for (int dk = min_dk; dk <= max_dk; ++dk) {
                            Block &currentBlock = bloques[(k * ny * nx) + (j * nx) + i];
                            // Comprobar que es un bloque adyacente y no el mismo bloque
                            if (isValidIndices(i + di, j + dj, k + dk)) {
                                if (bloques[((k + dk) * ny * nx) + ((j + dj) * nx) + (i + di)].id >=
                                    currentBlock.id) {
                                    currentBlock.adyacentes.push_back(
                                            bloques[((k + dk) * ny * nx) + ((j + dj) * nx) + (i + di)].id);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
*/
void incrementar_densidades(std::vector<Block>& bloques, std::vector<Particle>& particles) {
    const double h2 = smooth * smooth;
    int contador = 0;

    int index = 0;
    for (int index = 0; index < nx*ny*nz; index ++) {
        Block &currentBlock = bloques[index];

        int i = index % nx;
        int j = (index / nx) % ny;
        int k = index / (nx * ny);

        // Determine the limits for iterations over adjacent blocks
        int min_di = (i == 0) ? 0 : -1;
        int max_di = (i == nx - 1) ? 0 : 1;
        int min_dj = (j == 0) ? 0 : -1;
        int max_dj = (j == ny - 1) ? 0 : 1;
        int min_dk = (k == 0) ? 0 : -1;
        int max_dk = (k == nz - 1) ? 0 : 1;

        for (int di = min_di; di <= max_di; ++di) {
            for (int dj = min_dj; dj <= max_dj; ++dj) {
                for (int dk = min_dk; dk <= max_dk; ++dk) {
                    Block &neighbor = bloques[(k + dk) * ny * nx + (j + dj) * nx + (i + di)];
                    {
                        Block &neighbor = bloques[(k + dk) * ny * nx + (j + dj) * nx + (i + di)];

                        if (neighbor.id <= currentBlock.id) {
                            for (Particle &particle_i: currentBlock.particles) {
                                for (Particle &particle_j: neighbor.particles) {
                                    if (particle_i.id != particle_j.id) {
                                        if (distanciaMenosH2(particle_i, particle_j)) {
                                            double dist = diferenciaDistancias(particle_i, particle_j);
                                            double dens_increment = (h2 - dist) * (h2 - dist) * (h2 - dist);

                                            particle_i.rho += dens_increment;
                                            particle_j.rho += dens_increment;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                contador++;
            }
        }

    }
    std::cout << contador << "\n";
}


//BING
/*void incrementar_densidades(std::vector<Block>& bloques,std::vector<Particle>& particles) {
    double incrementar_densidad = 0;
    double h2 = (smooth*smooth);
    double h6 = h2 * h2 * h2; // Precomputar h^6

    // Crear un vector para almacenar las distancias precomputadas
    std::vector<std::vector<double>> distancias(particles.size(), std::vector<double>(particles.size(), 0));

    // Precomputar todas las distancias
    for (int i = 0; i < particles.size(); ++i) {
        for (int j = i+1; j < particles.size(); ++j) {
            distancias[i][j] = distanciaMenosH2(particles[i], particles[j]);
            distancias[j][i] = distancias[i][j]; // Aprovechar la simetría para reducir los cálculos
        }
    }

    for (const Block &bloque: bloques){
        int currentBlock = bloque.id;
        Block& currenBlock = bloques[bloque.id];
        for (int adyacente: bloque.adyacentes) {
            Block& bloque_adyacente = bloques[adyacente];
            for (Particle& particula_i : currenBlock.particles){
                for (Particle& particula_j : bloque_adyacente.particles){
                    if (distancias[particula_i.id][particula_j.id]) {
                        incrementar_densidad = h6 - distancias[particula_i.id][particula_j.id];
                        incrementar_densidad *= incrementar_densidad * incrementar_densidad;
                    }
                    else{
                        incrementar_densidad = 0;
                    }

                    particula_i.rho += incrementar_densidad;
                    particula_j.rho += incrementar_densidad;

                }

            }
        }
    }
}*/



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
    //anotar_adyacentes(bloques);
    //recorrer_adyacentes(bloques);
    for(int i=0;i<20;i++){
        anadir_particulas(bloques,particles);
        incrementar_densidades(bloques,particles);
    }
}



void generar_malla(std::vector<Particle>& particles) {
    Grid grid;
    std::vector<Block> bloques;
    crearBloques(bloques,particles);
}