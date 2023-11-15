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

double ajustarLim(double value, double lowerLimit, double upperLimit) {
    if (value < lowerLimit) {
        return lowerLimit;
    } if (value >= upperLimit) {
        return upperLimit-1;
    }
    return value;
}
double calcularDistij(const Particle& particle_i, const Particle& particle_j) {
    double distx = particle_i.px - particle_j.px;
    double disty = particle_i.py - particle_j.py;
    double distz = particle_i.pz - particle_j.pz;

    // Calcula la distancia euclidiana entre las partículas
    double distanciaCuadrada = (distx * distx )+ (disty * disty )+ (distz * distz);

    double distij = std::max(distanciaCuadrada, epsilon);
    double distij2 = std::sqrt(distij);
    //std::cout <<"2: "<<distanciaCuadrada<<"     ";

    return distij2;
}
void anadir_particulas(std::vector<Block>& bloques , std::vector<Particle>& particles){

    for (int i=0;i<4800;i++) {

        // Calcular los índices de bloque para la partícula
        double pos_i = ajustarLim(std::floor((particles[i].px - xmin) / sx), 0, nx);
        double pos_j = ajustarLim(std::floor((particles[i].py - ymin) / sy), 0, ny);
        double pos_k = ajustarLim(std::floor((particles[i].pz - zmin) / sz), 0, nz);


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
                            if (isValidIndices(i + di, j + dj, k + dk)) {
                                //if (currentBlock.id<=bloques[((k + dk) * ny * nx) + ((j + dj) * nx) + (i + di)].id){
                                bloques[(k * ny * nx) + (j * nx) + i].adyacentes.push_back(bloques[((k + dk) * ny * nx) + ((j + dj) * nx) + (i + di)].id);}
                            //}
                        }
                    }
                }
            }
        }
    }
}



//Este da bien 25s
void incrementar_densidades(std::vector<Block>& bloques, std::vector<Particle>& particles) {
    const double h2 = smooth * smooth;
    int contador = 0;
    double dens_increment = 0;

    for (Block& currentBlock : bloques) {
        for (Particle& particle_i : currentBlock.particles) {
            for (double neighbor : currentBlock.adyacentes) {
                Block &siguiente = bloques[neighbor];
                for (Particle& particle_j : siguiente.particles) {

                    if (particle_i.id>particle_j.id){
                        if (distanciaMenosH2(particle_i, particle_j)) {
                            double dist = diferenciaDistancias(particle_i, particle_j);
                            dens_increment = (h2 - dist) * (h2 - dist) * (h2 - dist);
                        }else{
                            dens_increment= 0;
                        }
                        particles[particle_i.id].rho+= dens_increment;
                        particles[particle_j.id].rho +=dens_increment;

                        particle_i.rho += dens_increment;
                        particle_j.rho += dens_increment;
                    }

                }

            }
        }
    }
}


void transformar_densidades(std::vector<Particle>& particles){
    double parte2 = numero_dens1/(numero_dens2*pi*(smooth*smooth*smooth*smooth*smooth*smooth*smooth*smooth*smooth));
    for (Particle& particle : particles) {
        double parte1 = particle.rho + std::pow(smooth, 6);
        particle.rho = parte1 * parte2 * masa;
        //std::cout <<"densidad_particula: "<<particle.rho<<"\n ";
    }
}

void transferencia_aceleracion(std::vector<Block>& bloques, std::vector<Particle>& particles) {
    double fraccion2 = quarentaycinco / (pi * smooth*smooth*smooth*smooth*smooth*smooth);
    double fraccion3 = quince / (pi * smooth*smooth*smooth*smooth*smooth*smooth);
    double fraccion4 = (3 * masa * presion_de_rigidez) / 2;
    for (Block& currentBlock : bloques) {
        for (Particle& particle_i : currentBlock.particles) {
            for (double neighbor: currentBlock.adyacentes) {
                Block &siguiente = bloques[neighbor];
                for (Particle &particle_j: siguiente.particles) {
                    if (particle_i.id > particle_j.id) {
                        if (distanciaMenosH2(particle_i, particle_j)) {
                            double distij = calcularDistij(particle_i, particle_j);

                            double fraccion1 = ((smooth - distij) * (smooth - distij)) / distij;

                            double incremento_aceleracionx =
                                    (((particle_i.px - particle_j.px) * fraccion3 * fraccion4 * fraccion1 * (particle_i.rho + particle_j.rho - 2 * densidad)) +
                                    ((particle_j.vx - particle_i.vx) * fraccion2 * viscosidad * masa))/ (particle_i.rho * particle_j.rho);
                            double incremento_aceleraciony =
                                    (((particle_i.py - particle_j.py) * fraccion3 * fraccion4 * fraccion1 * (particle_i.rho + particle_j.rho - 2 * densidad)) +
                                    ((particle_j.vy - particle_i.vy) * fraccion2 * viscosidad * masa))/(particle_i.rho * particle_j.rho);

                            double incremento_aceleracionz = (((particle_i.pz - particle_j.pz) * fraccion3 * fraccion4 * fraccion1 * (particle_i.rho + particle_j.rho - 2 * densidad)) +
                                                               ((particle_j.vz - particle_i.vz) * fraccion2 * viscosidad * masa)) / (particle_i.rho * particle_j.rho);

                            particle_i.aceleracion_externa.x += incremento_aceleracionx;
                            particle_i.aceleracion_externa.y += incremento_aceleraciony;
                            particle_i.aceleracion_externa.z += incremento_aceleracionz;

                            particle_j.aceleracion_externa.x -= incremento_aceleracionx;
                            particle_j.aceleracion_externa.y -= incremento_aceleraciony;
                            particle_j.aceleracion_externa.z -= incremento_aceleracionz;

                            particles[particle_i.id].aceleracion_externa.x = particle_i.aceleracion_externa.x;
                            particles[particle_j.id].aceleracion_externa.x = particle_j.aceleracion_externa.x;
                            particles[particle_i.id].aceleracion_externa.y = particle_i.aceleracion_externa.y;
                            particles[particle_j.id].aceleracion_externa.y = particle_j.aceleracion_externa.y;
                            particles[particle_i.id].aceleracion_externa.z = particle_i.aceleracion_externa.z;
                            particles[particle_j.id].aceleracion_externa.z = particle_j.aceleracion_externa.z;
                        }
                    }
                }
            }
        }
    }for (Particle& particle : particles) {
        std::cout <<"aceleracionx: "<<particle.aceleracion_externa.x<<" ";
        std::cout <<"aceleraciony: "<<particle.aceleracion_externa.y<<" ";
        std::cout <<"aceleracionz: "<<particle.aceleracion_externa.z<<"\n ";

    }
}



void eliminar_particulas(std::vector<Block>& bloques) {
    for (Block &bloque : bloques) {
        //std::cout << "Bloque ID: " << bloque.id << std::endl;
        for (Particle &particle : bloque.particles) {
            //std::cout << "  Particle ID: " << particle.id << ", rho: " << particle.rho << std::endl;
        }
        bloque.particles.clear();
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
    //recorrer_adyacentes(bloques);
    for(int i=0;i<1;i++){
        anadir_particulas(bloques,particles);
        incrementar_densidades(bloques,particles);
        transformar_densidades(particles);
        transferencia_aceleracion(bloques,particles);
        eliminar_particulas(bloques);
    }
}


void generar_malla(std::vector<Particle>& particles) {
    Grid grid;
    std::vector<Block> bloques;
    crearBloques(bloques,particles);
}