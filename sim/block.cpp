#include "grid.hpp"
#include "../constantes.hpp"
#include <iostream>
#include <cmath>
#include "progargs.hpp"
#include "block.hpp"
#include <vector>
#include <set>

/*
struct ParticleEntry {
        int cx, cy, cz;
        Particle particle;
};

bool distanciaMenosH2(const Particle& particle_i, const Particle& particle_j) {
    double distx = particle_i.px - particle_j.px;
    double disty = particle_i.py - particle_j.py;
    double distz = particle_i.pz - particle_j.pz;
    double distanciaCuadrada = (distx * distx )+ (disty * disty )+ (distz * distz);

    return distanciaCuadrada < smooth * smooth;
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

double diferenciaDistancias(const Particle& particle_i, const Particle& particle_j) {
    double distx = particle_i.px - particle_j.px;
    double disty = particle_i.py - particle_j.py;
    double distz = particle_i.pz - particle_j.pz;

    // Calcula la distancia euclidiana entre las partículas
    double distancia = (distx * distx )+ (disty * disty )+ (distz * distz);

    return distancia;
}

int ajustarLimite(int value, int lowerLimit, int  upperLimit) {
    if (value < lowerLimit) {
        return lowerLimit;
    } if (value >= upperLimit) {
        return upperLimit-1;
    }
    return value;
}


int main_block(std::vector<Particle>& particles) {
    Grid grid;


    // Accede a los valores de nx, ny, nz a través del método getN()
    std::vector<int> n = grid.getN();
    std::vector<double> blockSize = grid.calculateBlockSize();

    std::vector<ParticleEntry> particleGrid;
    std::set<std::pair<int, int>> processed_pairs;
    const double nx = std::floor(n[0]);
    const double ny = std::floor(n[1]);
    const double nz = std::floor(n[2]);

    const double n_bloques = nx * ny * nz;


    const double sx = blockSize[0];
    const double sy = blockSize[1];
    const double sz = blockSize[2];

    std::cout << "n_bloques: " << n_bloques << "\n";
    std::cout << "sx: " << sx << "\n";
    std::cout << "sy: " << sy << "\n";
    std::cout << "sz: " << sz << "\n";

    for (Particle& particle : particles) {
        int i = static_cast<int>(std::floor((particle.px - xmin) / sx));
        int j = static_cast<int>(std::floor((particle.py - ymin) / sy));
        int k = static_cast<int>(std::floor((particle.pz - zmin) / sz));


        // Asegurar de que i, j y k estén dentro de los límites.
        i = ajustarLimite(i, 0, nx);
        j = ajustarLimite(j, 0, ny);
        k = ajustarLimite(k, 0, nz);

        // Agrega la partícula a la matriz tridimensional en las coordenadas i, j, k correspondientes.
        particleGrid.push_back({i, j, k, particle});
    }
    double incrementar_densidad = 0;

    //Incremento de densidades OK!!
    for (std::vector<ParticleEntry>::size_type i = 0; i < particleGrid.size(); ++i) {
        ParticleEntry& entry_i = particleGrid[i];
        Particle& particle_i = entry_i.particle;
        for (std::vector<ParticleEntry>::size_type j = i + 1; j < particleGrid.size(); ++j) {
            ParticleEntry& entry_j = particleGrid[j];
            Particle& particle_j = entry_j.particle;
            // Solo procesa el par de partículas si son contiguas.
            if (abs(entry_i.cx - entry_j.cx) <= 1 &&
                abs(entry_i.cy - entry_j.cy) <= 1 &&
                abs(entry_i.cz - entry_j.cz) <= 1) {*/
                /*if (contador==0){in
                    std::cout <<"contador: "<< entry_i.cx << entry_i.cy << entry_i.cz << ", " << entry_j.cx << entry_j.cy << entry_j.cz<< std::endl;}*/
                //Este if SOBRA
                /*if (processed_pairs.find({i, j}) == processed_pairs.end()) {
                    if (distanciaMenosH2(particle_i, particle_j)) {
                        incrementar_densidad = diferenciaDistancias(particle_i, particle_j);
                        //std::cout <<"2: "<<incrementar_densidad<<"     ";
                        double h2 = smooth*smooth;
                        incrementar_densidad = (h2- incrementar_densidad) * (h2- incrementar_densidad) * (h2- incrementar_densidad);
                    }
                    else{
                        incrementar_densidad = 0;
                    }
                    //std::cout <<"1: "<<particle_i.rho <<"\n";
                    particle_i.rho += incrementar_densidad;
                    particle_j.rho += incrementar_densidad;

                    particles[i].rho = particle_i.rho;
                    particles[j].rho = particle_j.rho;

                    // Agrega el par de partículas al conjunto de pares procesados.
                    processed_pairs.insert({i, j});
                    processed_pairs.insert({j, i});
                }
            }
        }
    }
    processed_pairs.clear();


    for (Particle& particle : particles) {
        std::cout <<"densidad_particula: "<<particle.rho<<"\n ";

    }

    //Tranformacion de densidad.OK!!
    for (Particle& particle : particles) {
        double parte2 = numero_dens1/(numero_dens2*pi*std::pow(smooth, 9));
        double parte1 = particle.rho + std::pow(smooth, 6);
        particle.rho = parte1 * parte2 * masa;
        //std::cout <<"densidad_particula: "<<particle.rho<<"\n ";
    }

    particleGrid.clear();
    for (Particle& particle : particles) {
        int i = static_cast<int>(std::floor((particle.px - xmin) / sx));
        int j = static_cast<int>(std::floor((particle.py - ymin) / sy));
        int k = static_cast<int>(std::floor((particle.pz - zmin) / sz));

        // Asegurar de que i, j y k estén dentro de los límites.
        i = ajustarLimite(i, 0, nx);
        j = ajustarLimite(j, 0, ny);
        k = ajustarLimite(k, 0, nz);

        // Agrega la partícula a la matriz tridimensional en las coordenadas i, j, k correspondientes.
        particleGrid.push_back({i, j, k, particle});
    }


    //Transferencia de aceleracion.OK!!
    for (std::vector<ParticleEntry>::size_type i = 0; i < particleGrid.size(); ++i) {
        ParticleEntry &entry_i = particleGrid[i];
        Particle &particle_i = entry_i.particle;*/
        /*
        if (i > 10) {
            std::cout << "posicion: " << particle_i.px << "," << particle_i.py << "," << particle_i.pz << std::endl;
        }*/
        /*for (std::vector<ParticleEntry>::size_type j = i + 1; j < particleGrid.size(); ++j) {
            ParticleEntry &entry_j = particleGrid[j];
            Particle &particle_j = entry_j.particle;
            // Solo procesa el par de partículas si son contiguas.
            if (abs(entry_i.cx - entry_j.cx) <= 1 &&
                abs(entry_i.cy - entry_j.cy) <= 1 &&
                abs(entry_i.cz - entry_j.cz) <= 1) {
                if (processed_pairs.find({i, j}) == processed_pairs.end()) {
                    if (distanciaMenosH2(particle_i, particle_j)) {
                        double distij = calcularDistij(particle_i,particle_j);
                        double aceleracion_px = particle_i.px - particle_j.px;
                        double aceleracion_py = particle_i.py - particle_j.py;
                        double aceleracion_pz = particle_i.pz - particle_j.pz;
                        //double aceleracion = aceleracion_px+aceleracion_py+aceleracion_pz;
                        double velocidad_px = particle_j.vx - particle_i.vx;
                        double velocidad_py = particle_j.vy - particle_i.vy;
                        double velocidad_pz = particle_j.vz - particle_i.vz;
                        //double velocidad= velocidad_px+velocidad_py+velocidad_pz;
                        double fraccion1 = ((smooth-distij)*(smooth-distij))/distij;
                        double fraccion2 = quarentaycinco/(pi*std::pow(smooth, 6));
                        double fraccion3 = quince/(pi*std::pow(smooth, 6));
                        double fraccion4 = (3*masa*presion_de_rigidez)/2;
                        double parte1 = particle_i.rho + particle_j.rho - 2*densidad;
                        double parte2 = particle_i.rho * particle_j.rho;

                        double incremento_aceleracionx = (aceleracion_px*fraccion3*fraccion4*fraccion1*parte1)+(velocidad_px*fraccion2*viscosidad*masa);
                        double incremento_aceleracion_x = incremento_aceleracionx/parte2;
                        double incremento_aceleraciony = (aceleracion_py*fraccion3*fraccion4*fraccion1*parte1)+(velocidad_py*fraccion2*viscosidad*masa);
                        double incremento_aceleracion_y = incremento_aceleraciony/parte2;
                        double incremento_aceleracionz = (aceleracion_pz*fraccion3*fraccion4*fraccion1*parte1)+(velocidad_pz*fraccion2*viscosidad*masa);
                        double incremento_aceleracion_z = incremento_aceleracionz/parte2;

                        particle_i.aceleracion_externa.x += incremento_aceleracion_x;
                        particle_i.aceleracion_externa.y += incremento_aceleracion_y;
                        particle_i.aceleracion_externa.z += incremento_aceleracion_z;
                        //std::cout <<"2: "<<incremento_aceleracion_x<<"     ";

                        particle_j.aceleracion_externa.x -= incremento_aceleracion_x;
                        particle_j.aceleracion_externa.y -= incremento_aceleracion_y;
                        particle_j.aceleracion_externa.z -= incremento_aceleracion_z;

                        particles[i].aceleracion_externa.x = particle_i.aceleracion_externa.x;
                        particles[j].aceleracion_externa.x = particle_j.aceleracion_externa.x;
                        particles[i].aceleracion_externa.y = particle_i.aceleracion_externa.y;
                        particles[j].aceleracion_externa.y = particle_j.aceleracion_externa.y;
                        particles[i].aceleracion_externa.z = particle_i.aceleracion_externa.z;
                        particles[j].aceleracion_externa.z = particle_j.aceleracion_externa.z;

                        processed_pairs.insert({i, j});
                        processed_pairs.insert({j, i});

                    }

                }
            }
        }
    }*/

    /*for (Particle& particle : particles) {
        std::cout <<"aceleracion: "<<particle.aceleracion_externa.z<<"\n ";
    }*/
    /*double contador = 1;
    for (Particle& particle : particles) {
        std::cout <<"adios: "<<contador<<" ";
        std::cout <<"aceleracionx: "<<particle.aceleracion_externa.x<<" ";
        std::cout <<"aceleraciony: "<<particle.aceleracion_externa.y<<" ";
        std::cout <<"aceleracionz: "<<particle.aceleracion_externa.z<<"\n ";
        contador++;
    }*/

//-----------------------------------------------------------------------------------------------------
    //Colisiones de partículas.OK!!
    /*double difx = 0;
    double dify = 0;
    double difz= 0;
    double contador = 1;
    //Versión optimizada
    for (std::vector<ParticleEntry>::size_type i = 0; i < particleGrid.size(); ++i) {
        ParticleEntry &entry_i = particleGrid[i];
        Particle &particle_i = entry_i.particle;
        double posicion_x = particle_i.px + particle_i.hvx*paso_de_tiempo;
        //Límite marcado por x
        if (entry_i.cx == 0) {
            difx = tamaño_de_particula - (posicion_x - xmin);
            if (difx > 10e-10) {
                particle_i.aceleracion_externa.x = particle_i.aceleracion_externa.x + ((colisiones_de_rigidez*difx) - (amortiguamiento*particle_i.vx));
                particles[i].aceleracion_externa.x = particle_i.aceleracion_externa.x;
                //std::cout <<"adios: "<<contador<<"\n";
            }
        }
        else if (entry_i.cx == nx-1) {
            difx = tamaño_de_particula - (xmax - posicion_x);
            if (difx > 10e-10) {
                particle_i.aceleracion_externa.x = particle_i.aceleracion_externa.x - ((colisiones_de_rigidez*difx) + (amortiguamiento*particle_i.vx));
                particles[i].aceleracion_externa.x = particle_i.aceleracion_externa.x;
                //std::cout <<"adios: "<<contador<<"\n";
            }
        }

        double posicion_y = particle_i.py + particle_i.hvy*paso_de_tiempo;
        //Límite marcado por y
        if (entry_i.cy == 0) {
            dify = tamaño_de_particula - (posicion_y  - ymin);
            if (dify > 10e-10) {
                particle_i.aceleracion_externa.y = particle_i.aceleracion_externa.y + (colisiones_de_rigidez*dify) - (amortiguamiento*particle_i.vy);
                particles[i].aceleracion_externa.y =  particle_i.aceleracion_externa.y;
            }
        }
        else if (entry_i.cy == ny-1) {
            dify = tamaño_de_particula - (ymax - posicion_y);
            if (dify > 10e-10) {
                particle_i.aceleracion_externa.y = particle_i.aceleracion_externa.y - (colisiones_de_rigidez*dify) + (amortiguamiento*particle_i.vy);
                particles[i].aceleracion_externa.y =particle_i.aceleracion_externa.y;
                //std::cout <<"adios: "<<contador<<"\n";
            }
        }
        contador++;

        double posicion_z = particle_i.pz + particle_i.hvz*paso_de_tiempo;
        //Límite marcado por z
        if (entry_i.cz == 0) {
            difz = tamaño_de_particula - (posicion_z - zmin);
            if (difz > 10e-10) {
                particle_i.aceleracion_externa.z = particle_i.aceleracion_externa.z + (colisiones_de_rigidez*difz )- (amortiguamiento*particle_i.vz);
                particles[i].aceleracion_externa.z  = particle_i.aceleracion_externa.z;
            }
        }
        else if (entry_i.cz == nz-1) {
            difz = tamaño_de_particula - (zmax - posicion_z);
            if (difz > 10e-10) {
                particle_i.aceleracion_externa.z = particle_i.aceleracion_externa.z - ((colisiones_de_rigidez*difz) + (amortiguamiento*particle_i.vz));
                particles[i].aceleracion_externa.z  = particle_i.aceleracion_externa.z;
            }
        }
    }*/

    /*for (Particle& particle : particles) {
        std::cout <<"aceleracionx: "<<particle.aceleracion_externa.x<<" ";
        std::cout <<"aceleraciony: "<<particle.aceleracion_externa.y<<" ";
        std::cout <<"aceleracionz: "<<particle.aceleracion_externa.z<<"\n ";
    }*/

    //Movimiento de particulas.OK!!!
    /*for (Particle& particle : particles) {
        particle.px = particle.px + (particle.hvx*paso_de_tiempo)+(particle.aceleracion_externa.x*std::pow(paso_de_tiempo, 2));
        particle.py = particle.py + (particle.hvy*paso_de_tiempo)+(particle.aceleracion_externa.y*std::pow(paso_de_tiempo, 2));
        particle.pz = particle.pz + (particle.hvz*paso_de_tiempo)+(particle.aceleracion_externa.z*std::pow(paso_de_tiempo, 2));

        particle.vx = particle.hvx + ((particle.aceleracion_externa.x*paso_de_tiempo)/2);
        particle.vy = particle.hvy + ((particle.aceleracion_externa.y*paso_de_tiempo)/2);
        particle.vz = particle.hvz + ((particle.aceleracion_externa.z*paso_de_tiempo)/2);

        particle.hvx = particle.hvx + (particle.aceleracion_externa.x*paso_de_tiempo);
        particle.hvy = particle.hvy + (particle.aceleracion_externa.y*paso_de_tiempo);
        particle.hvz = particle.hvz + (particle.aceleracion_externa.z*paso_de_tiempo);
    }*/

    /*for (Particle& particle : particles) {
        std::cout <<"aceleracionx: "<<particle.hvx<<" ";
        std::cout <<"aceleraciony: "<<particle.hvy<<" ";
        std::cout <<"aceleracionz: "<<particle.hvz<<"\n ";
    }*/
    /*particleGrid.clear();
    for (Particle& particle : particles) {
        int i = static_cast<int>(std::floor((particle.px - xmin) / sx));
        int j = static_cast<int>(std::floor((particle.py - ymin) / sy));
        int k = static_cast<int>(std::floor((particle.pz - zmin) / sz));


        // Asegurar de que i, j y k estén dentro de los límites.
        i = ajustarLimite(i, 0, nx);
        j = ajustarLimite(j, 0, ny);
        k = ajustarLimite(k, 0, nz);

        // Agrega la partícula a la matriz tridimensional en las coordenadas i, j, k correspondientes.
        particleGrid.push_back({i, j, k, particle});
    }


    //Interacciones con los límites del recinto.okkk.
    //double conntador = 1;
    for (std::vector<ParticleEntry>::size_type i = 0; i < particleGrid.size(); ++i) {
        ParticleEntry &entry_i = particleGrid[i];
        Particle &particle_i = entry_i.particle;
        double d_x = 20;
        //Límite marcado por x
        if (entry_i.cx == 0) {
            d_x = particle_i.px - xmin;
        }
        if (entry_i.cx==nx-1){
            d_x = xmax - particle_i.px;
        }
        if (d_x < 0) {
            if (entry_i.cx == 0) {
                particle_i.px = xmin - d_x;
                particles[i].px = particle_i.px;
                //std::cout <<"adios: "<<contador<<"\n";
            }
            if (entry_i.cx == nx-1) {
                particle_i.px = xmax + d_x;
                particles[i].px = particle_i.px;
                //std::cout <<"adios: "<<contador<<"\n";
            }
            particle_i.vx = -particle_i.vx;
            particles[i].vx = particle_i.vx;
            particle_i.hvx = -particle_i.hvx;
            particles[i].hvx = particle_i.hvx;
        }

        double d_y = 20;
        //Límite marcado por x
        if (entry_i.cy == 0) {
            d_y = particle_i.py - ymin;
        }
        if (entry_i.cy==ny-1){
            d_y = ymax - particle_i.py;
        }
        if (d_y < 0) {
            if (entry_i.cy == 0) {
                particle_i.py = ymin - d_y;
                particles[i].py = particle_i.py;
            }
            if (entry_i.cy == ny-1) {
                particle_i.py = ymax + d_y;
                particles[i].py = particle_i.py;
            }
            particle_i.vy = -particle_i.vy;
            particles[i].vy = particle_i.vy;
            particle_i.hvy = -particle_i.hvy;
            particles[i].hvy = particle_i.hvy;
        }

        double d_z = 20;
        //Límite marcado por x
        if (entry_i.cz == 0) {
            d_z = particle_i.pz - zmin;
        }
        if (entry_i.cz==nz-1){
            d_z = zmax - particle_i.pz;
        }
        if (d_z < 0) {
            if (entry_i.cz == 0) {
                particle_i.pz = zmin - d_z;
                particles[i].pz = particle_i.pz;
            }
            if (entry_i.cz == nz-1) {
                particle_i.pz = zmax + d_z;
                particles[i].pz = particle_i.pz;
            }
            particle_i.vz = -particle_i.vz;
            particles[i].vz = particle_i.vz;
            particle_i.hvz = -particle_i.hvz;
            particles[i].hvz = particle_i.hvz;
        }
    }
    //Reinicio
    for (Particle& particle : particles) {
        particle.rho = 0;
        particle.aceleracion_externa.x = aceleracion_x;
        particle.aceleracion_externa.y = aceleracion_y;
        particle.aceleracion_externa.z = aceleracion_z;
    }

    return 0;
}*/