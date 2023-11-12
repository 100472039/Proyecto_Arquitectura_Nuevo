//
// Created by natha on 30/10/2023.
//
//Versión de Alberto con git checkout

#ifndef PROYECTO_ARQUITECTURA_PROGARGS_HPP
#define PROYECTO_ARQUITECTURA_PROGARGS_HPP

#endif //PROYECTO_ARQUITECTURA_PROGARGS_HPP

#pragma once

#include <vector>
#include <string>
struct VectorAceleracion {
    double x;
    double y;
    double z;

};

struct Particle {
    double px, py, pz;
    double hvx, hvy, hvz;
    double vx, vy, vz;
    double rho;
    VectorAceleracion aceleracion_externa;
    int id;
};

struct LoadResult {
    bool success;
    double ppm;
};

bool loadData(const std::string& filename, std::vector<Particle>& particles, LoadResult& result);
void parametrosSimulacion(double ppm, double& m, double& h);
double smooth;
double masa;

