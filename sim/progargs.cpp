#include <iostream>
#include <fstream>
#include <vector>
#include "progargs.hpp"
#include "../constantes.hpp"
#include "block.cpp"
#include "bloque.cpp"

int ReadFile(std::vector<Particle> &particles, std::ifstream &file, int np, int contador);
template <typename T>
requires(std::is_integral_v<T> or std::is_floating_point_v<T>)
char * as_writable_buffer(T & value) {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
    return reinterpret_cast<char *>(&value);
}

template <typename T>
requires(std::is_integral_v<T> or std::is_floating_point_v<T>)
char const * as_buffer(T const & value) {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
    return reinterpret_cast<char const *>(&value);
}

template <typename T>
requires(std::is_integral_v<T> or std::is_floating_point_v<T>)
T read_binary_value(std::istream & is) {
    T value{};
    is.read(as_writable_buffer(value), sizeof(value));
    return value;
}


bool loadData(const std::string& filename, std::vector<Particle>& particles, float& ppm, int& np) {
    // Abre el archivo binario.
    std::ifstream file(filename, std::ios::binary);

    if (!file.is_open()) {
        std::cerr << "Error al abrir el archivo de entrada.\n";
        return false;
    }

    ppm = read_binary_value<float>(file);
    np = read_binary_value<int>(file);
    std::cout << "ppm: " << ppm << "\n";;
    std::cout << "np: " << np << "\n";;

    if (np == 0) {
        std::cerr << "Invalid number of particles: " << np << "\n";;
        file.close();
        return false;
    }

    int contador = ReadFile(particles, file, np, 0);

    if (contador != np) {
        std::cerr << "Error: Number of particles mismatch. Header:" << np << ", Found:" << contador << "\n";;
        file.close();
        return false;
    }
    file.close();
    return true;
}

int ReadFile(std::vector<Particle>& particles, std::ifstream& file, int np, int contador) {
    for (int i = 0; i < np; ++i) {
        Particle particle;

        // Lee los valores en formato float
        float px_float, py_float, pz_float, hvx_float, hvy_float, hvz_float, vx_float, vy_float, vz_float;

        px_float = read_binary_value<float>(file);
        py_float = read_binary_value<float>(file);
        pz_float = read_binary_value<float>(file);
        hvx_float = read_binary_value<float>(file);
        hvy_float = read_binary_value<float>(file);
        hvz_float = read_binary_value<float>(file);
        vx_float = read_binary_value<float>(file);
        vy_float = read_binary_value<float>(file);
        vz_float = read_binary_value<float>(file);


        // Convierte los valores float a double y los guarda en la estructura Particle
        particle.px = static_cast<double>(px_float);
        particle.py = static_cast<double>(py_float);
        particle.pz = static_cast<double>(pz_float);
        particle.hvx = static_cast<double>(hvx_float);
        particle.hvy = static_cast<double>(hvy_float);
        particle.hvz = static_cast<double>(hvz_float);
        particle.vx = static_cast<double>(vx_float);
        particle.vy = static_cast<double>(vy_float);
        particle.vz = static_cast<double>(vz_float);
        //Inicializacion de particula
        particle.rho = 0.0;
        particle.aceleracion_externa.x = aceleracion_x;
        particle.aceleracion_externa.y = aceleracion_y;
        particle.aceleracion_externa.z = aceleracion_z;
        particle.id = contador;

        ++contador;
        particles.push_back(particle);
    }

    return contador;
}

void parametrosSimulacion(double ppm, double& m, double& h) {
    m = densidad / (ppm * ppm * ppm);
    h = radio / ppm;
    std::cout << "m: " << m << "\n";
    std::cout << "h: " << h << "\n";
    ::smooth=h;
    ::masa=m;
}

using namespace std;

bool isNumeric(const string& str) {
    for (char c : str) {
        if (!isdigit(c)) {
            return false;
        }
    }
    return true;
}

void escribir_datos(const std::string& filename, float ppm, int np, std::string particulas) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error.\n";
        return;
    }

    file << "ppm: " << ppm << "\n";
    file << "np: " << np << "\n";
    file << "particulas:\n\n" << particulas << std::endl;

    file.close();
    std::cout << "Datos almacenados correctamente" << "\n";
}

int sim_main(int argc, char* argv[]) {
    if (argc == 1) {
        cerr << "Error: Invalid number of arguments: 0" << "\n";
        return -1;
    }
    if (argc == 2) {
        cerr << "Error: Invalid number of arguments: 1" << "\n";
        return -1;
    }
    if (argc == 3) {
        cerr << "Error: Invalid number of arguments: 2" << "\n";
        return -1;
    }
    if (argc > 4) {
        cerr << "Error: Invalid number of arguments: 4" << "\n";
        return -1;
    }

    // Obtener el número de pasos de tiempo desde el primer argumento
    if (!isNumeric(argv[1])) {
        cerr << "Error: Time steps must be numeric." << "\n";
        return -1;
    }

    int num_steps = atoi(argv[1]);
    // Verificar si el número de pasos de tiempo es válido
    if (num_steps <= 0) {
        cerr << "Error: Invalid number of time steps." << "\n";
        return -2;
    }

    const char* input_file = argv[2];
    const char* output_file = argv[3];

    // Abrir el archivo de entrada
    ifstream input(input_file, ios::binary);
    if (!input) {
        cerr << "Error: Cannot open " << input_file << " for reading." << "\n";
        return -3;
    }

    // Crear y abrir el archivo de salida (si no existe)
    ofstream output(output_file, ios::binary | ios::app);
    if (!output) {
        cerr << "Error: Cannot open " << output_file << " for writing." <<"\n";
        return -4;
    }

    std::vector<Particle> particles;
    float ppm;
    int np;

    if (!loadData(argv[2], particles, ppm, np)) {
        std::cerr << "Error al cargar el archivo de entrada.\n";
        return 1;
    }

    double m;
    double h;
    parametrosSimulacion(ppm, m, h);
    for (int i=0; i<1; i++) {
        generar_malla(particles);
    }
/*
 *
    //Intento de escritura de particulas en el archivo final
    Grid grid;
    std::vector<double> n = grid.getN();
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

    std::string particulas_string;
    int cont = 1;
    for (std::vector<ParticleEntry>::size_type i = 0; i < particleGrid.size(); ++i) {
        ParticleEntry &entry_i = particleGrid[i];
        Particle &particle_i = entry_i.particle;
        particulas_string = particulas_string + "partícula " + std::to_string(cont) + ":\n"
                            + "\tp: " + std::to_string(particle_i.px) + "," + std::to_string(particle_i.py) + "," + std::to_string(particle_i.pz) + "\n"
                            + "\thv: " + std::to_string(particle_i.hvx) + "," + std::to_string(particle_i.hvy) + "," + std::to_string(particle_i.hvz) + "\n"
                            + "\tv: " + std::to_string(particle_i.vx) + "," + std::to_string(particle_i.vy) + "," + std::to_string(particle_i.vz) + "\n";
        cont++;
    }

    // Cerrar los archivos
    input.close();
    output.close();

    //Escritura de datos
    escribir_datos("./final.fld", ppm, np, particulas_string);*/

    return 0;
}