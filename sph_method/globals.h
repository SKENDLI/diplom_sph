#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <codecvt>
#include <locale>
#include <cmath>
#include <Windows.h>
#include <random>
#include <sstream>
#include <numeric>
#include <algorithm>
#include <omp.h>

const double M_PI = 3.14159265358979323846;
const double G = 6.67430e-11;

int numParticles;
double radius;
double density0;
double pressure;
double tau;
double gamma;
double alpha;
double beta;
double eps;
double maximum;
double t_end;
double t;
double Dt;
double system_Mass;
double scale;
double mass_source;
double a_halo;
double b_halo;
double c_halo;
double massa_halo;
double scale_halo;
double shag_dt = 0.01;
double dt_out = 100;
double B_rho = 0.0;
double Mah_d;
double AMS;
double r_core;
double rmax_hl;
double A;

struct Particle {
    double x = 0.0;
    double y = 0.0;
    double velocityX = 0.0;
    double velocityY = 0.0;
    double density = 0.0;
    double pressure = 0.0;
    double energy = 0.0;
    double mass = 0.0;
    double smoothingLength = 0.0;
    double distanceFromCenter = 0.0;
    int neighbors = 0;
};

std::vector<Particle> particles;
std::vector<Particle> predictedParticles;
std::vector<Particle> initialParticles;

std::vector <int> neighborCount;
std::vector <double> sum_rho;
std::vector <double> sum_energy;
std::vector <double> sum_velocity_x;
std::vector <double> sum_velocity_y;
std::vector <double> sum_smoothingLength;
std::vector <double> force_halo_x;
std::vector <double> force_halo_y;

void saveParticlesToFile(double times, double Dt);
void initializeGasCloud();
void Predictor(double dt);
void Korrector(double dt);
void SPH();
void dt();
int configuration();
double computeKsiForHalo(double x, double y, double z);
double computeForceGrav(double ksi, double koord, double A, double halo_OXYZ);
double computeRho(double dist, double rad);
double dComputeRho(double dist, double rad);
double fun_mass(double r1, double r2, int Nmm, double pi2, double Lf);
double SoundSpeed(double pressure, double density);
double Viscosity(
    double x_i, double y_i,
    double x_j, double y_j,
    double velocity_ij_x, double velocity_ij_y,
    double h_ij,
    double rho_i, double rho_j,
    double pressure_i, double pressure_j,
    double r_ij
);
double dW(double r_ij, double h, double r);