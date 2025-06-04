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
#include <mutex>
#include <unordered_map>
#include <functional>
#include <limits>

using namespace std;
std::mutex logMutex;

const double M_PI = 3.14159265358979323846;
const double G = 6.67430e-11;

int numParticles;
double radius;
double tau;
double gamma;
double alpha;
double beta;
double eps;
double maximum = -INFINITY;
double t_end;
double t;
double Dt;
double system_Mass;
double scale;
double scale_halo;
double shag_dt = 0.1;
double dt_out = 10;
double B_rho = 0.0;
double Mah_d;
double AMS;
double r_core;
double rmax_hl;
double A;
double Ap;
double Lh;
double x_in, x_ex, y_in, y_ex;
double hp_max;
double h;
double h2;
int Ngx;
int Ngy;
double dtn;
double courant = 0.2;
double Omega_h;
double a_h;
double A_G;
double tau_h;
double initial_energy;

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

std::vector<int> neighborCount;
std::vector<double> Rho;
std::vector<double> Rhot;
std::vector<double> Fu_p;
std::vector<double> Fv_p;
std::vector<double> Fe_p;
std::vector<double> ht_p;
std::vector<double> force_halo_x;
std::vector<double> force_halo_y;

double max_x = -INFINITY;
double max_y = -INFINITY;
double check_max_mu = 9999999;

void saveParticlesToFile(double times, double Dt);
void initializeGasCloud();
void SPH();
int configuration();
double computeRho(double dist, double rad);
double dComputeRho(double dist, double rad);
double fun_mass(double r1, double r2, int Nmm, double pi2, double Lf);
double SoundSpeed(double pressure, double density);
double W(double s, double hij);
double gradW(double s, double rij, double hij);
double mu(double s, double urdot, double hij);
double find_radius(double rrk, double rrk_1, double target_mass, double B_rho, double m_p, double Lf, double pi2);
double mass_difference(double rr, double rrk, double rrk_1, double target_mass, double B_rho, double m_p, double Lf, double pi2);
template <typename T>
void logVariable(const T& var);
template <typename T>
void logVariable(const std::vector<T>& vec);
void computeForcesGravCool(std::vector<Particle>& particles, std::vector<double>& force_halo_x, std::vector<double>& force_halo_y);
void law_cons(double& Energy, const std::vector<Particle>& parts);
void save_laws_cons(double t, double& next_print_time, double initial_energy, vector<Particle>& particles);