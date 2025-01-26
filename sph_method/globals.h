#pragma once
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
double shag_dt = 0.0001;
double dt_out = 10000;

struct Particle
{
    double x, y;
    double pressure{ 0.0 };
    double density{ 0.0 };
    double mass{ 0.0 };
    double velocityX{ 0.0 };
    double velocityY{ 0.0 };
    double smoothingLength{ 0.0 };
    double energy{ 0.0 };
    double distanceFromCenter{ 0 };
    int neighbors = { 0 };
};

std::vector<Particle> particles;
std::vector<Particle> predictedParticles;

std::vector <int> neighborCount;
std::vector <double> sum_rho;
std::vector <double> sum_energy;
std::vector <double> sum_velocity_x;
std::vector <double> sum_velocity_y;
std::vector <double> sum_smoothingLength;
std::vector <double> force_halo_x;
std::vector <double> force_halo_y;
