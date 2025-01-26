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


class Grid
{
public:
    double cellSize;
    int gridSizeX, gridSizeY;
    double a, b;
    std::vector<std::vector<std::vector<Particle*>>> cells;

    Grid(double maxSmoothingLength, const std::vector<Particle>& particles, double minX, double maxX)
        : a(minX), b(maxX)
    {
        cellSize = 2 * maxSmoothingLength;

        gridSizeX = static_cast<int>(ceil((b - a) / cellSize));
        gridSizeY = static_cast<int>(ceil((b - a) / cellSize));

        cells.resize(gridSizeX);
        for (int i = 0; i < gridSizeX; ++i)
        {
            cells[i].resize(gridSizeY);
        }

        for (const auto& particle : particles)
        {
            if (particle.x >= a && particle.x <= b && particle.y >= a && particle.y <= b)
            {
                int cellX = static_cast<int>((particle.x - a) / cellSize);
                int cellY = static_cast<int>((particle.y - a) / cellSize);

                cellX = min(max(cellX, 0), gridSizeX - 1);
                cellY = min(max(cellY, 0), gridSizeY - 1);

                cells[cellX][cellY].push_back(const_cast<Particle*>(&particle));
            }
        }
    }

    ~Grid()
    {
        for (int i = 0; i < gridSizeX; ++i)
        {
            for (int j = 0; j < gridSizeY; ++j)
            {
                cells[i][j].clear();
            }
        }
        cells.clear();
    }

    const std::vector<Particle*>& getParticlesInCell(int cellX, int cellY) const
    {
        return cells[cellX][cellY];
    }

    int getGridSizeX() const
    {
        return gridSizeX;
    }

    int getGridSizeY() const
    {
        return gridSizeY;
    }
};


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

double computeKsiForHalo(double x, double y, double z)
{
    return sqrt((x / a_halo) * (x / a_halo) + (y / b_halo) * (y / b_halo) + (z / c_halo) * (z / c_halo));
}

double computeForceGrav(double ksi, double koord, double A, double halo_OXYZ)
{
    return -A * (1.0 / ksi - 1.0 / (ksi * ksi) * atan(ksi)) * (koord / (ksi * halo_OXYZ));
}

double computeForcePressure(double p, double density_disk, double r, double scale_radius)
{
    return -gamma * p / density_disk * pow(density_disk, gamma - 2.0) * (-2.0 * sinh(r / scale_radius) / scale_radius);
}

double computeRho(double dist, double rad)
{
    double rho = (1.0 / (cosh(dist / rad) * cosh(dist / rad)));
    return rho;
}
