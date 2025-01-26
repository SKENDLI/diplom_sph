#include "globals.h"

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
