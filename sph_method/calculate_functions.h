#pragma once
#include "globals.h"
#include "grid.h"
double computeKsiForHalo(double x, double y, double z)
{
    return sqrt((x / a_halo) * (x / a_halo) + (y / b_halo) * (y / b_halo) + (z / c_halo) * (z / c_halo));
}

double computeForceGrav(double ksi, double koord, double A, double halo_OXYZ) {
    return -A * (1.0 / ksi - 1.0 / (ksi * ksi) * atan(ksi)) * (koord / (ksi * halo_OXYZ));
}

double computeRho(double distance, double scale_radius)
{
    double cosh_val = cosh(distance / scale_radius);
    return 1.0 / (cosh_val);
}

double dComputeRho(double distance, double scale_radius)
{
    return 1.0 / scale_radius * sinh(distance / scale_radius) / pow(cosh(distance / scale_radius), 2.0);
}

double fun_mass(double r1, double r2, int Nmm, double pi2, double Lf) {
    double dr = (r2 - r1) / Nmm;
    double fun_mass = 0.0;
    for (int k = 1; k <= Nmm; ++k) {
        double rr = r1 + dr * (k - 1);
        fun_mass += 0.5 * (rr * computeRho(rr, Lf) + (rr + dr) * computeRho(rr+dr, Lf)) * dr;
    }
    return pi2 * fun_mass;
}

void ComputeForces(const std::vector<Particle>& particles,
    std::vector<double>& sum_rho,
    std::vector<double>& sum_vx,
    std::vector<double>& sum_vy,
    std::vector<double>& sum_energy,
    double radius_limit)
{
    double max_h = 0.0;
    for (const auto& p : particles) {
        if (p.distanceFromCenter < radius_limit) {
            max_h = max(max_h, p.smoothingLength);
        }
    }

    GasDiskGrid grid(radius_limit, max_h);
    grid.rebuild(particles);
#pragma omp parallel for
    for (int i = 0; i < particles.size(); ++i) {
        const Particle& pi = particles[i];
        if (pi.distanceFromCenter >= radius_limit) continue;
        std::vector<Particle*> neighbors;
        grid.get_neighbors(&pi, neighbors);

        double sum_rho_i = 0.0;
        double sum_vx_i = 0.0;
        double sum_vy_i = 0.0;
        double sum_energy_i = 0.0;

        for (Particle* pj : neighbors) {
            if (&pi == pj) continue;

            const double h_avg = 0.5 * (pi.smoothingLength + pj->smoothingLength);
            const double dx = pi.x - pj->x;
            const double dy = pi.y - pj->y;
            const double r = std::sqrt(dx * dx + dy * dy);

            if (r > 2.0 * h_avg) continue;

            const double dWx = dW(dx, h_avg, r);
            const double dWy = dW(dy, h_avg, r);

            // Производная плотности
            sum_rho_i += pj->mass * (
                (pj->velocityX - pi.velocityX) * dWx +
                (pj->velocityY - pi.velocityY) * dWy);

            // Сила давления
            const double P_term = (pi.pressure / (pi.density * pi.density) +
                pj->pressure / (pj->density * pj->density));
            sum_vx_i += -pj->mass * P_term * dWx;
            sum_vy_i += -pj->mass * P_term * dWy;

            // Вязкость
            const double visc = Viscosity(
                pi.x, pi.y,
                pj->x, pj->y,
                pi.velocityX - pj->velocityX,
                pi.velocityY - pj->velocityY,
                h_avg,
                pi.density, pj->density,
                pi.pressure, pj->pressure,
                r
            );
            sum_vx_i += -pj->mass * visc * dWx;
            sum_vy_i += -pj->mass * visc * dWy;

            // Энергия
            sum_energy_i += 0.5 * pj->mass * (P_term + visc) * (
                (pi.velocityX - pj->velocityX) * dWx +
                (pi.velocityY - pj->velocityY) * dWy
                );
        }

        sum_rho[i] = sum_rho_i;
        sum_vx[i] = sum_vx_i;
        sum_vy[i] = sum_vy_i;
        sum_energy[i] = sum_energy_i;
    }
}

void Predictor(double dt) {
    for (int i = 0; i < numParticles; ++i) {
        predictedParticles[i] = particles[i];

        predictedParticles[i].velocityX += 0.5 * dt * sum_velocity_x[i];
        predictedParticles[i].velocityY += 0.5 * dt * sum_velocity_y[i];


        predictedParticles[i].x += predictedParticles[i].velocityX * dt;
        predictedParticles[i].y += predictedParticles[i].velocityY * dt;

        predictedParticles[i].density += sum_rho[i] * dt;
        predictedParticles[i].pressure = (gamma - 1.0) *
            predictedParticles[i].energy * predictedParticles[i].density;
    }
}

// Корректор (окончательное интегрирование)
void Korrector(double dt) {
    for (int i = 0; i < numParticles; ++i) {
        particles[i].velocityX += sum_velocity_x[i] * dt;
        particles[i].velocityY += sum_velocity_y[i] * dt;

        particles[i].x = predictedParticles[i].x +
            0.5 * dt * (particles[i].velocityX + predictedParticles[i].velocityX);
        particles[i].y = predictedParticles[i].y +
            0.5 * dt * (particles[i].velocityY + predictedParticles[i].velocityY);

        particles[i].density = predictedParticles[i].density + sum_rho[i] * dt;
        particles[i].pressure = (gamma - 1.0) * particles[i].energy * particles[i].density;
    }
}

double dW(double r_ij, double h, double r) {
    const double sigma = 15.0 / (7.0 * M_PI * h * h * h); // Правильная нормировка
    double q = r / h;
    double z = 0.0;

    if (q <= 1.0) {
        z = (-3.0 + 2.25 * q) * q;
    }
    else if (q <= 2.0) {
        z = -0.75 * pow(2.0 - q, 2);
    }

    return z * sigma * r_ij / r; // Градиент по компонентам
}

double Viscosity(
    double x_i, double y_i,
    double x_j, double y_j,
    double velocity_ij_x, double velocity_ij_y,
    double h_ij,
    double rho_i, double rho_j,
    double pressure_i, double pressure_j,
    double r_ij
) {
    double dx = x_i - x_j;
    double dy = y_i - y_j;
    double r = std::sqrt(dx * dx + dy * dy);
    if (r == 0.0) return 0.0;

    double v_rel = velocity_ij_x * dx + velocity_ij_y * dy;
    if (v_rel > 0.0) return 0.0;

    double rho_ij = (rho_i + rho_j) * 0.5;
    double c_ij = SoundSpeed(pressure_i, rho_i) + SoundSpeed(pressure_j, rho_j);
    double mu = -v_rel / (r + 1e-5);
    if (mu > maximum) {
        maximum = mu;
    }
    return (alpha * c_ij * mu + beta * mu * mu) / rho_ij;
}

double SoundSpeed(double p_i, double rho_i)
{
    return sqrt(gamma * p_i / rho_i);
}

void dt()
{
    double dt_cv, dt_f, dt;
    for (int i = 0; i < numParticles; i++)
    {
        if (particles[i].distanceFromCenter > radius)
        {
            continue;
        }

        double C = SoundSpeed(particles[i].pressure, particles[i].density);
        double velocity_magnitude = sqrt(particles[i].velocityX * particles[i].velocityX + particles[i].velocityY * particles[i].velocityY);

        dt_cv = 0.3 * particles[i].smoothingLength / (C + 0.6 * (alpha * C + beta * maximum));
        dt_f = 0.3 * sqrt(particles[i].smoothingLength / velocity_magnitude);

        dt = min(dt_cv, dt_f);
        tau = min(tau, dt);
    }
}

inline double particleDistance(const Particle& a, const Particle& b, double softening = 1e-4) {
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    return sqrt(dx * dx + dy * dy + softening * softening);
}

// Вычисление полной энергии системы
double computeTotalEnergy(const std::vector<Particle>& particles, double G, double softening = 1e-4) {
    double totalEnergy = 0.0;
    const size_t n = particles.size();

   #pragma omp parallel for reduction(+:totalEnergy)
    for (size_t i = 0; i < n; ++i) {
        const Particle& p = particles[i];
        const double v_sq = p.velocityX * p.velocityX + p.velocityY * p.velocityY;
        totalEnergy += p.mass * (0.5 * v_sq + p.energy);    }

   #pragma omp parallel for reduction(-:totalEnergy)
    for (size_t i = 0; i < n; ++i) {
        const Particle& pi = particles[i];
        for (size_t j = i + 1; j < n; ++j) {
            const Particle& pj = particles[j];
            const double r = particleDistance(pi, pj, softening);
            totalEnergy -= G * pi.mass * pj.mass / r;
        }
    }

    return totalEnergy;
}

// Проверка сохранения энергии
void checkEnergyConservation(
    std::vector<Particle>& particlesBefore,
    std::vector<Particle>& particlesAfter,
    double G,
    double tolerance = 1e-3,
    double softening = 1e-4
) {
    const double energyBefore = computeTotalEnergy(particlesBefore, G, softening);
    const double energyAfter = computeTotalEnergy(particlesAfter, G, softening);
    const double abs_diff = abs(energyAfter - energyBefore);
    const double rel_diff = abs_diff / max(abs(energyBefore), abs(energyAfter));

    std::cout << "\n=== Energy Conservation Check ==="
        << "\nInitial energy: " << energyBefore
        << "\nFinal energy:   " << energyAfter
        << "\nAbsolute difference: " << abs_diff
        << "\nRelative difference: " << rel_diff * 100 << "%"
        << "\nTolerance: " << tolerance * 100 << "%"
        << std::endl;
}