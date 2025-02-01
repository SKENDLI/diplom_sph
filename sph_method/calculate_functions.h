#pragma once
#include "globals.h"
#include "grid.h"
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
    double cosh_val = cosh(dist / rad);
    return 1.0 / (cosh_val * cosh_val);
}


void Predictor()
{
    std::vector<Particle*> neighbors;

    // Создаём сетку
    double minX = particles[0].x;
    double minY = particles[0].y;
    double maxX = minX;
    double maxY = minY;

    for (const auto& particle : particles)
    {
        minX = min(minX, particle.x);
        minY = min(minY, particle.y);
        maxX = max(maxX, particle.x);
        maxY = max(maxY, particle.y);
    }

    // Создаём объект Grid с нужными параметрами
    Grid grid(particles[0].smoothingLength, particles, minX, maxX);

    // Шаг 3: Для каждой частицы вычисляем необходимые величины
    for (int i = 0; i < numParticles; i++)
    {
        Particle& particle = particles[i];
        sum_rho[i] = 0.0;
        sum_energy[i] = 0.0;
        sum_velocity_x[i] = 0.0;
        sum_velocity_y[i] = 0.0;
        sum_smoothingLength[i] = 0.0;
        neighborCount[i] = 0;

        if (particle.distanceFromCenter < radius)
        {
            neighbors.clear();

            // Получаем соседей из соседних ячеек
            int cellX = static_cast<int>((particle.x - minX) / grid.cellSize);
            int cellY = static_cast<int>((particle.y - minY) / grid.cellSize);

            // Проверяем соседние ячейки
            for (int dx = -1; dx <= 1; dx++)
            {
                for (int dy = -1; dy <= 1; dy++)
                {
                    int neighborCellX = cellX + dx;
                    int neighborCellY = cellY + dy;

                    if (neighborCellX >= 0 && neighborCellX < grid.getGridSizeX() &&
                        neighborCellY >= 0 && neighborCellY < grid.getGridSizeY())
                    {
                        // Добавляем соседей из соседней ячейки
                        const auto& particlesInCell = grid.getParticlesInCell(neighborCellX, neighborCellY);
                        neighbors.insert(neighbors.end(), particlesInCell.begin(), particlesInCell.end());
                    }
                }
            }

            // Обрабатываем соседей
            for (Particle* neighbor : neighbors)
            {
                double h = (particle.smoothingLength + neighbor->smoothingLength) * 0.5;
                double r_x = particle.x - neighbor->x;
                double r_y = particle.y - neighbor->y;
                double r = std::sqrt(r_x * r_x + r_y * r_y);

                if ((r / h) <= 2.0 && (r / h) > 0.0)
                {
                    double velocity_ij_x = particle.velocityX - neighbor->velocityX;
                    double velocity_ij_y = particle.velocityY - neighbor->velocityY;
                    double dKernelX = dW(r_x, h, r);
                    double dKernelY = dW(r_y, h, r);
                    double dKernelXY = dKernelX * velocity_ij_x + dKernelY * velocity_ij_y;

                    double Pressure_rho_ij = neighbor->pressure / (neighbor->density * neighbor->density)
                        + particle.pressure / (particle.density * particle.density)
                        + Viscosity(particle.x, particle.y, neighbor->x, neighbor->y,
                            velocity_ij_x, velocity_ij_y, h, particle.density,
                            neighbor->density, particle.pressure, neighbor->pressure, r);

                    sum_rho[i] += particle.mass * dKernelXY;
                    if (particle.distanceFromCenter < radius)
                    {
                        sum_velocity_x[i] -= particle.mass * Pressure_rho_ij * dKernelX;
                        sum_velocity_y[i] -= particle.mass * Pressure_rho_ij * dKernelY;
                    }
                    sum_energy[i] += particle.mass * Pressure_rho_ij * dKernelXY * 0.5;

                }
                neighborCount[i] += 1;
            }

            double dDensity_dt = sum_rho[i] / particle.density;

            double nu = 2.0;
            sum_smoothingLength[i] = -particle.smoothingLength * (1.0 / (nu * particle.density)) * dDensity_dt;
            particles[i].neighbors = neighborCount[i];
        }
    }
}


void Korrector()
{
    std::vector<Particle*> neighbors;

    double minX = predictedParticles[0].x;
    double minY = predictedParticles[0].y;
    double maxX = minX;
    double maxY = minY;

    for (const auto& particle : predictedParticles)
    {
        minX = min(minX, particle.x);
        minY = min(minY, particle.y);
        maxX = max(maxX, particle.x);
        maxY = max(maxY, particle.y);
    }

    Grid grid(predictedParticles[0].smoothingLength, predictedParticles, minX, maxX);

    for (int i = 0; i < numParticles; i++)
    {
        Particle& particle = predictedParticles[i];
        sum_rho[i] = 0.0;
        sum_energy[i] = 0.0;
        sum_velocity_x[i] = 0.0;
        sum_velocity_y[i] = 0.0;
        sum_smoothingLength[i] = 0.0;

        if (particle.distanceFromCenter < radius)
        {
            neighbors.clear();

            int cellX = static_cast<int>((particle.x - minX) / grid.cellSize);
            int cellY = static_cast<int>((particle.y - minY) / grid.cellSize);

            for (int dx = -1; dx <= 1; dx++)
            {
                for (int dy = -1; dy <= 1; dy++)
                {
                    int neighborCellX = cellX + dx;
                    int neighborCellY = cellY + dy;

                    if (neighborCellX >= 0 && neighborCellX < grid.getGridSizeX() &&
                        neighborCellY >= 0 && neighborCellY < grid.getGridSizeY())
                    {
                        const auto& particlesInCell = grid.getParticlesInCell(neighborCellX, neighborCellY);
                        neighbors.insert(neighbors.end(), particlesInCell.begin(), particlesInCell.end());
                    }
                }
            }

            for (Particle* neighbor : neighbors)
            {
                double h = (particle.smoothingLength + neighbor->smoothingLength) * 0.5;
                double r_x = particle.x - neighbor->x;
                double r_y = particle.y - neighbor->y;
                double r = std::sqrt(r_x * r_x + r_y * r_y);

                if ((r / h) <= 2.0 && (r / h) > 0.0)
                {
                    double velocity_ij_x = particle.velocityX - neighbor->velocityX;
                    double velocity_ij_y = particle.velocityY - neighbor->velocityY;
                    double dKernelX = dW(r_x, h, r);
                    double dKernelY = dW(r_y, h, r);
                    double dKernelXY = dKernelX * velocity_ij_x + dKernelY * velocity_ij_y;

                    double Pressure_rho_ij = neighbor->pressure / (neighbor->density * neighbor->density)
                        + particle.pressure / (particle.density * particle.density)
                        + Viscosity(particle.x, particle.y, neighbor->x, neighbor->y,
                            velocity_ij_x, velocity_ij_y, h, particle.density,
                            neighbor->density, particle.pressure, neighbor->pressure, r);

                    sum_rho[i] += particle.mass * dKernelXY;
                    if (particle.distanceFromCenter < radius)
                    {
                        sum_velocity_x[i] -= particle.mass * Pressure_rho_ij * dKernelX;
                        sum_velocity_y[i] -= particle.mass * Pressure_rho_ij * dKernelY;
                    }
                    sum_energy[i] += particle.mass * Pressure_rho_ij * dKernelXY * 0.5;
                }
            }

            double dDensity_dt = sum_rho[i] / particle.density;

            double nu = 2.0;

            sum_smoothingLength[i] = -particle.smoothingLength * (1.0 / (nu * particle.density)) * dDensity_dt;
        }
    }
}

double dW(double r_ij, double h, double r)
{
    double z = 0.0;
    double q = r / h;
    double sigma = 10.0 / (7.0 * M_PI);

    if ((q > 0.0) && (q <= 1.0))
    {
        z = (-3.0 + 2.25 * q) * q;
    }
    else if ((q > 1.0) && (q <= 2.0))
    {
        double q1 = 2.0 - q;
        z = -0.75 * q1 * q1;
    }

    if (q > 2.0)
    {
        z = 0.0;
    }

    z = z * sigma * r_ij / (r * h * h * h);
    return z;
}

double Viscosity(double x_i, double y_i, double x_j, double y_j, double velocity_ij_x, double velocity_ij_y, double h_ij, double rho_i, double rho_j, double pressure_i, double pressure_j, double r_ij)
{
    double x_ij = x_i - x_j;
    double y_ij = y_i - y_j;
    double c_i = SoundSpeed(pressure_i, rho_i);  // Звуковая скорость для частицы i
    double c_j = SoundSpeed(pressure_j, rho_j);  // Звуковая скорость для частицы j
    double rho_ij = (rho_i + rho_j) * 0.5;      // Средняя плотность
    double c_ij = (c_i + c_j) * 0.5;            // Средняя звуковая скорость
    // Вычисление вязкости (mu)
    double mu = (h_ij * (velocity_ij_x * x_ij + velocity_ij_y * y_ij)) / (r_ij * r_ij + eps * eps * h_ij * h_ij);

    // Инициализация результата
    double result = 0.0;

    // Если частицы движутся друг относительно друга (отрицательная компонента скорости)
    if ((velocity_ij_x * x_ij + velocity_ij_y * y_ij) < 0.0)
    {
        // Формула для вязкостного термина, с учетом коэффициентов alpha и beta
        result = (-(mu)*alpha * c_ij + beta * mu * mu) / rho_ij;
    }

    // Обновление максимальной вязкости (если необходимо)
    if (mu > maximum)
    {
        maximum = mu;
    }

    return result;
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

        dt_cv = 0.24 * 0.4 * particles[i].smoothingLength / (C + 0.6 * (alpha * C + beta * maximum));
        dt_f = 0.25 * 0.25 * sqrt(particles[i].smoothingLength / velocity_magnitude);

        dt = min(dt_cv, dt_f);
        tau = min(tau, dt);
    }
}

double computeTotalEnergy(const std::vector<Particle>& particles, double G) {
    double totalEnergy = 0.0;

    for (int i = 0; i < particles.size(); ++i) {
        Particle pi = particles[i];
        // Кинетическая энергия
        double kineticEnergy = 0.5 * pi.mass * (pi.velocityX * pi.velocityX + pi.velocityY * pi.velocityY);

        // Внутренняя энергия
        double internalEnergy = pi.energy * pi.mass;

        // Потенциальная энергия (гравитационная)
        double potentialEnergy = 0.0;
        for (int j = 0; j < particles.size(); ++j) {
            if (i == j) continue;
            Particle pj = particles[j];
            double dx = pi.x - pj.x;
            double dy = pi.y - pj.y;
            double distance = sqrt(dx * dx + dy * dy);
            potentialEnergy -= G * pi.mass * pj.mass / distance;
        }

        // Суммируем все виды энергии
        totalEnergy += kineticEnergy + internalEnergy + potentialEnergy;
    }

    return totalEnergy;
}

void checkEnergyConservation(const std::vector<Particle>& particlesBefore, const std::vector<Particle>& particlesAfter, double G, double tolerance = 1e-6) 
{
    double energyBefore = computeTotalEnergy(particlesBefore, G);
    double energyAfter = computeTotalEnergy(particlesAfter, G);

    double energyDifference = std::abs(energyAfter - energyBefore);
    std::cout << "ENERGY CHECK: " << energyDifference << std::endl;
}
