#include "calculate_functions.h"
#include "grid.h"
#include "initial_conditions.h"
#include "in_out_data.h"

int main()
{
    int checkConfiguration = configuration();
    if (checkConfiguration != 0)
    {
        particles.resize(numParticles);
        predictedParticles.resize(numParticles);

        neighborCount.resize(numParticles);
        sum_rho.resize(numParticles);
        sum_energy.resize(numParticles);
        sum_velocity_x.resize(numParticles);
        sum_velocity_y.resize(numParticles);
        sum_smoothingLength.resize(numParticles);
        force_halo_x.resize(numParticles);
        force_halo_y.resize(numParticles);
    }

    std::string path = "";
    int choice;
    while (true)
    {
        std::cout << "\n1 - start the program with initialization of new values\n2 - launch the program with initialization of values from the file\n";
        std::cin >> choice;
        switch (choice)
        {
        case 1:
            system("CLS");
            initializeGasCloud();
            saveParticlesToFile(t, Dt);
            std::cout << Dt << std::endl;
            SPH();
            return 0;
        case 2:
            system("CLS");
            std::cout << "\nEnter the folder and the file name\n";
            std::cin >> path;
            path = "data\\" + path;
            //FillParticles(path);
            //SPH();
            return 0;
        default:
            system("CLS");
            std::cout << "Something went wrong.";
            return 0;
        }
    }
    return 0;
}

void SPH()
{
    double temp = 0.0;
    int check = 0;
    int check_print = 0;
    int count = 0;
    int count1 = 0;
    while (t <= t_end)
    {
        Predictor();
        dt();
        double A = 4.0 * M_PI * G * density0 * scale_halo * scale_halo;
        for (int i = 0; i < numParticles; ++i)
        {
            if (particles[i].mass > 0)
            {
                double r = particles[i].distanceFromCenter;
                double ksi = computeKsiForHalo(particles[i].x, particles[i].y, 0.0);
                double force_halo_x = computeForceGrav(ksi, particles[i].x, A, 1.0);
                double force_halo_y = computeForceGrav(ksi, particles[i].y, A, 1.0);
                predictedParticles[i].velocityX = particles[i].velocityX + (sum_velocity_x[i] + force_halo_x) * 0.5 * tau;
                predictedParticles[i].velocityY = particles[i].velocityY + (sum_velocity_y[i] + force_halo_y) * 0.5 * tau;
                predictedParticles[i].x = particles[i].x + predictedParticles[i].velocityX * 0.5 * tau;
                predictedParticles[i].y = particles[i].y + predictedParticles[i].velocityY * 0.5 * tau;
                predictedParticles[i].distanceFromCenter = sqrt(particles[i].x * particles[i].x + particles[i].y * particles[i].y);
                if (predictedParticles[i].distanceFromCenter < radius)
                {
                    predictedParticles[i].mass = particles[i].mass;
                    predictedParticles[i].density = particles[i].density + sum_rho[i] * 0.5 * tau;
                    predictedParticles[i].energy = particles[i].energy + sum_energy[i] * 0.5 * tau;
                    predictedParticles[i].pressure = (gamma - 1.0) * predictedParticles[i].energy * predictedParticles[i].density;
                    predictedParticles[i].smoothingLength = particles[i].smoothingLength + sum_smoothingLength[i] * 0.5 * tau;
                }
            }
        }
        Korrector();
        for (int i = 0; i < numParticles; ++i)
        {
            if (predictedParticles[i].mass > 0)
            {
                double r = predictedParticles[i].distanceFromCenter;
                double ksi = computeKsiForHalo(predictedParticles[i].x, predictedParticles[i].y, 0.0);
                double force_halo_x = computeForceGrav(ksi, predictedParticles[i].x, A, 1.0);
                double force_halo_y = computeForceGrav(ksi, predictedParticles[i].y, A, 1.0);
                particles[i].velocityX = predictedParticles[i].velocityX + (sum_velocity_x[i] + force_halo_x) * 0.5 * tau;
                particles[i].velocityY = predictedParticles[i].velocityY + (sum_velocity_y[i] + force_halo_y) * 0.5 * tau;
                particles[i].x = predictedParticles[i].x + particles[i].velocityX * 0.5 * tau;
                particles[i].y = predictedParticles[i].y + particles[i].velocityY * 0.5 * tau;
                particles[i].distanceFromCenter = sqrt(particles[i].x * particles[i].x + particles[i].y * particles[i].y);
                if (particles[i].distanceFromCenter < radius)
                {
                    particles[i].density = predictedParticles[i].density + sum_rho[i] * 0.5 * tau;
                    particles[i].energy = predictedParticles[i].energy + sum_energy[i] * 0.5 * tau;
                    particles[i].pressure = (gamma - 1.0) * particles[i].energy * particles[i].density;
                    particles[i].smoothingLength = predictedParticles[i].smoothingLength + sum_smoothingLength[i] * 0.5 * tau;
                }
            }
        }
        t += tau;
        std::cout << t << " tau: " << tau << std::endl;
        if (t >= Dt)
        {
            if (check_print == 0)
            {
                check_print = 1;
            }
            else
            {
                saveParticlesToFile(t, Dt);
            }
            std::cout << t << " DT: " << Dt << std::endl;
            Dt += shag_dt;
        }
    }
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
