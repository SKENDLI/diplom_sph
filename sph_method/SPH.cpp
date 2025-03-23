#include "calculate_functions.h"
#include "initial_conditions.h"
#include "in_out_data_functions.h"

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
        initialParticles.resize(numParticles);
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
            initialParticles = particles;
            //checkEnergyConservation(initialParticles, particles, G);
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

void SPH() {
    // Временные массивы для хранения ускорений и производных
    std::vector<double> sum_velocity_x_old(numParticles);
    std::vector<double> sum_velocity_y_old(numParticles);
    std::vector<double> f_halo_x_old(numParticles);
    std::vector<double> f_halo_y_old(numParticles);
    std::vector<double> sum_rho_old(numParticles);
    std::vector<double> sum_energy_old(numParticles);

    int check_print = 0;
    while (t <= t_end) {
        // Шаг 1: Вычисляем силы для текущего состояния
        ComputeForces(particles, sum_rho, sum_velocity_x, sum_velocity_y, sum_energy, radius);

        // Сохраняем текущие ускорения и производные
        sum_velocity_x_old = sum_velocity_x;
        sum_velocity_y_old = sum_velocity_y;
        sum_rho_old = sum_rho;
        sum_energy_old = sum_energy;

        // Вычисляем гравитацию гало для текущих частиц
        for (int i = 0; i < numParticles; ++i) {
            double ksi = computeKsiForHalo(particles[i].x, particles[i].y, 0.0);
            f_halo_x_old[i] = computeForceGrav(ksi, particles[i].x, A, a_halo);
            f_halo_y_old[i] = computeForceGrav(ksi, particles[i].y, A, b_halo);
        }

        // Обновляем временной шаг
        dt();

        // Шаг 2: Предиктор
        for (int i = 0; i < numParticles; ++i) {
            predictedParticles[i] = particles[i];

            // Предсказываем скорости (на полный шаг для простоты)
            double ax_old = sum_velocity_x_old[i] + f_halo_x_old[i];
            double ay_old = sum_velocity_y_old[i] + f_halo_y_old[i];
            predictedParticles[i].velocityX += ax_old * tau;
            predictedParticles[i].velocityY += ay_old * tau;

            // Предсказываем позиции (используем исходные скорости + полушаг ускорения)
            predictedParticles[i].x += particles[i].velocityX * tau + 0.5 * ax_old * tau * tau;
            predictedParticles[i].y += particles[i].velocityY * tau + 0.5 * ay_old * tau * tau;
            predictedParticles[i].distanceFromCenter = std::hypot(predictedParticles[i].x, predictedParticles[i].y);

            // Предсказываем гидродинамические величины
            predictedParticles[i].density += sum_rho_old[i] * tau;
            predictedParticles[i].energy += sum_energy_old[i] * tau;
            predictedParticles[i].pressure = (gamma - 1.0) * predictedParticles[i].energy * predictedParticles[i].density;
        }

        // Шаг 3: Вычисляем силы для предсказанного состояния
        ComputeForces(predictedParticles, sum_rho, sum_velocity_x, sum_velocity_y, sum_energy, radius);

        // Шаг 4: Корректор
        for (int i = 0; i < numParticles; ++i) {
            // Вычисляем гравитацию гало для предсказанных частиц
            double ksi_pred = computeKsiForHalo(predictedParticles[i].x, predictedParticles[i].y, 0.0);
            double f_halo_x_pred = computeForceGrav(ksi_pred, predictedParticles[i].x, A, a_halo);
            double f_halo_y_pred = computeForceGrav(ksi_pred, predictedParticles[i].y, A, b_halo);

            // Ускорения: старые и предсказанные
            double ax_old = sum_velocity_x_old[i] + f_halo_x_old[i];
            double ay_old = sum_velocity_y_old[i] + f_halo_y_old[i];
            double ax_pred = sum_velocity_x[i] + f_halo_x_pred;
            double ay_pred = sum_velocity_y[i] + f_halo_y_pred;

            // Корректируем скорости (усреднение ускорений)
            particles[i].velocityX = particles[i].velocityX + 0.5 * (ax_old + ax_pred) * tau;
            particles[i].velocityY = particles[i].velocityY + 0.5 * (ay_old + ay_pred) * tau;

            // Корректируем позиции (усреднение скоростей)
            particles[i].x = particles[i].x + 0.5 * tau * (particles[i].velocityX + predictedParticles[i].velocityX);
            particles[i].y = particles[i].y + 0.5 * tau * (particles[i].velocityY + predictedParticles[i].velocityY);
            particles[i].distanceFromCenter = std::hypot(particles[i].x, particles[i].y);

            // Корректируем гидродинамические величины
            particles[i].density = particles[i].density + 0.5 * (sum_rho_old[i] + sum_rho[i]) * tau;
            particles[i].energy = particles[i].energy + 0.5 * (sum_energy_old[i] + sum_energy[i]) * tau;
            particles[i].pressure = (gamma - 1.0) * particles[i].energy * particles[i].density;
            // Если есть sum_smoothingLength, то:
            // particles[i].smoothingLength = particles[i].smoothingLength + 0.5 * (sum_smoothingLength_old[i] + sum_smoothingLength[i]) * tau;
        }

        // Шаг 5: Обновление времени и вывод
        t += tau;
        std::cout << "Time: " << t << " | dt: " << tau << std::endl;

        if (t >= Dt) {
            if (check_print == 0) {
                check_print = 1;
            }
            else {
                saveParticlesToFile(t, Dt);
            }
            std::cout << t << " DT: " << Dt << std::endl;
            Dt += shag_dt;
        }
    }
}