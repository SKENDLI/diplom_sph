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
            checkEnergyConservation(initialParticles, particles, G);
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
    // Добавляем временные массивы для хранения старых ускорений
    std::vector<double> sum_velocity_x_old(numParticles);
    std::vector<double> sum_velocity_y_old(numParticles);
    std::vector<double> f_halo_x_old(numParticles);
    std::vector<double> f_halo_y_old(numParticles);

    int check_print = 0;
    while (t <= t_end) {
        ComputeForces(particles, sum_rho, sum_velocity_x, sum_velocity_y, sum_energy, radius);

        // Сохраняем ускорения исходных частиц
        sum_velocity_x_old = sum_velocity_x;
        sum_velocity_y_old = sum_velocity_y;

        // Вычисляем гравитацию гало для исходных частиц
#pragma omp parallel for
        for (int i = 0; i < numParticles; ++i) {
            //if (particles[i].distanceFromCenter >= radius) continue;
            double ksi = computeKsiForHalo(particles[i].x, particles[i].y, 0.0);
            f_halo_x_old[i] = computeForceGrav(ksi, particles[i].x, A, a_halo);
            f_halo_y_old[i] = computeForceGrav(ksi, particles[i].y, A, b_halo);
        }
        //TODO УБРАТЬ ГРАНИЦЫ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        dt();  // Обновляем tau

        // Этап 3: Предиктор (полушаг)
#pragma omp parallel for
        for (int i = 0; i < numParticles; ++i) {
            //if (particles[i].distanceFromCenter >= radius) continue;

            // Используем сохранённые ускорения исходных частиц
            predictedParticles[i] = particles[i];
            predictedParticles[i].velocityX += (sum_velocity_x_old[i] + f_halo_x_old[i]) * 0.5 * tau;
            predictedParticles[i].velocityY += (sum_velocity_y_old[i] + f_halo_y_old[i]) * 0.5 * tau;
            predictedParticles[i].x += predictedParticles[i].velocityX * 0.5 * tau;
            predictedParticles[i].y += predictedParticles[i].velocityY * 0.5 * tau;
            predictedParticles[i].distanceFromCenter = std::hypot(predictedParticles[i].x, predictedParticles[i].y);
        }

        ComputeForces(predictedParticles, sum_rho, sum_velocity_x, sum_velocity_y, sum_energy, radius);

        // Этап 5: Корректор (полный шаг)
#pragma omp parallel for
        for (int i = 0; i < numParticles; ++i) {
            //if (particles[i].distanceFromCenter >= radius) continue;

            // Вычисляем гравитацию гало для предсказанных частиц
            double ksi_pred = computeKsiForHalo(predictedParticles[i].x, predictedParticles[i].y, 0.0);
            double f_halo_x_pred = computeForceGrav(ksi_pred, predictedParticles[i].x, A, a_halo);
            double f_halo_y_pred = computeForceGrav(ksi_pred, predictedParticles[i].y, A, b_halo);

            // Вычисляем ускорения для корректора
            double a_old_x = sum_velocity_x_old[i] + f_halo_x_old[i];
            double a_pred_x = sum_velocity_x[i] + f_halo_x_pred;
            double a_old_y = sum_velocity_y_old[i] + f_halo_y_old[i];
            double a_pred_y = sum_velocity_y[i] + f_halo_y_pred;

            // Обновляем скорость с учётом среднего ускорения
            particles[i].velocityX = predictedParticles[i].velocityX + (a_pred_x + a_old_x) * tau / 2;
            particles[i].velocityY = predictedParticles[i].velocityY + (a_pred_y + a_old_y) * tau / 2;

            // Обновляем положение с усреднением скоростей
            particles[i].x = predictedParticles[i].x + 0.5 * tau * (particles[i].velocityX + predictedParticles[i].velocityX);
            particles[i].y = predictedParticles[i].y + 0.5 * tau * (particles[i].velocityY + predictedParticles[i].velocityY);
            particles[i].distanceFromCenter = std::hypot(particles[i].x, particles[i].y);
        }

        // Этап 6: Обновление параметров
#pragma omp parallel for
        for (int i = 0; i < numParticles; ++i) {
            particles[i].density += sum_rho[i] * tau;
            particles[i].energy += sum_energy[i] * tau;
            particles[i].pressure = (gamma - 1.0) * particles[i].energy * particles[i].density;
            particles[i].smoothingLength += sum_smoothingLength[i] * tau;
        }
        
        t += tau;
        std::cout << "Time: " << t << " | dt: " << tau << std::endl;
        //checkEnergyConservation(initialParticles, particles, G);
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