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

void SPH()
{
    double temp = 0.0;
    int check = 0;
    int check_print = 0;
    int count = 0;
    int count1 = 0;
    double A = 4.0 * M_PI * G * B_rho * scale_halo * scale_halo;
    while (t <= t_end)
    {
        Predictor();
        dt();
        for (int i = 0; i < numParticles; ++i)
        {
            if (particles[i].mass > 0)
            {
                double ksi = computeKsiForHalo(particles[i].x, particles[i].y, 0.0);
                double force_halo_x = computeForceGrav(ksi, particles[i].x, A, 1.0);
                double force_halo_y = computeForceGrav(ksi, particles[i].y, A, 1.0);

                // Обновляем скорость и положение (предиктор)
                predictedParticles[i].velocityX = particles[i].velocityX + (sum_velocity_x[i] + force_halo_x) * 0.5 * tau;
                predictedParticles[i].velocityY = particles[i].velocityY + (sum_velocity_y[i] + force_halo_y) * 0.5 * tau;
                predictedParticles[i].x = particles[i].x + predictedParticles[i].velocityX * 0.5 * tau;
                predictedParticles[i].y = particles[i].y + predictedParticles[i].velocityY * 0.5 * tau;

                // Исправлено: вычисляем расстояние от центра по новым координатам
                predictedParticles[i].distanceFromCenter = sqrt(predictedParticles[i].x * predictedParticles[i].x + predictedParticles[i].y * predictedParticles[i].y);

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
                double ksi = computeKsiForHalo(predictedParticles[i].x, predictedParticles[i].y, 0.0);
                double force_halo_x = computeForceGrav(ksi, predictedParticles[i].x, A, 1.0);
                double force_halo_y = computeForceGrav(ksi, predictedParticles[i].y, A, 1.0);

                // Обновляем скорость и положение (корректор)
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
        checkEnergyConservation(initialParticles, particles, G);
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
