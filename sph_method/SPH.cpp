#include "calculate_functions.h"
#include "initial_conditions.h"
#include "in_out_data_functions.h"

int main()
{
    initial_energy = 0.0;
    int checkConfiguration = configuration();
    if (checkConfiguration != 0)
    {
        particles.resize(numParticles);
        predictedParticles.resize(numParticles);
        neighborCount.resize(numParticles);
        Rho.resize(numParticles);            // Начальная плотность
        Rhot.resize(numParticles);           // Вычисляемая плотность
        Fu_p.resize(numParticles);           // Сила по X
        Fv_p.resize(numParticles);           // Сила по Y
        Fe_p.resize(numParticles);           // Изменение энергии
        ht_p.resize(numParticles);           // Обновляемая сглаживающая длина
        force_halo_x.resize(numParticles);   // Сила от гало по X
        force_halo_y.resize(numParticles);   // Сила от гало по Y
    }

    string path = "";
    int choice;
    while (true)
    {
        cout << "\n1 - start the program with initialization of new values\n2 - launch the program with initialization of values from the file\n";
        cin >> choice;
        switch (choice)
        {
        case 1:
            system("CLS");
            initializeGasCloud();
            saveParticlesToFile(t, Dt);
            cout << Dt << endl;
            law_cons(initial_energy, particles);
            cout << "Initial Energy: " << initial_energy << endl;
            SPH();
            return 0;
        case 2:
            system("CLS");
            cout << "\nEnter the folder and the file name\n";
            cin >> path;
            path = "data\\" + path;
            //FillParticles(path);
            //SPH();
            return 0;
        default:
            system("CLS");
            cout << "Something went wrong.";
            return 0;
        }
    }
    return 0;
}

void SPH() {
    vector<double> Fu_p_old(numParticles);
    vector<double> Fv_p_old(numParticles);
    vector<double> Fe_p_old(numParticles);
    double dt_x = 1000, dt_y = 1000;
    int check_print = 0;
    double next_print_time = 0.0;
    tau = 1e-3;

    while (t <= t_end) {
        // Шаг 1: Вычисляем силы для текущего состояния
        ComputeForces(particles);  // Вычисляет Fu_p, Fv_p, Fe_p, Rhot, ht_p

        // Сохраняем текущие силы
        Fu_p_old = Fu_p;
        Fv_p_old = Fv_p;
        Fe_p_old = Fe_p;

        // Вычисляем гравитационные силы гало
        computeForcesGravCool(particles, force_halo_x, force_halo_y);

        // Шаг 2: Предиктор (t + dt/2)
        double dt05 = 0.5 * tau;
        for (int i = 0; i < numParticles; ++i) {
            predictedParticles[i] = particles[i];

            // Полярные координаты (предиктор)
            double xi = predictedParticles[i].x;
            double yi = predictedParticles[i].y;
            double ri = sqrt(xi * xi + yi * yi);
            double cxy = xi / ri;
            double sxy = yi / ri;
            double Vx = predictedParticles[i].velocityX;
            double Vy = predictedParticles[i].velocityY;
            double Vfi = Vy * cxy - Vx * sxy;

            // Силы
            double Fx = force_halo_x[i] + Fu_p_old[i];
            double Fy = force_halo_y[i] + Fv_p_old[i];
            double Fr = Fx * cxy + Fy * sxy + Vfi * Vfi / ri;
            double Ffi = (Fy * cxy - Fx * sxy) * ri;

            // Корректор (текущее состояние)
            xi = particles[i].x;
            yi = particles[i].y;
            ri = sqrt(xi* xi + yi * yi);
            cxy = xi / ri;
            sxy = yi / ri;
            double fi = atan2(yi, xi);
            Vx = particles[i].velocityX;
            Vy = particles[i].velocityY;
            double Vr = Vx * cxy + Vy * sxy;
            Vfi = Vy * cxy - Vx * sxy;
            double Lrv = ri * Vfi;

            // Предиктор на t + dt/2
            double Vr_t = Vr + dt05 * Fr;
            double ri_t = ri + dt05 * (Vr_t + Vr);
            double Lrv_t = Lrv + dt05 * Ffi;
            double Vfi_t = Lrv_t / ri_t;
            double fi_t = fi + dt05 * (Vfi / ri + Vfi_t / ri_t);
            double cfi = cos(fi_t);
            double sfi = sin(fi_t);

            // Обновляем предсказанные значения
            predictedParticles[i].x = ri_t * cfi;
            predictedParticles[i].y = ri_t * sfi;
            predictedParticles[i].velocityX = Vr_t * cfi - Vfi_t * sfi;
            predictedParticles[i].velocityY = Vfi_t * cfi + Vr_t * sfi;
            predictedParticles[i].energy = particles[i].energy + dt05 * Fe_p_old[i];
            predictedParticles[i].distanceFromCenter = sqrt(predictedParticles[i].x * predictedParticles[i].x + predictedParticles[i].y * predictedParticles[i].y);
            predictedParticles[i].pressure = (gamma - 1.0) * predictedParticles[i].energy * Rhot[i];
        }

        // Шаг 3: Вычисляем силы для предсказанного состояния (t + dt/2)
        ComputeForces(predictedParticles);  // Обновляет Fu_p, Fv_p, Fe_p, Rhot, ht_p
        computeForcesGravCool(predictedParticles, force_halo_x, force_halo_y);

        // Шаг 4: Корректор (t + dt) и вычисление временного шага

        for (int i = 0; i < numParticles; ++i) {
            // Полярные координаты (предиктор)
            double xi = predictedParticles[i].x;
            double yi = predictedParticles[i].y;
            double ri = sqrt(xi * xi + yi * yi);
            double cxy = xi / ri;
            double sxy = yi / ri;
            double Vx = predictedParticles[i].velocityX;
            double Vy = predictedParticles[i].velocityY;
            double Vfi = Vy * cxy - Vx * sxy;

            // Силы
            double Fx = force_halo_x[i] + Fu_p[i];
            double Fy = force_halo_y[i] + Fv_p[i];
            double Fr = Fx * cxy + Fy * sxy + Vfi * Vfi / ri;
            double Ffi = (Fy * cxy - Fx * sxy) * ri;

            // Корректор (текущее состояние)
            xi = particles[i].x;
            yi = particles[i].y;
            ri = sqrt(xi * xi + yi * yi);
            cxy = xi / ri;
            sxy = yi / ri;
            double fi = atan2(yi, xi);
            Vx = particles[i].velocityX;
            Vy = particles[i].velocityY;
            double Vr = Vx * cxy + Vy * sxy;
            Vfi = Vy * cxy - Vx * sxy;
            double Lrv = ri * Vfi;

            // Корректор на t + dt
            double Vr_tt = Vr + dt05 * Fr;
            double ri_tt = ri + dt05 * Vr_tt;
            double Vr_t = 2.0 * Vr_tt - Vr;
            double ri_t = 2.0 * ri_tt - ri;
            double Lrv_tt = Lrv + dt05 * Ffi;
            double Lrv_t = 2.0 * Lrv_tt - Lrv;
            double Vfi_t = Lrv_t / ri_t;
            double fi_t = fi + tau * Lrv_tt / (ri_tt * ri_tt);
            double cfi = cos(fi_t);
            double sfi = sin(fi_t);

            // Обновляем значения
            particles[i].x = ri_t * cfi;
            particles[i].y = ri_t * sfi;
            particles[i].velocityX = Vr_t * cfi - Vfi_t * sfi;
            particles[i].velocityY = Vfi_t * cfi + Vr_t * sfi;
            particles[i].energy = particles[i].energy + tau * Fe_p[i];
            particles[i].distanceFromCenter = sqrt(particles[i].x * particles[i].x + particles[i].y * particles[i].y);
            particles[i].pressure = (gamma - 1.0) * particles[i].energy * Rhot[i];

            // Обновление временного шага (как в Fortran)
            double hpi = ht_p[i];
            double ax = sqrt(hpi / abs(Fx));
            double ay = sqrt(hpi / abs(Fy));
            if (dt_x > ax) dt_x = ax;
            if (dt_y > ay) dt_y = ay;
        }
        double min = min(dtn, dt_x);
        tau = courant * min(min, dt_y);

        // Шаг 5: Обновление времени и вывод
        t += tau;
        cout << "Time: " << t << " | dt: " << tau << endl;

        save_laws_cons(t, next_print_time, initial_energy, particles);

        if (t >= Dt) {
            if (check_print == 0) {
                check_print = 1;
            }
            else {
                saveParticlesToFile(t, Dt);
            }
            cout << t << " DT: " << Dt << endl;
            Dt += shag_dt;
        }
    }
}