#pragma once
#pragma once
void initializeGasCloud() {
    std::vector<double> rk(1000, 0.0);  // Радиусы колец
    std::vector<int> N_fi(1000, 0);     // Общее число частиц до данного кольца
    std::vector<int> Nk_fi(1000, 0);    // Число частиц в текущем кольце
    double particleMass = system_Mass / numParticles;
    double scale_radius = radius * scale;
    double drr = 1e-5;
    int Nm = 20;
    double pi_2 = 2.0 * M_PI;

    double RCORE1 = rmax_hl / r_core;
    double A_G = AMS / (RCORE1 - atan(RCORE1));
    B_rho = system_Mass / fun_mass(0.0, scale_radius, 1000, pi_2, scale_halo);
    double rd_mah = 2.0;
    double B_p = B_rho / gamma * rd_mah * rd_mah * A_G / (Mah_d * Mah_d);
    A = 4.0 * M_PI * G * B_rho * scale_halo * scale_halo;

    // Инициализация первого кольца
    N_fi[1] = 0;
    Nk_fi[1] = 0;

    double rr = 0.01;
    double ff_rr = B_rho * fun_mass(0.0, rr, Nm, pi_2, scale_halo) / particleMass - pi_2;
    double rr_1 = rr + drr;
    double ff_rr1 = B_rho * fun_mass(0.0, rr_1, Nm, pi_2, scale_halo) / particleMass - pi_2;
    double dff = drr * ff_rr / (ff_rr1 - ff_rr);

    while (std::abs(ff_rr) > 1e-6 || std::abs(dff) > 1e-6) {
        rr -= dff;
        rr_1 = rr + drr;
        ff_rr = B_rho * fun_mass(0.0, rr, Nm, pi_2, scale_halo) / particleMass - pi_2;
        ff_rr1 = B_rho * fun_mass(0.0, rr_1, Nm, pi_2, scale_halo) / particleMass - pi_2;
        dff = drr * ff_rr / (ff_rr1 - ff_rr);
    }

    rk[2] = rr;
    Nk_fi[2] = 3;
    N_fi[2] = 3;

    int k = 2;

    while (rk[k] <= radius && N_fi[k] < numParticles) {
        double rrk = rk[k];
        double rrk_1 = rk[k - 1];
        rr = rrk;
        rr_1 = rr + drr;
        ff_rr = B_rho * fun_mass(rrk, rr, Nm, pi_2, scale_halo) / particleMass - pi_2 * (rr + rrk) / (rr - rrk_1);
        ff_rr1 = B_rho * fun_mass(rrk, rr_1, Nm, pi_2, scale_halo) / particleMass - pi_2 * (rr_1 + rrk) / (rr_1 - rrk_1);
        dff = drr * ff_rr / (ff_rr1 - ff_rr);

        while (std::abs(ff_rr) > 1e-5 || std::abs(dff) > 1e-6) {
            rr -= dff;
            rr_1 = rr + drr;
            ff_rr = B_rho * fun_mass(rrk, rr, Nm, pi_2, scale_halo) / particleMass - pi_2 * (rr + rrk) / (rr - rrk_1);
            ff_rr1 = B_rho * fun_mass(rrk, rr_1, Nm, pi_2, scale_halo) / particleMass - pi_2 * (rr_1 + rrk) / (rr_1 - rrk_1);
            dff = drr * ff_rr / (ff_rr1 - ff_rr);
        }

        rk[k + 1] = rr;
        Nk_fi[k + 1] = static_cast<int>(2.0 * M_PI * (rr + rrk) / (rr - rrk_1) + 0.5);
        N_fi[k + 1] = N_fi[k] + Nk_fi[k + 1];
        k++;
    }

    // Размещение частиц в кольцах
    particles.clear();
    int particleIndex = 0;

    for (int ringIdx = 2; ringIdx < k; ++ringIdx) {
        double r_min = rk[ringIdx - 1];
        double r_max = rk[ringIdx];
        double r_avg = (r_min + r_max) / 2.0;
        int numParticlesInRing = Nk_fi[ringIdx];
        double phi_step = pi_2 / numParticlesInRing;

        for (int j = 0; j < numParticlesInRing; ++j) {
            if (particleIndex >= numParticles) break;

            double phi = j * phi_step;
            double x = r_avg * cos(phi);
            double y = r_avg * sin(phi);
            double distance = sqrt(x * x + y * y);

            // Вычисление плотности и давления
            double density = B_rho * computeRho(scale_halo, distance);
            double pressure_local = B_p * pow(density / B_rho, gamma);

            // Создание частицы
            Particle p;
            p.x = x;
            p.y = y;
            p.distanceFromCenter = distance;
            p.density = density;
            p.mass = particleMass;
            p.smoothingLength = 0.04;
            p.pressure = pressure_local;
            p.energy = Ap * pow(density, gamma - 1.0) / (gamma - 1.0);
            p.velocityX = 0.0;
            p.velocityY = 0.0;
            particles.push_back(p);
            ++particleIndex;
        }
    }
    std::cout << particleIndex << std::endl;
    particles.resize(particleIndex);
    predictedParticles.resize(particleIndex);
    neighborCount.resize(particleIndex);
    sum_rho.resize(particleIndex);
    sum_energy.resize(particleIndex);
    sum_velocity_x.resize(particleIndex);
    sum_velocity_y.resize(particleIndex);
    sum_smoothingLength.resize(particleIndex);
    force_halo_x.resize(particleIndex);
    force_halo_y.resize(particleIndex);
    initialParticles.resize(particleIndex);
    numParticles = particleIndex;


    // Поиск максимального x
    for (const auto& p : particles) {
        if (p.x > max_x) {
            max_x = p.x;
        }
    }
    std::cout << "Max X: " << max_x << std::endl;

    // Вычисление сил
    ComputeForces(particles, sum_rho, sum_velocity_x, sum_velocity_y, sum_energy, radius);
    int i = 0;
    // Установка вращательных скоростей
    for (auto& p : particles) {
        double x = p.x;
        double y = p.y;
        double distance = sqrt(x * x + y * y);
        double ff = atan2(y, x);

        double ksi = computeKsiForHalo(x, y, 0.0);
        double force_halo_x = computeForceGrav(ksi, x, A, a_halo);
        double force_halo_y = computeForceGrav(ksi, y, A, b_halo);
        double F_r = (force_halo_x * x + force_halo_y * y) / distance;
        double Fuv_r = (sum_velocity_x[i] * x + sum_velocity_y[i] * y) / distance;

        double v_rot = sqrt(distance * (F_r + Fuv_r));
        p.velocityX = -v_rot * sin(ff);
        p.velocityY = v_rot * cos(ff);
        i++;
    }
}
