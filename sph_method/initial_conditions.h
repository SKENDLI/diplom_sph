﻿#pragma once
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
    A_G = AMS / (RCORE1 - atan(RCORE1));
    B_rho = system_Mass / fun_mass(0.0, scale_radius, 1000, pi_2, scale_halo);
    double rd_mah = 2.0;
    double B_p = B_rho / gamma * rd_mah * rd_mah * A_G / (Mah_d * Mah_d);
    A = 4.0 * M_PI * G * B_rho * scale_halo * scale_halo;
    double gamma1 = gamma - 1.0;
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
    int Np = N_fi[k];
    int Nk = k;
    std::cout << "Np = " << Np << "\n";
    std::cout << "k_max = " << Nk << "\n";
    std::cout << "rk_max = " << rk[Nk] << "\n";
    particles.resize(Np);
    std::vector<double> h_p(Np);

    hp_max = Lh * (rk[Nk] - rk[Nk - 1]);
    double hp_min = hp_max;
    double Bh = (2.0 * Nk - 0.5 * (Nk - 5)) / 5.0;
    double Ah = (0.5 - Bh) / Nk;

    int particleIndex = 0;
    for (int ringIdx = 2; ringIdx <= Nk; ++ringIdx) {
        double r_min = rk[ringIdx - 1];
        double r_max = rk[ringIdx];
        double r_avg = (r_min + r_max) / 2.0;
        int numParticlesInRing = Nk_fi[ringIdx];
        double phi_step = pi_2 / numParticlesInRing;

        for (int j = 0; j < numParticlesInRing && particleIndex < Np; ++j) {
            double phi = (j - 1 - N_fi[ringIdx]) * phi_step;
            double x = r_avg * cos(phi);
            double y = r_avg * sin(phi);
            double distance = sqrt(x * x + y * y);

            h_p[particleIndex] = Lh * (r_max - r_min);
            if (h_p[particleIndex] > hp_max) h_p[particleIndex] = hp_max;
            if (h_p[particleIndex] < hp_min) hp_min = h_p[particleIndex];

            double density = B_rho * computeRho(scale_halo, distance);
            double pressure_local = B_p * pow(density / B_rho, gamma);
            double energy = Ap * pow(density, gamma1) / gamma1;

            particles[particleIndex].density = B_rho * computeRho(scale_halo, distance);
            particles[particleIndex].pressure = B_p * pow(particles[particleIndex].density / B_rho, gamma);
            particles[particleIndex].energy = Ap * pow(particles[particleIndex].density, gamma1) / gamma1;
            particles[particleIndex].mass = particleMass;
            particles[particleIndex].smoothingLength = h_p[particleIndex];
            particles[particleIndex].velocityX = 0.0;
            particles[particleIndex].velocityY = 0.0;
            particles[particleIndex].neighbors = 0;
            particles[particleIndex].x = x;
            particles[particleIndex].y = y;

            particleIndex++;
        }
    }

    h = hp_min;
    h2 = 2.0 * h;

    Ngx = static_cast<int>((x_ex - x_in) / h2 + 0.5);
    Ngy = static_cast<int>((y_ex - y_in) / h2 + 0.5);

    std::cout << "Ngx = " << Ngx << ", Ngy = " << Ngy << "\n";

    std::cout << Np << std::endl;
    particles.resize(Np);
    predictedParticles.resize(Np);
    neighborCount.resize(Np);      // Оставляем, если нужно считать соседей
    Rho.resize(Np);                // Начальная плотность
    Rhot.resize(Np);               // Вычисляемая плотность
    Fu_p.resize(Np);               // Сила по X
    Fv_p.resize(Np);               // Сила по Y
    Fe_p.resize(Np);               // Изменение энергии
    ht_p.resize(Np);               // Обновляемая сглаживающая длина
    force_halo_x.resize(Np);       // Сила от гало по X (опционально)
    force_halo_y.resize(Np);       // Сила от гало по Y (опционально)
    initialParticles.resize(Np);
    numParticles = Np;

    int i = 0;
    for (const auto& p : particles) {
        max_x = max(p.x, max_x);
        ht_p[i] = particles[i].smoothingLength;
        Rho[i] = particles[i].smoothingLength;
        i++;
    }
    std::cout << "Max X: " << max_x << std::endl;
    
    // Вычисляем гидродинамические силы
    ComputeForces(particles);

    // Вычисляем гравитационные силы гало
    computeForcesGravCool(particles, force_halo_x, force_halo_y);

    i = 0;
    for (auto& p : particles) {
        double x = p.x;
        double y = p.y;
        double distance = std::hypot(x, y);  // Используем hypot для большей точности
        double ff = std::atan2(y, x);

        // Радиальная составляющая гравитационной силы
        double F_r = (force_halo_x[i] * x + force_halo_y[i] * y) / distance;

        // Радиальная составляющая гидродинамической силы
        double Fuv_r = (Fu_p[i] * x + Fv_p[i] * y) / distance;

        // Скорость вращения
        double v_rot = 0.0;
        double Fr_total = F_r + Fuv_r;
        v_rot = std::sqrt(-distance * Fr_total);

        // Задаём скорости
        p.velocityX = -v_rot * std::sin(ff);
        p.velocityY = v_rot * std::cos(ff);

        // Давление (если нужно для вывода или дальнейших расчётов)
        double Ppi = (gamma - 1.0) * Rho[i] * p.energy;

        i++;
    }
    ComputeForces(particles);
}
