#pragma once
void initializeGasCloud() {
    particles.resize(numParticles);
    double v_rot;
    double particleMass = system_Mass / numParticles;
    double scale_radius = radius / scale;
    int numRings = 100;
    double ringWidth = radius / numRings;

    double RCORE1 = rmax_hl / r_core;
    double A_G = AMS / (RCORE1 - atan(RCORE1));
    B_rho = system_Mass / fun_mass(r_core, rmax_hl, 1000, M_PI, scale_radius);
    double rd_mah = radius;
    double B_p = B_rho / gamma * rd_mah * rd_mah * A_G / (Mah_d * Mah_d);

    double A = 4.0 * M_PI * G * B_rho * scale_halo * scale_halo;


    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    int particleIndex = 0;

    for (int ringIdx = 0; ringIdx < numRings; ++ringIdx) {
        double r_min = ringIdx * ringWidth;
        double r_max = r_min + ringWidth;
        double r_avg = (r_min + r_max) / 2.0;

        double numParticlesInRing = (numParticles * exp(-r_avg / scale_radius)) /
            (1.0 - (1.0 + 1.0 / scale_radius) * exp(-1.0 / scale_radius)) *
            ((1 + (ringIdx - 1) * (ringWidth / scale_radius)) *
                (exp(ringWidth / scale_radius) - 1) - (ringWidth / scale_radius));

        numParticlesInRing = std::round(numParticlesInRing);

        for (int j = 0; j < numParticlesInRing; ++j) {
            if (particleIndex >= numParticles) break;

            double r = r_min + sqrt(dis(gen)) * (r_max - r_min);
            double phi = 2.0 * M_PI * dis(gen);
            double x = r * cos(phi);
            double y = r * sin(phi);
            double distance = sqrt(x * x + y * y);
            // Вычисляем плотность и давление
            double density = B_rho * computeRho(distance, scale_radius);
            double pressure_local = B_p * pow(density / B_rho, gamma);

            double dRho_dr = -2.0 * sinh(distance / scale_radius) / (scale_radius * pow(cosh(distance / scale_radius), 3));
            // Градиент давления
            double dP_dr = gamma * pressure * pow(density, gamma - 2.0) * dRho_dr;

            // Сила давления
            double force_pressure_x = -dP_dr * (x / distance);
            double force_pressure_y = -dP_dr * (y / distance);

            // Гравитационные силы
            double ksi = computeKsiForHalo(x, y, 0.0);
            double force_halo_x = computeForceGrav(ksi, x, A, a_halo);
            double force_halo_y = computeForceGrav(ksi, y, A, b_halo);

            double force_gravity_x = force_halo_x;
            double force_gravity_y = force_halo_y;

            // Общая сила (гравитация + давление)
            double F_r = (force_gravity_x * x + force_gravity_y * y) / distance;
            double Fuv_r = (force_pressure_x * x + force_pressure_y * y) / distance;

            // Расчет скорости вращения
            v_rot = sqrt(distance * (F_r + Fuv_r));

            // Угол для вычисления компонент скорости
            double ff = atan2(y, x);
            particles[particleIndex].velocityX = -v_rot * sin(ff);
            particles[particleIndex].velocityY = v_rot * cos(ff);

            // Заполнение данных частицы
            particles[particleIndex].x = x;
            particles[particleIndex].y = y;
            particles[particleIndex].distanceFromCenter = distance;
            particles[particleIndex].density = density;
            particles[particleIndex].mass = particleMass;
            particles[particleIndex].smoothingLength = sqrt(particleMass / density);
            particles[particleIndex].pressure = pressure_local;
            particles[particleIndex].energy = (pressure_local * pow(density, gamma - 1.0)) / (gamma - 1.0);

            ++particleIndex;
        }
    }
    // Заполнение оставшихся частиц в центре
    double r_center = r_core * radius; // Радиус центральной области
    while (particleIndex < numParticles) {
        // Распределение частиц в пределах малого радиуса
        double r = sqrt(dis(gen)) * r_center;

        double phi = 2.0 * M_PI * dis(gen);
        double x = r * cos(phi);
        double y = r * sin(phi);
        double distance = sqrt(x * x + y * y);
        // Вычисляем плотность и давление
        double density = B_rho * computeRho(distance, scale_radius);
        double pressure_local = B_p * pow(density / B_rho, gamma);

        double dRho_dr = -2.0 * sinh(distance / scale_radius) / (scale_radius * pow(cosh(distance / scale_radius), 3));
        // Градиент давления
        double dP_dr = gamma * pressure * pow(density, gamma - 2.0) * dRho_dr;

        // Сила давления
        double force_pressure_x = -dP_dr * (x / distance);
        double force_pressure_y = -dP_dr * (y / distance);

        // Гравитационные силы
        double ksi = computeKsiForHalo(x, y, 0.0);
        double force_halo_x = computeForceGrav(ksi, x, A, a_halo);
        double force_halo_y = computeForceGrav(ksi, y, A, b_halo);

        double force_gravity_x = force_halo_x;
        double force_gravity_y = force_halo_y;

        // Общая сила (гравитация + давление)
        double F_r = (force_gravity_x * x + force_gravity_y * y) / distance;
        double Fuv_r = (force_pressure_x * x + force_pressure_y * y) / distance;

        // Расчет скорости вращения
        v_rot = sqrt(distance * (F_r + Fuv_r));

        // Угол для вычисления компонент скорости
        double ff = atan2(y, x);
        particles[particleIndex].velocityX = -v_rot * sin(ff);
        particles[particleIndex].velocityY = v_rot * cos(ff);

        // Заполнение данных частицы
        particles[particleIndex].x = x;
        particles[particleIndex].y = y;
        particles[particleIndex].distanceFromCenter = distance;
        particles[particleIndex].density = density;
        particles[particleIndex].mass = particleMass;
        particles[particleIndex].smoothingLength = sqrt(particleMass / density);
        particles[particleIndex].pressure = pressure_local;
        particles[particleIndex].energy = (pressure_local * pow(density, gamma - 1.0)) / (gamma - 1.0);

        ++particleIndex;
    }
}
