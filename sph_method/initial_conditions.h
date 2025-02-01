#pragma once
void initializeGasCloud() {
    particles.resize(numParticles);

    double particleMass = system_Mass / numParticles;
    double scale_radius = radius / scale;
    int numRings = 100;
    double ringWidth = radius / numRings;

    double A = 4.0 * M_PI * G * density0 * scale_halo * scale_halo;
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

            double ksi = computeKsiForHalo(x, y, 0.0);
            double force_halo_x = computeForceGrav(ksi, x, A, a_halo);
            double force_halo_y = computeForceGrav(ksi, y, A, b_halo);

            double force_gravity_x = force_halo_x;
            double force_gravity_y = force_halo_y;

            double density_disk = computeRho(distance, scale_radius);
            double density = density_disk * density0;
            double pressure_local = pressure * pow(density, gamma);

            double dP_dr = -gamma * pressure * pow(density, gamma - 1.0) * (density / distance);
            double force_pressure_x = dP_dr * (x / distance);
            double force_pressure_y = dP_dr * (y / distance);

            double force_total_x = force_gravity_x + force_pressure_x;
            double force_total_y = force_gravity_y + force_pressure_y;

            particles[particleIndex].x = x;
            particles[particleIndex].y = y;
            particles[particleIndex].distanceFromCenter = distance;
            particles[particleIndex].density = density;
            particles[particleIndex].mass = particleMass;
            particles[particleIndex].smoothingLength = sqrt(particleMass / density);
            particles[particleIndex].pressure = pressure_local;
            particles[particleIndex].energy = (pressure * pow(density, gamma - 1.0)) / (gamma - 1.0);

            // Расчет скорости вращения
            double v_rot = sqrt(distance * sqrt(force_total_x * force_total_x + force_total_y * force_total_y));
            double ff = atan2(y, x);

            // Составление компонент скорости
            particles[particleIndex].velocityX = -v_rot * sin(ff);
            particles[particleIndex].velocityY = v_rot * cos(ff);

            ++particleIndex;
        }
    }

    // Заполнение оставшихся частиц равномерно в центре
    double r_center = 0.1 * radius; // Радиус области для равномерного распределения
    while (particleIndex < numParticles) {
        double r = sqrt(dis(gen)) * r_center; // Радиус равномерно по площади
        double phi = 2.0 * M_PI * dis(gen);
        double x = r * cos(phi);
        double y = r * sin(phi);
        double distance = sqrt(x * x + y * y);

        // Вычисление плотности и давления для частиц в центре
        double density_disk = computeRho(distance, scale_radius);
        double density = density_disk * density0;
        double pressure_local = pressure * pow(density, gamma);

        // Расчет градиента давления
        double dP_dr = -gamma * pressure * pow(density, gamma - 1.0) * (density / distance);
        double force_pressure_x = dP_dr * (x / distance);
        double force_pressure_y = dP_dr * (y / distance);

        // Вычисление гравитационных сил
        double ksi = computeKsiForHalo(x, y, 0.0); // z = 0, если 2D модель
        double force_halo_x = computeForceGrav(ksi, x, A, a_halo); // Учитываем ξ при вычислении силы
        double force_halo_y = computeForceGrav(ksi, y, A, b_halo); // Учитываем ξ при вычислении силы

        // Итоговая гравитационная сила
        double force_gravity_x = force_halo_x;
        double force_gravity_y = force_halo_y;

        // Общая сила
        double force_total_x = force_gravity_x + force_pressure_x;
        double force_total_y = force_gravity_y + force_pressure_y;

        // Заполнение параметров частицы
        particles[particleIndex].x = x;
        particles[particleIndex].y = y;
        particles[particleIndex].distanceFromCenter = distance;
        particles[particleIndex].density = density;
        particles[particleIndex].mass = particleMass;
        particles[particleIndex].smoothingLength = sqrt(particleMass / density);
        particles[particleIndex].pressure = pressure_local;
        particles[particleIndex].energy = pressure_local / (density * (gamma - 1.0));

        // Расчет скорости вращения
        double v_rot = sqrt(distance * sqrt(force_total_x * force_total_x + force_total_y * force_total_y));
        double ff = atan2(y, x);

        // Составление компонент скорости
        particles[particleIndex].velocityX = -v_rot * sin(ff);
        particles[particleIndex].velocityY = v_rot * cos(ff);

        ++particleIndex;
    }
}
