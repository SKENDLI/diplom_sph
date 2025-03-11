#pragma once
void initializeGasCloud() {
    particles.resize(numParticles);
    double v_rot;
    double particleMass = system_Mass / numParticles;
    double scale_radius = radius / scale;
    int numRings = 1000;
    double ringWidth = radius / numRings;

    double RCORE1 = rmax_hl / r_core;
    double A_G = AMS / (RCORE1 - atan(RCORE1));
    B_rho = system_Mass / fun_mass(r_core, rmax_hl, 1000, M_PI, scale_radius);
    double rd_mah = 2.0;
    double B_p = B_rho / gamma * rd_mah * rd_mah * A_G / (Mah_d * Mah_d);

    A = 4.0 * M_PI * G * B_rho * scale_halo * scale_halo;

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
            double Fp_r = -gamma * B_p / B_rho * pow(computeRho(distance, scale_radius), (gamma - 2.0)) * dComputeRho(distance, scale_radius);

            double ksi = computeKsiForHalo(x, y, 0.0);
            double force_halo_x = computeForceGrav(ksi, x, A, a_halo);
            double force_halo_y = computeForceGrav(ksi, y, A, b_halo);
            double F_r = (force_halo_x * x + force_halo_y * y) / distance;
            //double Fuv_r = (sum_velocity_x[i] * x + sum_velocity_y[i] * y) / distance;
            double ff = atan2(y, x);
            v_rot = sqrt(-distance * (F_r + Fp_r));
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

    // Уменьшаем размер вектора до фактически использованных частиц
    if (particleIndex < numParticles) {
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
    }

    //ComputeForces(particles, sum_rho, sum_velocity_x, sum_velocity_y, sum_energy, radius);

    //for (int i = 0; i < particleIndex; i++) {

    //    double x = particles[i].x;
    //    double y = particles[i].y;
    //    double ff = atan2(y, x);
    //    double distance = sqrt(x * x + y * y);

    //    double ksi = computeKsiForHalo(x, y, 0.0);
    //    double force_halo_x = computeForceGrav(ksi, x, A, a_halo);
    //    double force_halo_y = computeForceGrav(ksi, y, A, b_halo);

    //    double F_r = (force_halo_x * x + force_halo_y * y) / distance;
    //    double Fuv_r = (sum_velocity_x[i] * x + sum_velocity_y[i] * y) / distance;

    //    v_rot = sqrt(-distance * (F_r + Fuv_r));
    //    particles[i].velocityX = -v_rot * sin(ff);
    //    particles[i].velocityY = v_rot * cos(ff);
    //}
}