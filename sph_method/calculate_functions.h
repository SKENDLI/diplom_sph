#pragma once
#include "globals.h"
#include "grid.h"
using namespace std;
void computeForcesGravCool(vector<Particle>& particles, vector<double>& force_halo_x, vector<double>& force_halo_y) {
    for (size_t i = 0; i < particles.size(); ++i) {
        double xi = particles[i].x;
        double yi = particles[i].y;
        double ri = hypot(xi, yi);

        double at;
        if (t <= tau_h) {
            at = 1.0 - (1.0 - a_h) * t / tau_h;
        }
        else {
            at = a_h;
        }
        double at2 = at * at;

        double cw = cos(Omega_h * t);
        double sw = sin(Omega_h * t);
        double cw2 = cw * cw;
        double sw2 = sw * sw;

        double ksi = sqrt(
            xi * xi * (cw2 + sw2 / at2) +
            yi * yi * (sw2 + cw2 / at2) +
            2.0 * xi * yi * cw * sw * (1.0 - 1.0 / at2)
        );

        double ksicore = ksi / r_core;
        double forcehalo = -A_G * (ksicore - atan(ksicore)) / (ksi * ksi * ksi);

        force_halo_x[i] = forcehalo * (xi * (cw2 + sw2 / at2) + yi * cw * sw * (1.0 - 1.0 / at2));
        force_halo_y[i] = forcehalo * (yi * (sw2 + cw2 / at2) + xi * cw * sw * (1.0 - 1.0 / at2));
    }
}

double computeRho(double Lf, double rr)
{
    double cosh_val = cosh(rr / Lf);
    return 1.0 / (cosh_val);
}

double dComputeRho(double distance, double scale_radius)
{
    return 1.0 / scale_radius * sinh(distance / scale_radius) / pow(cosh(distance / scale_radius), 2.0);
}

double fun_mass(double r1, double r2, int Nmm, double pi2, double Lf) {
    double fun_mass = 0.0;
    double dr = (r2 - r1) / static_cast<double>(Nmm);

    for (int k = 1; k <= Nmm; ++k) {
        double rr = r1 + dr * static_cast<double>(k - 1);
        fun_mass += 0.5 * (rr * computeRho(Lf, rr) + (rr + dr) * computeRho(Lf, rr + dr)) * dr;
    }

    fun_mass *= pi2;
    return fun_mass;
}


double mass_difference(double rr, double rrk, double rrk_1, double target_mass, double B_rho, double m_p, double Lf, double pi2) {
    double current_mass = B_rho * fun_mass(rrk, rr, 1000, pi2, Lf);
    return current_mass - target_mass;
}

double find_radius(double rrk, double rrk_1, double target_mass, double B_rho, double m_p, double Lf, double pi2) {
    double rr = rrk + (rrk - rrk_1); // Начальное приближение
    double drr = 1e-5;                // Шаг для численной производной
    double tolerance = 1e-6;          // Точность
    int max_iter = 100;               // Максимум итераций
    for (int iter = 0; iter < max_iter; ++iter) {
        double ff_rr = mass_difference(rr, rrk, rrk_1, target_mass, B_rho, m_p, Lf, pi2);
        double ff_rr1 = mass_difference(rr + drr, rrk, rrk_1, target_mass, B_rho, m_p, Lf, pi2);
        double dff = drr * ff_rr / (ff_rr1 - ff_rr); // Защита от деления на 0
        rr -= dff;
        if (abs(ff_rr) < tolerance && abs(dff) < tolerance) {
            break;
        }
    }
    return rr;
}

void ComputeForces(vector<Particle>& particles) {
    if (particles.empty()) return;

    size_t Np = particles.size();
    Fu_p.assign(Np, 0.0);
    Fv_p.assign(Np, 0.0);
    Fe_p.assign(Np, 0.0);
    Rhot.assign(Np, 0.0);
    neighborCount.assign(Np, 0);

    // Обновление hp_min
    double hp_min = hp_max;
    for (size_t i = 0; i < Np; ++i) {
        if (ht_p[i] < hp_min) hp_min = ht_p[i];
    }
    h = hp_min;
    h2 = 2.0 * h;

    // Сортировка частиц по ячейкам
    vector<vector<vector<int>>> box(Ngx, vector<vector<int>>(Ngy));
    vector<vector<int>> Nbox(Ngx, vector<int>(Ngy, 0));
    vector<vector<double>> hmaxbox(Ngx, vector<double>(Ngy, 0.0));
    int Np_box_max = 0;

    // Первый проход: определяем максимальное количество частиц в ячейке
    for (size_t i = 0; i < Np; ++i) {
        int xbox = static_cast<int>((particles[i].x - x_in) / h2 + 1.0);
        int ybox = static_cast<int>((particles[i].y - y_in) / h2 + 1.0);
        if (xbox >= 1 && xbox <= Ngx && ybox >= 1 && ybox <= Ngy) {
            xbox--; ybox--;
            Nbox[xbox][ybox]++;
            if (Nbox[xbox][ybox] > Np_box_max) Np_box_max = Nbox[xbox][ybox];
        }
    }

    // Выделяем память под box
    for (int x = 0; x < Ngx; ++x) {
        for (int y = 0; y < Ngy; ++y) {
            box[x][y].resize(Np_box_max, 0);
        }
    }

    // Второй проход: распределяем частицы по ячейкам
    Nbox.assign(Ngx, vector<int>(Ngy, 0));
    for (size_t i = 0; i < Np; ++i) {
        int xbox = static_cast<int>((particles[i].x - x_in) / h2 + 1.0);
        int ybox = static_cast<int>((particles[i].y - y_in) / h2 + 1.0);
        if (xbox >= 1 && xbox <= Ngx && ybox >= 1 && ybox <= Ngy) {
            xbox--; ybox--;
            Nbox[xbox][ybox]++;
            box[xbox][ybox][Nbox[xbox][ybox] - 1] = i;
            if (ht_p[i] > hmaxbox[xbox][ybox]) hmaxbox[xbox][ybox] = ht_p[i];
        }
    }


    // Вычисление плотности (Rhot)
    for (size_t i = 0; i < Np; ++i) {
        int xbox = static_cast<int>((particles[i].x - x_in) / h2 + 1.0);
        int ybox = static_cast<int>((particles[i].y - y_in) / h2 + 1.0);
        if (xbox >= 1 && xbox <= Ngx && ybox >= 1 && ybox <= Ngy) {

            xbox--; ybox--;
            double xi = particles[i].x;
            double yi = particles[i].y;
            double hpi = ht_p[i];

            int kbox = static_cast<int>(hmaxbox[xbox][ybox] / h + 1.0);

            for (int ybox_1 = ybox - kbox; ybox_1 <= ybox + kbox; ++ybox_1) {
                if (ybox_1 >= 0 && ybox_1 < Ngy) {
                    for (int xbox_1 = xbox - kbox; xbox_1 <= xbox + kbox; ++xbox_1) {
                        if (xbox_1 >= 0 && xbox_1 < Ngx) {
                            for (int l = 0; l < Nbox[xbox_1][ybox_1]; ++l) {
                                size_t j = box[xbox_1][ybox_1][l];
                                if (j >= i) {
                                    double s = sqrt(pow(xi - particles[j].x, 2) + pow(yi - particles[j].y, 2));
                                    double hij = 0.5 * (hpi + ht_p[j]);
                                    if (s <= 2.0 * hij) {
                                        double mpW = particles[i].mass * W(s, hij);
                                        Rhot[i] += mpW;
                                        if (j != i) Rhot[j] += mpW;
                                        neighborCount[i]++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    for (size_t i = 0; i < Np; ++i) {
        ht_p[i] = particles[i].smoothingLength * sqrt(Rho[i] / Rhot[i]);
        if (ht_p[i] > hp_max) ht_p[i] = hp_max;
    }

    // Вычисление сил и временного шага
    double dt_1 = 1.0;
    double gamma1 = gamma - 1.0;

    for (size_t i = 0; i < Np; ++i) {
        int xbox = static_cast<int>((particles[i].x - x_in) / h2 + 1.0);
        int ybox = static_cast<int>((particles[i].y - y_in) / h2 + 1.0);
        maximum = 0.0;

        if (xbox >= 1 && xbox <= Ngx && ybox >= 1 && ybox <= Ngy) {
            xbox--; ybox--;
            double xi = particles[i].x;
            double yi = particles[i].y;
            double ui = particles[i].velocityX;
            double vi = particles[i].velocityY;
            double Rhoi = Rhot[i];
            double Ppi = gamma1 * Rhoi * particles[i].energy;
            double hpi = ht_p[i];
            int kbox = static_cast<int>(hmaxbox[xbox][ybox] / h + 1.0);

            for (int ybox_1 = ybox - kbox; ybox_1 <= ybox + kbox; ++ybox_1) {
                if (ybox_1 >= 0 && ybox_1 < Ngy) {
                    for (int xbox_1 = xbox - kbox; xbox_1 <= xbox + kbox; ++xbox_1) {
                        if (xbox_1 >= 0 && xbox_1 < Ngx) {
                            for (int l = 0; l < Nbox[xbox_1][ybox_1]; ++l) {
                                size_t j = box[xbox_1][ybox_1][l];
                                if (j > i) {
                                    double xj = particles[j].x;
                                    double yj = particles[j].y;
                                    double xij = xi - xj;
                                    double yij = yi - yj;
                                    double s = sqrt(xij * xij + yij * yij);
                                    double hpj = ht_p[j];
                                    double hij = 0.5 * (hpi + hpj);
                                    if (s <= 2.0 * hij) {
                                        double uj = particles[j].velocityX;
                                        double vj = particles[j].velocityY;
                                        double uij = ui - uj;
                                        double vij = vi - vj;
                                        double urdot = uij * xij + vij * yij;
                                        double Rhoj = Rhot[j];
                                        double Ppj = gamma1 * Rhoj * particles[j].energy;

                                        double Rho_ij = 0.5 * (Rhoi + Rhoj);
                                        double cij = 0.5 * (sqrt(gamma * Ppi / Rhoi) + sqrt(gamma * Ppj / Rhoj));
                                        double mu_ij = mu(s, urdot, hij);
                                        if (maximum < mu_ij) {
                                            maximum = mu_ij;
                                        }

                                        double Aij = (beta * mu_ij - alpha * cij) * mu_ij / Rho_ij;
                                        double Dij = Ppi / (Rhoi * Rhoi) + Ppj / (Rhoj * Rhoj) + Aij;
                                        double dWx = gradW(s, xij, hij);
                                        double dWy = gradW(s, yij, hij);
                                        double dWxyuv = dWx * (ui - uj) + dWy * (vi - vj);

                                        Fu_p[i] -= particles[i].mass * Dij * dWx;
                                        Fv_p[i] -= particles[i].mass * Dij * dWy;
                                        Fe_p[i] += 0.5 * particles[i].mass * Dij * dWxyuv;

                                        Fu_p[j] += particles[i].mass * Dij * dWx;
                                        Fv_p[j] += particles[i].mass * Dij * dWy;
                                        Fe_p[j] += 0.5 * particles[i].mass * Dij * dWxyuv;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // Временной шаг
            double ci = sqrt(gamma * Ppi / Rhoi);
            double a1 = hpi / (ci * (1.0 + 1.2 * alpha) + 1.2 * beta * maximum);
            if (dt_1 > a1) dt_1 = a1;
        }
    }
    dtn = dt_1;
}

double W(double s, double hij) {
    double q = s / hij;
    double result;

    if (q <= 1.0) {
        result = 1.0 + (0.75 * q - 1.5) * q * q;
    }
    else if (q <= 2.0 && q > 1.0) {
        result = 0.25 * pow(2.0 - q, 3);
    }
    else {
        result = 0.0;
    }

    const double A_W = 10.0 / (7.0 * M_PI);
    result = A_W * result / (hij * hij);

    return result;
}

double gradW(double s, double rij, double hij) {
    double q = s / hij;
    double q1 = 2.0 - q;
    double result;

    if (q <= 1.0) {
        result = (-3.0 + 2.25 * q) * q;
    }
    else if (q <= 2.0 && q > 1.0) {
        result = -0.75 * q1 * q1;
    }
    else {
        result = 0.0;
    }

    const double A_W = 10.0 / (7.0 * M_PI);
    result = A_W * result * (rij / s) / (hij * hij * hij);

    return result;
}

double mu(double s, double urdot, double hij) {
    double q = s / hij;
    double result;
    if (urdot < 0.0) {
        double denominator = hij * (q * q + eps * eps);

        return result = urdot / denominator;
    }
    return result = 0.0;
}

double SoundSpeed(double p_i, double rho_i)
{
    return sqrt(gamma * p_i / rho_i);
}

inline double particleDistance(const Particle& a, const Particle& b, double softening = 1e-4) {
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    return sqrt(dx * dx + dy * dy + softening * softening);
}

void law_cons(double& Energy, const std::vector<Particle>& parts) {
    double Psi;
    Energy = 0.0;

    for (int i = 0; i < numParticles; ++i) {
        Psi = -1000.0 + 0.5 * A_G * (pow(parts[i].x, 2) + pow(parts[i].y, 2));
        Energy +=
            parts[i].mass * 0.5 * (pow(parts[i].velocityX, 2)
            + pow(parts[i].velocityY, 2))
            + parts[i].mass * parts[i].energy
            + Psi * parts[i].mass;
    }
}