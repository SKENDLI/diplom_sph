#pragma once
using namespace std;
void saveParticlesToFile(double times, double Dt)
{
    string folderName = "data";
    if (!CreateDirectoryA(folderName.c_str(), NULL) && ERROR_ALREADY_EXISTS != GetLastError()) {
        cerr << "Error creating the folder!" << endl;
        return;
    }

    // Формирование имени файла
    wstring_convert<codecvt_utf8_utf16<wchar_t>> converter;
    wstring path_wide = converter.from_bytes("data\\" + to_string(static_cast<int>(round(Dt * dt_out))) + ".bin");

    // Открытие файла в бинарном режиме
    ofstream fout(path_wide, ios::binary);
    if (!fout) {
        cerr << "Error opening the file!" << endl;
        return;
    }

    fout.write(reinterpret_cast<const char*>(&Dt), sizeof(Dt));
    fout.write(reinterpret_cast<const char*>(&times), sizeof(times));
    int np_count = particles.size();
    fout.write(reinterpret_cast<const char*>(&np_count), sizeof(np_count));

    // Запись данных частиц
    for (const auto& particle : particles)
    {
        fout.write(reinterpret_cast<const char*>(&particle.x), sizeof(particle.x));
        fout.write(reinterpret_cast<const char*>(&particle.y), sizeof(particle.y));
        fout.write(reinterpret_cast<const char*>(&particle.density), sizeof(particle.density));
        fout.write(reinterpret_cast<const char*>(&particle.mass), sizeof(particle.mass));
        fout.write(reinterpret_cast<const char*>(&particle.pressure), sizeof(particle.pressure));
        fout.write(reinterpret_cast<const char*>(&particle.velocityX), sizeof(particle.velocityX));
        fout.write(reinterpret_cast<const char*>(&particle.velocityY), sizeof(particle.velocityY));
        fout.write(reinterpret_cast<const char*>(&particle.energy), sizeof(particle.energy));
        fout.write(reinterpret_cast<const char*>(&particle.smoothingLength), sizeof(particle.smoothingLength));
        fout.write(reinterpret_cast<const char*>(&particle.neighbors), sizeof(particle.neighbors));
    }

    fout.close();
}

int configuration()
{
    ifstream file("configuration.csv");
    if (!file.is_open())
    {
        cerr << "Ошибка: Не удалось открыть файл configuration.csv!" << endl;
        return 0;
    }

    string line;
    getline(file, line);

    if (getline(file, line))
    {
        istringstream iss(line);
        string value;

        getline(iss, value, ',');
        numParticles = stoi(value);

        getline(iss, value, ',');
        radius = stod(value);

        getline(iss, value, ',');
        tau = stod(value);

        getline(iss, value, ',');
        gamma = stod(value);

        getline(iss, value, ',');
        alpha = stod(value);

        getline(iss, value, ',');
        beta = stod(value);

        getline(iss, value, ',');
        eps = stod(value);

        getline(iss, value, ',');
        t_end = stod(value);

        getline(iss, value, ',');
        t = stod(value);

        getline(iss, value, ',');
        Dt = stod(value);

        getline(iss, value, ',');
        system_Mass = stod(value);

        getline(iss, value, ',');
        scale = stod(value);

        getline(iss, value, ',');
        scale_halo = stod(value);

        getline(iss, value, ',');
        Mah_d = stod(value);

        getline(iss, value, ',');
        r_core = stod(value);

        getline(iss, value, ',');
        rmax_hl = stod(value);

        getline(iss, value, ',');
        AMS = stod(value);
        
        getline(iss, value, ',');
        Ap = stod(value);

        getline(iss, value, ',');
        Lh = stod(value);

        getline(iss, value, ',');
        x_in = stod(value);

        getline(iss, value, ',');
        x_ex = stod(value);

        getline(iss, value, ',');
        y_in = stod(value);

        getline(iss, value, ',');
        y_ex = stod(value);

        getline(iss, value, ',');
        Omega_h = stod(value);

        getline(iss, value, ',');
        a_h = stod(value);

        getline(iss, value, ',');
        tau_h = stod(value);

        cout << "Variables value:" << endl;
        cout << "numParticles = " << numParticles << endl;
        cout << "radius = " << radius << endl;
        cout << "tau = " << tau << endl;
        cout << "gamma = " << gamma << endl;
        cout << "alpha = " << alpha << endl;
        cout << "beta = " << beta << endl;
        cout << "eps = " << eps << endl;
        cout << "maximum = " << maximum << endl;
        cout << "t_end = " << t_end << endl;
        cout << "t = " << t << endl;
        cout << "Dt = " << Dt << endl;
        cout << "system_Mass = " << system_Mass << endl;
        cout << "scale_halo = " << scale_halo << endl;
        cout << "mah_d = " << Mah_d << endl;
        cout << "r_core = " << r_core << endl;
        cout << "rmax_hl = " << rmax_hl << endl;
        cout << "AMS = " << AMS << endl;
        cout << "Ap = " << Ap << endl;
        cout << "Lh = " << Lh << endl;
        cout << "x_in = " << x_in << endl;
        cout << "x_ex = " << x_ex << endl;
        cout << "y_in = " << y_in << endl;
        cout << "y_ex = " << y_ex << endl;
        cout << "Omega_h = " << Omega_h << endl;
        cout << "a_h = " << a_h << endl;
        cout << "tau_h = " << tau_h << endl;
    }
    else
    {
        cerr << "Ошибка: Файл configuration.csv не содержит данных во второй строке!" << endl;
        return 0;
    }

    file.close();
    return 1;
}

template <typename T>
void logVariable(const T& var) {
    lock_guard<mutex> lock(logMutex);
    ofstream logFile("C:\\Users\\SKENDLI\\Desktop\\diplom\\diplom_sph\\sph_method\\debug.log", ios::app);
    logFile << var << endl;
    logFile.flush();
}

template <typename T>
void logVariable(const vector<T>& vec) {
    ofstream logFile("debug.log", ios::app);
    ostringstream oss;
    oss << "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        oss << vec[i];
        if (i != vec.size() - 1) oss << ", ";
    }
    oss << "]";
    logFile << oss.str() << endl;
}

void save_laws_cons(double t, double& next_print_time, double initial_energy, vector<Particle>& particles) {
    if (t >= next_print_time) {
        double current_energy = 0.0;
        law_cons(current_energy, particles);

        ofstream file("laws_cons.txt", ios::app);
        if (file.is_open()) {
            file << "Initial energy: " << initial_energy
                << " Current energy: " << current_energy
                << " Difference: " << abs(initial_energy - current_energy)
                << " Time: " << next_print_time
                << endl;
            file.close();
        }
        else {
            cerr << "Error opening file laws_cons.txt" << endl;
        }

        cout << "Initial energy: " << initial_energy
            << " Current energy: " << current_energy
            << " Difference: " << abs(initial_energy - current_energy)
            << " Time: " << next_print_time
            << endl;

        next_print_time += 0.1;
    }
}
