#pragma once
void saveParticlesToFile(double times, double Dt)
{
    std::string folderName = "data";
    if (!CreateDirectoryA(folderName.c_str(), NULL) && ERROR_ALREADY_EXISTS != GetLastError()) {
        std::cerr << "Error creating the folder!" << std::endl;
        return;
    }

    // Формирование имени файла
    std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;
    std::wstring path_wide = converter.from_bytes("data\\" + std::to_string(static_cast<int>(round(Dt * dt_out))) + ".bin");

    // Открытие файла в бинарном режиме
    std::ofstream fout(path_wide, std::ios::binary);
    if (!fout) {
        std::cerr << "Error opening the file!" << std::endl;
        return;
    }

    // Запись Dt
    fout.write(reinterpret_cast<const char*>(&Dt), sizeof(Dt));
    // Запись times
    fout.write(reinterpret_cast<const char*>(&times), sizeof(times));
    // Запись количества частиц
    int np_count = particles.size();
    fout.write(reinterpret_cast<const char*>(&np_count), sizeof(np_count));

    // Запись данных частиц
    for (const auto& particle : particles)
    {
        //if (particle.x > 3.0 || particle.x < -3.0 || particle.y > 3.0 || particle.y < -3.0) continue;
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
    std::ifstream file("configuration.csv");
    if (!file.is_open())
    {
        std::cerr << "Ошибка: Не удалось открыть файл configuration.csv!" << std::endl;
        return 0;
    }

    std::string line;
    std::getline(file, line);

    if (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string value;

        std::getline(iss, value, ',');
        numParticles = std::stoi(value);

        std::getline(iss, value, ',');
        radius = std::stod(value);

        std::getline(iss, value, ',');
        tau = std::stod(value);

        std::getline(iss, value, ',');
        gamma = std::stod(value);

        std::getline(iss, value, ',');
        alpha = std::stod(value);

        std::getline(iss, value, ',');
        beta = std::stod(value);

        std::getline(iss, value, ',');
        eps = std::stod(value);

        std::getline(iss, value, ',');
        t_end = std::stod(value);

        std::getline(iss, value, ',');
        t = std::stod(value);

        std::getline(iss, value, ',');
        Dt = std::stod(value);

        std::getline(iss, value, ',');
        system_Mass = std::stod(value);

        std::getline(iss, value, ',');
        scale = std::stod(value);

        std::getline(iss, value, ',');
        a_halo = std::stod(value);

        std::getline(iss, value, ',');
        b_halo = std::stod(value);

        std::getline(iss, value, ',');
        c_halo = std::stod(value);

        std::getline(iss, value, ',');
        scale_halo = std::stod(value);

        std::getline(iss, value, ',');
        Mah_d = std::stod(value);

        std::getline(iss, value, ',');
        r_core = std::stod(value);

        std::getline(iss, value, ',');
        rmax_hl = std::stod(value);

        std::getline(iss, value, ',');
        AMS = std::stod(value);
        
        std::getline(iss, value, ',');
        Ap = std::stod(value);

        std::getline(iss, value, ',');
        Lh = std::stod(value);

        std::getline(iss, value, ',');
        x_in = std::stod(value);

        std::getline(iss, value, ',');
        x_ex = std::stod(value);

        std::getline(iss, value, ',');
        y_in = std::stod(value);

        std::getline(iss, value, ',');
        y_ex = std::stod(value);

        std::cout << "Variables value:" << std::endl;
        std::cout << "numParticles = " << numParticles << std::endl;
        std::cout << "radius = " << radius << std::endl;
        std::cout << "tau = " << tau << std::endl;
        std::cout << "gamma = " << gamma << std::endl;
        std::cout << "alpha = " << alpha << std::endl;
        std::cout << "beta = " << beta << std::endl;
        std::cout << "eps = " << eps << std::endl;
        std::cout << "maximum = " << maximum << std::endl;
        std::cout << "t_end = " << t_end << std::endl;
        std::cout << "t = " << t << std::endl;
        std::cout << "Dt = " << Dt << std::endl;
        std::cout << "system_Mass = " << system_Mass << std::endl;
        std::cout << "a_halo = " << a_halo << std::endl;
        std::cout << "b_halo = " << b_halo << std::endl;
        std::cout << "c_halo = " << c_halo << std::endl;
        std::cout << "scale_halo = " << scale_halo << std::endl;
        std::cout << "mah_d = " << Mah_d << std::endl;
        std::cout << "r_core = " << r_core << std::endl;
        std::cout << "rmax_hl = " << rmax_hl << std::endl;
        std::cout << "AMS = " << AMS << std::endl;
        std::cout << "Ap = " << Ap << std::endl;
        std::cout << "Lh = " << Lh << std::endl;
        std::cout << "x_in = " << x_in << std::endl;
        std::cout << "x_ex = " << x_ex << std::endl;
        std::cout << "y_in = " << y_in << std::endl;
        std::cout << "y_ex = " << y_ex << std::endl;
    }
    else
    {
        std::cerr << "Ошибка: Файл configuration.csv не содержит данных во второй строке!" << std::endl;
        return 0;
    }

    file.close();
    return 1;
}

template <typename T>
void logVariable(const T& var) {
    std::lock_guard<std::mutex> lock(logMutex);
    std::ofstream logFile("C:\\Users\\SKENDLI\\Desktop\\diplom\\diplom_sph\\sph_method\\debug.log", std::ios::app);

    if (!logFile) {
        std::cerr << "Ошибка: не могу открыть debug.log" << std::endl;
        return;
    }

    logFile << var << std::endl;
    logFile.flush();
}

template <typename T>
void logVariable(const std::vector<T>& vec) {
    std::ofstream logFile("debug.log", std::ios::app);
    if (!logFile) {
        std::cerr << "Ошибка открытия файла debug.log" << std::endl;
        return;
    }

    std::ostringstream oss;
    oss << "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        oss << vec[i];
        if (i != vec.size() - 1) oss << ", ";
    }
    oss << "]";
    logFile << oss.str() << std::endl;
}
