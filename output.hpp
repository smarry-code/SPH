void output(std::vector<Particle>& Elem, size_t ntotal, double time, size_t itimestep, double dt)
{
#ifdef _WIN32
    SetConsoleOutputCP(65001); // Активация UTF-8 в Windows
#endif
    // Создание имени файла с временной меткой
    std::ostringstream filename;
    filename << "sph_output_" << std::setw(5) << std::setfill('0') << itimestep << ".txt";

    // Открытие файла для записи
    std::ofstream file(filename.str());
    if (!file.is_open())
    {
        throw std::runtime_error("Не удалось открыть файл для вывода данных");
    }

    // Запись заголовка
    /*file << "SPH Simulation Output\n";
    file << "Timestep: " << itimestep << "\n";
    file << "Time: " << time << "\n";
    file << "Total particles: " << ntotal << "\n";
    file << "-----------------------------------------\n";

    // Запись заголовков колонок
    file << "ID Type Mass";
    for (int d = 0; d < dim; ++d)
    {
        file << " X" << d + 1;
    }
    for (int d = 0; d < dim; ++d)
    {
        file << " V" << d + 1;
    }
    file << " Density Pressure Energy SoundSpeed HSML\n";      */

    // Запись данных о частицах
    file << "N" << " " << "Nmat" << " " << "M" << " ";
    if (dim == 1)
    {
        file << "x" << " ";
    }
    else if (dim == 2)
    {
         file << "x" << " " << "y" << " ";
    }
    else if (dim == 3)
    {
        file << "x" << " " << "y" << " " << "z" << " ";
    }
    if (dim == 1)
    {
        file << "vx" << " ";
    }
    else if (dim == 2)
    {
         file << "vx" << " " << "vy" << " ";
    }
    else if (dim == 3)
    {
        file << "vx" << " " << "vy" << " " << "vz" << " ";
    }
    file << "rho" << " " << "p" << " " << "u" << " " << "c" << " " << "sigma_i" << " " << "e_pl" << " " << "D" << " " << "hsml" << " "
         << "sxx" << " " << "syy" << " " << "szz" << " " << "sxy" << " " << "syz" << " " << "sxz" << "\n";
    for (size_t i = 0; i < ntotal; ++i)
    {
        // Запись основных данных
        file << i + 1;  // ID частицы (в Fortran индексы начинаются с 1)
        Elem[i].vuvod(file);
    }

    file.close();
    //--------------------------
    std::ostringstream filename1;
    filename1 << "sph_gp" << ".csv";
    std::ofstream file1(filename1.str(), std::ios::app);
    for (size_t i = 0; i < ntotal; ++i)
    {
        file1 << time << " ";
        for (size_t d = 0; d < dim; ++d)
        {
            file1 << Elem[i].x[d] << " ";
        }
        for (size_t d = 0; d < dim; ++d)
        {
            file1 << Elem[i].vx[d] << " ";
        }
        file1 << "\n";
    }
    file1.close();
    //--------------------------
    // Вывод информации о сохранении данных
    std::wcout << "Данные сохранены в файл: " << filename.str().c_str() << "\n";
    std::wcout << "Шаг по времени: " << dt << " Время расчета: " << time << "\n";
    std::wcout << "Текущий цикл: " << itimestep << "\n";
}
