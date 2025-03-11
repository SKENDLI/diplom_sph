#pragma once
class GasDiskGrid {
    // Конфигурация сетки
    double domain_radius;  // Радиус моделируемой области
    double cell_size;      // Размер ячейки (обычно ~max_h)

    // Границы сетки
    double x_min, x_max;
    double y_min, y_max;

    // Топология сетки
    int num_cells_x;
    int num_cells_y;

    // Данные
    std::vector<std::vector<std::vector<Particle*>>> cells;

public:
    GasDiskGrid(double disk_radius, double max_h)
        : domain_radius(disk_radius),
        cell_size(max_h * 1.2)  // Запас 20% для перекрытия
    {
        // Инициализация границ
        x_min = -domain_radius - cell_size * 2;
        x_max = domain_radius + cell_size * 2;
        y_min = -domain_radius - cell_size * 2;
        y_max = domain_radius + cell_size * 2;

        // Расчёт числа ячеек
        num_cells_x = int((x_max - x_min) / cell_size) + 1;
        num_cells_y = int((y_max - y_min) / cell_size) + 1;

        // Выделение памяти
        cells.resize(num_cells_x);
        for (auto& col : cells) {
            col.resize(num_cells_y);
        }
    }

    // Обновление сетки с новыми позициями частиц
    void rebuild(const std::vector<Particle>& particles) {
        // Очистка предыдущих данных
        for (auto& col : cells) {
            for (auto& cell : col) {
                cell.clear();
            }
        }

        // Распределение частиц
        for (const auto& p : particles) {
            // Только частицы внутри рабочей области
            if (p.distanceFromCenter > domain_radius) continue;

            // Вычисление индексов ячейки
            int i = int((p.x - x_min) / cell_size);
            int j = int((p.y - y_min) / cell_size);

            // Ограничение индексов
            if (i < 0) i = 0;
            if (j < 0) j = 0;
            if (i >= num_cells_x) i = num_cells_x - 1;
            if (j >= num_cells_y) j = num_cells_y - 1;

            cells[i][j].push_back(const_cast<Particle*>(&p));
        }
    }

    // Получение кандидатов в соседи
    void get_neighbors(const Particle* p, std::vector<Particle*>& neighbors) const {
        neighbors.clear();

        // Определение центральной ячейки
        const int ic = int((p->x - x_min) / cell_size);
        const int jc = int((p->y - y_min) / cell_size);

        // Поиск в окрестности 5x5 ячеек
        for (int di = -2; di <= 2; ++di) {
            for (int dj = -2; dj <= 2; ++dj) {
                const int i = ic + di;
                const int j = jc + dj;

                // Проверка границ сетки
                if (i < 0 || i >= num_cells_x) continue;
                if (j < 0 || j >= num_cells_y) continue;

                // Фильтрация ячеек вне диска
                if (!cell_intersects_disk(i, j)) continue;

                // Сбор частиц
                for (auto* candidate : cells[i][j]) {
                    // Проверка радиального расстояния
                    const double dx = p->x - candidate->x;
                    const double dy = p->y - candidate->y;
                    if (dx * dx + dy * dy < (2.0 * p->smoothingLength) * (2.0 * p->smoothingLength)) {
                        if (candidate != p) {
                            neighbors.push_back(candidate);
                        }
                    }
                }
            }
        }
    }

private:
    // Проверка пересечения ячейки с диском
    bool cell_intersects_disk(int i, int j) const {
        // Координаты углов ячейки
        const double x1 = x_min + i * cell_size;
        const double y1 = y_min + j * cell_size;
        const double x2 = x1 + cell_size;
        const double y2 = y1 + cell_size;

        const double cx = (x1 <= 0 && x2 >= 0) ? 0 : (fabs(x1) < fabs(x2) ? x1 : x2);
        const double cy = (y1 <= 0 && y2 >= 0) ? 0 : (fabs(y1) < fabs(y2) ? y1 : y2);

        // Расстояние ближайшей точки
        const double r_min = sqrt(cx * cx + cy * cy);

        // Дальняя точка
        const double fx = (fabs(x1) > fabs(x2)) ? x1 : x2;
        const double fy = (fabs(y1) > fabs(y2)) ? y1 : y2;
        const double r_max = sqrt(fx * fx + fy * fy);

        return !(r_min > domain_radius || r_max < (domain_radius - cell_size * 1.414));
    }
};