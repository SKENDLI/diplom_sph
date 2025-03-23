#pragma once
class Grid {
public:
    Grid(double minX, double maxX, double minY, double maxY, double cellSize)
        : minX(minX), maxX(maxX), minY(minY), maxY(maxY), cellSize(cellSize) {
        numCellsX = static_cast<int>(std::ceil((maxX - minX) / cellSize));
        numCellsY = static_cast<int>(std::ceil((maxY - minY) / cellSize));
        cells.resize(numCellsX * numCellsY);
    }

    void build(const std::vector<Particle>& particles) {
        for (auto& cell : cells) {
            cell.clear();
        }
        for (size_t i = 0; i < particles.size(); ++i) {
            int cellX = static_cast<int>(std::floor((particles[i].x - minX) / cellSize));
            int cellY = static_cast<int>(std::floor((particles[i].y - minY) / cellSize));
            if (cellX >= 0 && cellX < numCellsX && cellY >= 0 && cellY < numCellsY) {
                int cellIndex = cellY * numCellsX + cellX;
                cells[cellIndex].push_back(i);
            }
        }
    }

    std::vector<size_t> getNeighbors(size_t particleIndex, const std::vector<Particle>& particles, double searchRadius) {
        std::vector<size_t> neighbors;
        const Particle& p = particles[particleIndex];
        // Вычисляем координаты ячейки для текущей частицы
        int cellX = static_cast<int>(std::floor((p.x - minX) / cellSize));
        int cellY = static_cast<int>(std::floor((p.y - minY) / cellSize));

        // Определяем границы поиска: текущая ячейка ±1
        int minCellX = max(0, cellX - 1);
        int maxCellX = min(numCellsX - 1, cellX + 1);
        int minCellY = max(0, cellY - 1);
        int maxCellY = min(numCellsY - 1, cellY + 1);

        // Проходим по всем соседним ячейкам (включая текущую)
        for (int cy = minCellY; cy <= maxCellY; ++cy) {
            for (int cx = minCellX; cx <= maxCellX; ++cx) {
                int cellIndex = cy * numCellsX + cx;
                // Проверяем все частицы в текущей ячейке
                for (size_t idx : cells[cellIndex]) {
                    if (idx != particleIndex) { // Исключаем саму частицу
                        const Particle& q = particles[idx];
                        double dx = p.x - q.x;
                        double dy = p.y - q.y;
                        double distanceSquared = dx * dx + dy * dy;
                        // Добавляем только те частицы, что находятся в пределах searchRadius
                        if (distanceSquared < searchRadius * searchRadius) {
                            neighbors.push_back(idx);
                        }
                    }
                }
            }
        }
        return neighbors;
    }

    void printCells() const {
        std::cout << "CELLS CONTENT:\n";
        for (int y = 0; y < numCellsY; ++y) {
            for (int x = 0; x < numCellsX; ++x) {
                int cellIndex = y * numCellsX + x;
                if (!cells[cellIndex].empty()) {
                    std::cout << "Cell (" << x << ", " << y << ") [cellCount " << cells[cellIndex].size() << "]: ";
                    std::cout << "\n";
                }
            }
        }
        std::cout << "END CELL CONTENT.\n";
    }

private:
    double minX, maxX, minY, maxY, cellSize;
    int numCellsX, numCellsY;
    std::vector<std::vector<size_t>> cells;
};
