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
        int cellX = static_cast<int>(std::floor((p.x - minX) / cellSize));
        int cellY = static_cast<int>(std::floor((p.y - minY) / cellSize));

        int minCellX = max(0, cellX - 1);
        int maxCellX = min(numCellsX - 1, cellX + 1);
        int minCellY = max(0, cellY - 1);
        int maxCellY = min(numCellsY - 1, cellY + 1);

        for (int cy = minCellY; cy <= maxCellY; ++cy) {
            for (int cx = minCellX; cx <= maxCellX; ++cx) {
                int cellIndex = cy * numCellsX + cx;
                for (size_t idx : cells[cellIndex]) {
                    if (idx != particleIndex) {
                        const Particle& q = particles[idx];
                        double dx = p.x - q.x;
                        double dy = p.y - q.y;
                        double distanceSquared = dx * dx + dy * dy;
                        if (distanceSquared < searchRadius * searchRadius) {
                            neighbors.push_back(idx);
                        }
                    }
                }
            }
        }
        return neighbors;
    }

private:
    double minX, maxX, minY, maxY, cellSize;
    int numCellsX, numCellsY;
    std::vector<std::vector<size_t>> cells;
};
