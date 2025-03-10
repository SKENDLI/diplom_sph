#pragma once
class Grid
{
public:
    double cellSize;
    int gridSizeX, gridSizeY;
    double a, b;
    std::vector<std::vector<std::vector<Particle*>>> cells;

    Grid(double maxSmoothingLength, const std::vector<Particle>& particles, double minX, double maxX)
        : a(minX), b(maxX)
    {
        cellSize = 2 * maxSmoothingLength;

        gridSizeX = static_cast<int>(ceil((b - a) / cellSize));
        gridSizeY = static_cast<int>(ceil((b - a) / cellSize));

        cells.resize(gridSizeX);
        for (int i = 0; i < gridSizeX; ++i)
        {
            cells[i].resize(gridSizeY);
        }

        for (const auto& particle : particles)
        {
            if (particle.x >= a && particle.x <= b && particle.y >= a && particle.y <= b)
            {
                int cellX = static_cast<int>((particle.x - a) / cellSize);
                int cellY = static_cast<int>((particle.y - a) / cellSize);

                cellX = min(max(cellX, 0), gridSizeX - 1);
                cellY = min(max(cellY, 0), gridSizeY - 1);

                cells[cellX][cellY].push_back(const_cast<Particle*>(&particle));
            }
        }
    }

    ~Grid()
    {
        for (int i = 0; i < gridSizeX; ++i)
        {
            for (int j = 0; j < gridSizeY; ++j)
            {
                cells[i][j].clear();
            }
        }
        cells.clear();
    }

    const std::vector<Particle*>& getParticlesInCell(int cellX, int cellY) const
    {
        return cells[cellX][cellY];
    }

    int getGridSizeX() const
    {
        return gridSizeX;
    }

    int getGridSizeY() const
    {
        return gridSizeY;
    }
};
