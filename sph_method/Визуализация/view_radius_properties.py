import struct
import numpy as np
import matplotlib.pyplot as plt

def read_particle_data(file_path):
    with open(file_path, "rb") as file:
        # Считывание Dt и times
        Dt = struct.unpack('d', file.read(8))[0]
        times = struct.unpack('d', file.read(8))[0]
        # Считывание количества частиц
        num_particles = struct.unpack('i', file.read(4))[0]
        particles = []
        for _ in range(num_particles):
            particle = {
                "x": struct.unpack('d', file.read(8))[0],
                "y": struct.unpack('d', file.read(8))[0],
                "density": struct.unpack('d', file.read(8))[0],
                "mass": struct.unpack('d', file.read(8))[0],
                "pressure": struct.unpack('d', file.read(8))[0],
                "velocityX": struct.unpack('d', file.read(8))[0],
                "velocityY": struct.unpack('d', file.read(8))[0],
                "energy": struct.unpack('d', file.read(8))[0],
                "smoothingLength": struct.unpack('d', file.read(8))[0],
                "neighbors": struct.unpack('i', file.read(4))[0],
            }
            particle["radius"] = np.sqrt(particle["x"]**2 + particle["y"]**2)  # Добавляем радиус
            particle["velocity"] = np.sqrt(particle["velocityX"]**2 + particle["velocityY"]**2)  # Модуль скорости
            particles.append(particle)
        return particles

def visualize_particle_data(particles, x_key, y_key):
    """
    Визуализирует данные частиц.
    
    :param particles: Список данных частиц
    :param x_key: Ключ для данных по оси X
    :param y_key: Ключ для данных по оси Y
    """
    x_data = [particle[x_key] for particle in particles]
    y_data = [particle[y_key] for particle in particles]
    
    # Для симметрии отразим радиусы, если x_key — "radius"
    if x_key == "radius":
        x_data = np.concatenate([[-r for r in x_data], x_data])
        y_data = np.concatenate([y_data, y_data])
    
    plt.figure(figsize=(8, 6))
    plt.scatter(x_data, y_data, s=1, label=f'Распределение {y_key} по радиусу диска', color='blue')
    plt.xlabel(x_key, fontsize=14)
    plt.ylabel(y_key, fontsize=14)
    plt.title(f'Распределение {y_key} по радиусу диска', fontsize=16)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)
    plt.show()

def visualize_velocity_by_radius(particles):
    """
    Визуализирует зависимость скорости от радиуса.
    
    :param particles: Список данных частиц
    """
    radius_data = [particle["radius"] for particle in particles]
    velocity_data = [particle["velocity"] for particle in particles]
    
    plt.figure(figsize=(8, 6))
    plt.scatter(radius_data, velocity_data, s=1, label='Скорость по радиусу', color='blue')
    plt.xlabel('Радиус', fontsize=14)
    plt.ylabel('Модуль скорости', fontsize=14)
    plt.title('Зависимость скорости от радиуса частиц', fontsize=16)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)
    plt.show()

# Прочитаем данные из файла
particles = read_particle_data("C:/Users/SKENDLI/Desktop/diplom/diplom_sph/sph_method/data/0.bin")

# Выбор данных для построения графика
x_column = "radius"  # Используем радиус как X
y_column = "pressure"  # Можно выбрать любую характеристику, например, скорость по X

# Визуализируем выборочные данные
#visualize_particle_data(particles, x_column, y_column)

# Визуализируем зависимость скорости от радиуса
visualize_velocity_by_radius(particles)
