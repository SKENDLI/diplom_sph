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
    Визуализирует зависимость скорости от радиуса с усреднением (только положительные радиусы).
    
    :param particles: Список данных частиц
    """
    # Извлечём радиусы и скорости
    radius_data = [particle["radius"] for particle in particles]
    velocity_data = [particle["velocity"] for particle in particles]

    # Уберём выбросы (например, скорости > 10 или радиусы < 0)
    filtered_data = [(r, v) for r, v in zip(radius_data, velocity_data) if r >= 0 and v < 10]
    if not filtered_data:
        print("Нет данных для визуализации после фильтрации!")
        return
    radius_data, velocity_data = zip(*filtered_data)

    # Усреднение скорости по радиусам
    # Создаём словарь для группировки
    radius_velocity_dict = {}
    for r, v in zip(radius_data, velocity_data):
        # Округляем радиус до 2 знаков для группировки
        r_rounded = round(r, 2)
        if r_rounded not in radius_velocity_dict:
            radius_velocity_dict[r_rounded] = []
        radius_velocity_dict[r_rounded].append(v)

    # Вычисляем средние значения
    avg_radius = []
    avg_velocity = []
    for r, velocities in sorted(radius_velocity_dict.items()):
        avg_radius.append(r)
        avg_velocity.append(np.mean(velocities))

    # Оставляем только положительные радиусы (убираем симметрию)
    x_data = avg_radius
    y_data = avg_velocity

    plt.figure(figsize=(8, 6))
    plt.plot(x_data, y_data, label='Скорость по радиусу', color='blue', linestyle='--')  # Используем plot для непрерывной линии
    plt.xlabel('Радиус', fontsize=14)
    plt.ylabel('Модуль скорости', fontsize=14)
    plt.title('Зависимость скорости от радиуса частиц', fontsize=16)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)
    plt.legend()
    plt.show()

# Прочитаем данные из файла
particles = read_particle_data("C:/Users/SKENDLI/Desktop/diplom/diplom_sph/sph_method/data/0.bin")

# Выбор данных для построения графика
x_column = "radius"  # Используем радиус как X
y_column = "pressure"  # Можно выбрать любую характеристику, например, давление

# Визуализируем выборочные данные
visualize_particle_data(particles, x_column, y_column)

# Визуализируем зависимость скорости от радиуса
visualize_velocity_by_radius(particles)
