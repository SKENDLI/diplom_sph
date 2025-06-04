import struct
import numpy as np
import matplotlib.pyplot as plt


def read_particle_data(file_path):
    with open(file_path, "rb") as file:
        Dt = struct.unpack('d', file.read(8))[0]
        times = struct.unpack('d', file.read(8))[0]
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
            particle["r"] = np.sqrt(particle["x"]**2 + particle["y"]**2)
            particle["v"] = np.sqrt(particle["velocityX"]**2 + particle["velocityY"]**2)
            particles.append(particle)
        return particles


def get_axis_limits(data, padding=0.02):
    """Вычисляет границы оси с небольшим отступом"""
    min_val = min(data)
    max_val = max(data)
    span = max_val - min_val
    if span == 0:  # если все значения одинаковые
        span = 0.1 * abs(min_val) if min_val != 0 else 0.1
    padding = span * padding
    return min_val - padding, max_val + padding


def visualize_particle_data(particles, x_key, y_key):
    x_data = np.array([particle[x_key] for particle in particles])
    y_data = np.array([particle[y_key] for particle in particles])
    
    if x_key == "r":
        x_data = np.concatenate([-x_data, x_data])
        y_data = np.concatenate([y_data, y_data])
    
    # Фильтрация NaN и бесконечных значений
    mask = np.isfinite(x_data) & np.isfinite(y_data)
    x_data = x_data[mask]
    y_data = y_data[mask]
    
    if len(x_data) == 0 or len(y_data) == 0:
        print("Нет данных для отображения после фильтрации!")
        return
    
    plt.figure(figsize=(10, 8))
    
    # Определяем подписи осей
    x_label = "r" if x_key == "r" else x_key
    y_label = {
        "density": "ρ",
        "v": "V",
        "mass": "m",
        "energy": "E"
    }.get(y_key, y_key)
    
    plt.scatter(x_data, y_data, s=2, color='blue')
    plt.xlabel(x_label, fontsize=24)
    plt.ylabel(y_label, fontsize=24)
    plt.title(f'{ "Плотность ρ" if y_key == "density" else "Скорость V" } ({x_label})', fontsize=24)
    
    # Устанавливаем границы осей только по фактическим данным
    x_min, x_max = get_axis_limits(x_data)
    y_min, y_max = get_axis_limits(y_data)
    plt.xlim(0, 3.0)
    plt.ylim(y_min, y_max)
    
    # Увеличение размера числовых меток на осях
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.grid(True)
    plt.tight_layout(pad=2.0)  # Увеличиваем отступ
    plt.show()


def visualize_velocity_by_radius(particles):
    radius_data = np.array([particle["r"] for particle in particles])
    velocity_data = np.array([particle["v"] for particle in particles])

    # Фильтрация данных
    mask = (radius_data >= 0) & (velocity_data < 10) & np.isfinite(radius_data) & np.isfinite(velocity_data)
    radius_data = radius_data[mask]
    velocity_data = velocity_data[mask]
    
    if len(radius_data) == 0:
        print("Нет данных для визуализации после фильтрации!")
        return
    
    # Усреднение скорости по радиусам
    radius_velocity_dict = {}
    for r, v in zip(radius_data, velocity_data):
        r_rounded = round(r, 2)
        if r_rounded not in radius_velocity_dict:
            radius_velocity_dict[r_rounded] = []
        radius_velocity_dict[r_rounded].append(v)

    avg_radius = []
    avg_velocity = []
    for r, velocities in sorted(radius_velocity_dict.items()):
        avg_radius.append(r)
        avg_velocity.append(np.mean(velocities))
    
    plt.figure(figsize=(10, 8))
    plt.plot(avg_radius, avg_velocity, color='blue', linestyle='--', linewidth=2)
    plt.xlabel('r', fontsize=24)
    plt.ylabel('V', fontsize=24)
    plt.title('Скорость V(r)', fontsize=24)
    
    # Устанавливаем границы осей
    x_min, x_max = get_axis_limits(avg_radius)
    y_min, y_max = get_axis_limits(avg_velocity)
    plt.xlim(0, 3.0)
    plt.ylim(0, y_max)
    
    # Увеличение размера числовых меток на осях
    plt.xticks(fontsize=20)  # Увеличено с 16 до 20
    plt.yticks(fontsize=20)  # Увеличено с 16 до 20
    plt.grid(True)
    plt.tight_layout(pad=2.0)
    plt.show()


# Пример использования
particles = read_particle_data("C:/Users/SKENDLI/Desktop/diplom/diplom_sph/sph_method/data/0.bin")
visualize_particle_data(particles, "r", "density")
visualize_velocity_by_radius(particles)
