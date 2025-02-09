import os
import struct
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

# Функция для загрузки данных о частицах из бинарного файла
def load_particles_from_binary(file_path):
    particles = []
    with open(file_path, 'rb') as f:
        # Считываем Dt
        Dt = struct.unpack('d', f.read(8))[0]
        # Считываем times
        times = struct.unpack('d', f.read(8))[0]
        # Считываем количество частиц
        np_count = struct.unpack('i', f.read(4))[0]

        # Считываем данные каждой частицы
        for _ in range(np_count):
            x = struct.unpack('d', f.read(8))[0]
            y = struct.unpack('d', f.read(8))[0]
            density = struct.unpack('d', f.read(8))[0]
            mass = struct.unpack('d', f.read(8))[0]
            pressure = struct.unpack('d', f.read(8))[0]
            velocityX = struct.unpack('d', f.read(8))[0]
            velocityY = struct.unpack('d', f.read(8))[0]
            energy = struct.unpack('d', f.read(8))[0]
            smoothingLength = struct.unpack('d', f.read(8))[0]
            neighbors  = struct.unpack('i', f.read(4))[0]  # Число соседей
            
            particles.append({
                'x': x, 'y': y, 'density': density, 'mass': mass,
                'pressure': pressure, 'velocityX': velocityX,
                'velocityY': velocityY, 'energy': energy,
                'smoothingLength': smoothingLength, 'neighbors': neighbors 
            })
    return Dt, times, particles

# Функция для пакетной обработки файлов
def process_files(input_folder, output_folder):
    # Получаем список всех файлов в папке, отфильровав только .bin файлы
    bin_files = [f for f in os.listdir(input_folder) if f.endswith('.bin')]
    bin_files.sort(key=lambda x: int(x.split('.')[0]))  # Сортируем по числовым именам файлов

    # Обрабатываем каждый файл
    for bin_file in bin_files:
        file_path = os.path.join(input_folder, bin_file)
        
        # Загрузка данных о частицах
        Dt, times, particles = load_particles_from_binary(file_path)

        # Визуализация
        fig, axes = plt.subplots(1, 2, figsize=(19.2, 10.8), dpi=96)

        # Получаем соответствующие данные для выбранной характеристики
        x_coords = [p['x'] for p in particles]
        y_coords = [p['y'] for p in particles]
        
        # Давление (с палитрой magma)
        pressure_values = [p['pressure'] for p in particles]
        sc1 = axes[0].scatter(x_coords, y_coords, c=pressure_values, cmap='magma', s=size_scat, edgecolors="none")
        axes[0].set_title("Давление p", fontsize=20, pad=20)
        axes[0].set_xlabel("x", fontsize=24, labelpad=10)
        axes[0].set_ylabel("y", fontsize=24, labelpad=10)
        axes[0].set_aspect("equal", adjustable='box')  # Это гарантирует одинаковое соотношение сторон

        # Увеличиваем шрифт значений на осях
        axes[0].tick_params(axis='x', labelsize=24)
        axes[0].tick_params(axis='y', labelsize=24)

        # Add colorbar for pressure
        divider1 = make_axes_locatable(axes[0])
        cax1 = divider1.append_axes("right", size="5%", pad=0.1)
        cbar1 = plt.colorbar(sc1, cax=cax1)
        cbar1.ax.tick_params(labelsize=16)

        # Плотность (с палитрой viridis)
        density_values = [p['density'] for p in particles]
        sc2 = axes[1].scatter(x_coords, y_coords, c=density_values, cmap='viridis', s=size_scat, edgecolors="none")
        axes[1].set_title(r'Плотность $\rho$', fontsize=20, pad=20)
        axes[1].set_xlabel("x", fontsize=24, labelpad=10)
        axes[1].set_ylabel("y", fontsize=24, labelpad=10)
        axes[1].set_aspect("equal", adjustable='box')  # Это гарантирует одинаковое соотношение сторон

        # Увеличиваем шрифт значений на осях
        axes[1].tick_params(axis='x', labelsize=24)
        axes[1].tick_params(axis='y', labelsize=24)

        # Add colorbar for density
        divider2 = make_axes_locatable(axes[1])
        cax2 = divider2.append_axes("right", size="5%", pad=0.1)
        cbar2 = plt.colorbar(sc2, cax=cax2)
        cbar2.ax.tick_params(labelsize=16)

        # Добавляем вектора скорости для частицы
        # Выбираем 1% частиц для отображения векторов скорости
        sample_particles = np.random.choice(particles, size=max(1, len(particles) // 10000), replace=False)

        for p in sample_particles:
            # Нормализуем вектор скорости
            velocity_magnitude = np.sqrt(p['velocityX']**2 + p['velocityY']**2)
            if velocity_magnitude > 0:
                velocityX_normalized = p['velocityX'] / velocity_magnitude
                velocityY_normalized = p['velocityY'] / velocity_magnitude
            else:
                velocityX_normalized = 0
                velocityY_normalized = 0
            
            # Масштабируем вектор скорости (например, умножаем на 0.1 для визуализации)
            scale_factor = 0.1
            axes[1].quiver(p['x'], p['y'], velocityX_normalized * scale_factor, velocityY_normalized * scale_factor, 
                           angles='xy', scale_units='xy', scale=1, color='r', width=0.003)

        # Сохранение изображения в разрешении 1920x1080
        plt.tight_layout()
        output_file = os.path.join(output_folder, f'{bin_file.split(".")[0]}.jpeg')
        plt.savefig(output_file, dpi=96, bbox_inches='tight', pad_inches=0.1, transparent=True)
        plt.close()  # Закрыть фигуру, чтобы не занимать память

# Папка с бинарными файлами
input_folder = 'C:/Users/SKENDLI/Desktop/diplom/diplom_sph/sph_method/data'  # Путь к папке с .bin файлами
size_scat = 3
# Папка для сохранения картинок
output_folder = 'C:/Users/SKENDLI/Desktop/diplom/diplom_sph/sph_method/Визуализация/Images'  # Путь к папке для сохранения картинок

# Создаем папку для изображений, если она не существует
os.makedirs(output_folder, exist_ok=True)

# Запуск обработки файлов с разными палитрами
process_files(input_folder, output_folder)
