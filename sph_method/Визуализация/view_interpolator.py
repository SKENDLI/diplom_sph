import os
import struct
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from scipy.interpolate import griddata

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
            neighbors = struct.unpack('i', f.read(4))[0]  # Число соседей
            
            particles.append({
                'x': x, 'y': y, 'density': density, 'mass': mass,
                'pressure': pressure, 'velocityX': velocityX,
                'velocityY': velocityY, 'energy': energy,
                'smoothingLength': smoothingLength, 'neighbors': neighbors
            })
    return Dt, times, particles

def process_files(input_folder, output_folder):
    bin_files = [f for f in os.listdir(input_folder) if f.endswith('.bin')]
    bin_files.sort(key=lambda x: int(x.split('.')[0]))

    for bin_file in bin_files:
        file_path = os.path.join(input_folder, bin_file)
        
        # Загрузка данных
        Dt, times, particles = load_particles_from_binary(file_path)

        # Подготовка данных
        x = np.array([p['x'] for p in particles])
        y = np.array([p['y'] for p in particles])
        pressure = np.array([p['pressure'] for p in particles])
        density = np.array([p['density'] for p in particles])

        # Создание сетки для интерполяции
        grid_resolution = 1000
        xi = np.linspace(x.min(), x.max(), grid_resolution)
        yi = np.linspace(y.min(), y.max(), grid_resolution)
        xi, yi = np.meshgrid(xi, yi)

        # Интерполяция давления
        pressure_grid = griddata(
            (x, y), pressure, 
            (xi, yi), method='linear', fill_value=0
        )

        # Интерполяция плотности
        density_grid = griddata(
            (x, y), density,
            (xi, yi), method='linear', fill_value=0
        )

        # Визуализация
        fig, axes = plt.subplots(1, 2, figsize=(19.2, 10.8), dpi=96)
        extent = [x.min(), x.max(), y.min(), y.max()]

        # Давление
        im1 = axes[0].imshow(
            pressure_grid.T, extent=extent, origin='lower',
            cmap='magma', aspect='auto', interpolation='bicubic'
        )
        axes[0].set_title("Давление", fontsize=20, pad=15)
        axes[0].tick_params(labelsize=14)

        # Плотность
        im2 = axes[1].imshow(
            density_grid.T, extent=extent, origin='lower',
            cmap='viridis', aspect='auto', interpolation='bicubic'
        )
        axes[1].set_title("Плотность", fontsize=20, pad=15)
        axes[1].tick_params(labelsize=14)

        # Цветовые бары
        for ax, im in zip(axes, [im1, im2]):
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            plt.colorbar(im, cax=cax).ax.tick_params(labelsize=12)

        # Вектора скорости
        sample = np.random.choice(particles, size=min(500, len(particles)), replace=False)
        vx = np.array([p['velocityX'] for p in sample])
        vy = np.array([p['velocityY'] for p in sample])
        xs = np.array([p['x'] for p in sample])
        ys = np.array([p['y'] for p in sample])
        
        speed = np.sqrt(vx**2 + vy**2)
        mask = speed > 0.1 * speed.max()  # Фильтр малых скоростей
        
        axes[1].quiver(
            xs[mask], ys[mask], 
            vx[mask]/speed[mask], vy[mask]/speed[mask],
            scale_units='xy', scale=10, 
            width=0.003, color='white'
        )

        # Сохранение
        plt.tight_layout(pad=2.0)
        output_file = os.path.join(output_folder, f"{bin_file.split('.')[0]}.png")
        plt.savefig(output_file, bbox_inches='tight', dpi=150)
        plt.close()

# Конфигурация
input_folder = 'C:/Users/SKENDLI/Desktop/diplom/diplom_sph/sph_method/data'
output_folder = 'C:/Users/SKENDLI/Desktop/diplom/diplom_sph/sph_method/Визуализация/Images'
os.makedirs(output_folder, exist_ok=True)

# Запуск обработки
process_files(input_folder, output_folder)
