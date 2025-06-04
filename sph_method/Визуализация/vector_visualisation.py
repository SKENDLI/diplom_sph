import os
import struct
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# Настройки путей
input_dir = "C:/Users/SKENDLI/Desktop/diplom/diplom_sph/sph_method/data"
output_dir = "C:/Users/SKENDLI/Desktop/diplom/diplom_sph/sph_method/Визуализация/Images/Vectors1"
os.makedirs(output_dir, exist_ok=True)

# Функция загрузки данных
def load_particles_from_binary(file_path):
    particles = []
    with open(file_path, 'rb') as f:
        Dt = struct.unpack('d', f.read(8))[0]
        time = struct.unpack('d', f.read(8))[0]
        np_count = struct.unpack('i', f.read(4))[0]
        
        for _ in range(np_count):
            data = struct.unpack('dddddddddi', f.read(76))
            particles.append({
                'x': data[0], 'y': data[1],
                'vx': data[5], 'vy': data[6]
            })
    return time, pd.DataFrame(particles)

# Функция обработки одного файла
def process_single_file(file_path, output_dir):
    try:
        time, df = load_particles_from_binary(file_path)
        
        plt.figure(figsize=(12, 12))
        
        # 1. Отрисовка всех частиц (точки)
        plt.scatter(df['x'], df['y'], s=0.1, alpha=0.3, color='black')
        
        # 2. Создание сетки для интерполяции (уменьшенный размер)
        grid_x, grid_y = np.meshgrid(
            np.linspace(df['x'].min(), df['x'].max(), 50),  # Уменьшено до 50
            np.linspace(df['y'].min(), df['y'].max(), 50)   # Уменьшено до 50
        )
        
        # 3. Интерполяция скоростей на сетку
        grid_vx = griddata((df['x'], df['y']), df['vx'], (grid_x, grid_y), method='linear')
        grid_vy = griddata((df['x'], df['y']), df['vy'], (grid_x, grid_y), method='linear')
        
        # 4. Отрисовка потоковых линий
        plt.streamplot(grid_x, grid_y, grid_vx, grid_vy, density=1.5, color='blue', linewidth=1)
        
        plt.title(f"Time = {time:.3f} (Particles: {len(df)})")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.grid(alpha=0.2)
        
        # Автоматическое масштабирование с небольшим отступом
        x_margin = 0.1 * (df['x'].max() - df['x'].min())
        y_margin = 0.1 * (df['y'].max() - df['y'].min())
        plt.xlim(df['x'].min()-x_margin, df['x'].max()+x_margin)
        plt.ylim(df['y'].min()-y_margin, df['y'].max()+y_margin)
        
        output_path = os.path.join(output_dir, f"plot_{time:.3f}.png")
        plt.savefig(output_path, dpi=200, bbox_inches='tight')
        plt.close()
    except Exception as e:
        print(f"Error in {file_path}: {str(e)}")

# Основная функция с последовательной обработкой
def process_files_sequentially(input_dir, output_dir):
    files = [f for f in os.listdir(input_dir) if f.endswith(".bin")]
    files.sort(key=lambda x: float(x.split('_')[-1].replace('.bin', '')))
    
    # Обрабатываем все файлы
    for f in files:
        file_path = os.path.join(input_dir, f)
        process_single_file(file_path, output_dir)

if __name__ == "__main__":
    process_files_sequentially(input_dir, output_dir)
    print(f"Visualization complete. Results saved in {output_dir}")
