import os
import struct
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count

# Настройки путей
input_dir = "C:/Users/SKENDLI/Desktop/diplom/diplom_sph/sph_method/data"
output_dir = "C:/Users/SKENDLI/Desktop/diplom/diplom_sph/sph_method/Визуализация/Images/Vectors"
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
def process_single_file(args):
    file_path, output_dir = args
    try:
        time, df = load_particles_from_binary(file_path)
        
        plt.figure(figsize=(12, 12))
        
        # 1. Отрисовка всех частиц (точки)
        plt.scatter(df['x'], df['y'], s=0.1, alpha=0.3, color='black')
        
        # 2. Отрисовка векторов (подвыборка)
        sample_size = min(500, len(df))  # Фиксированное количество векторов
        sample = df.sample(n=sample_size)
        
        # Масштабирование векторов для лучшей визуализации
        velocity_scale = 0.001 * (df['x'].max() - df['x'].min())
        plt.quiver(sample['x'], sample['y'], 
                   sample['vx'], sample['vy'], 
                   scale=1/velocity_scale, width=0.002, 
                   color='blue', alpha=0.7)
        
        plt.title(f"Time = {time:.3f} (Particles: {len(df)}, Vectors: {sample_size})")
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

# Основная функция с распараллеливанием
def process_files_parallel(input_dir, output_dir):
    files = [f for f in os.listdir(input_dir) if f.endswith(".bin")]
    files.sort(key=lambda x: float(x.split('_')[-1].replace('.bin', '')))
    
    with Pool(processes=cpu_count()) as pool:
        args_list = [
            (os.path.join(input_dir, f), output_dir) 
            for f in files
        ]
        pool.map(process_single_file, args_list)

if __name__ == "__main__":
    process_files_parallel(input_dir, output_dir)
    print(f"Visualization complete. Results saved in {output_dir}")
