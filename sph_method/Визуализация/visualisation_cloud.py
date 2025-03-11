import os
import struct
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from multiprocessing import Pool, cpu_count

# Функция для загрузки данных о частицах из бинарного файла
def load_particles_from_binary(file_path):
    particles = []
    with open(file_path, 'rb') as f:
        Dt = struct.unpack('d', f.read(8))[0]
        times = struct.unpack('d', f.read(8))[0]
        np_count = struct.unpack('i', f.read(4))[0]
        
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
            neighbors = struct.unpack('i', f.read(4))[0]
            
            particles.append({
                'x': x, 'y': y, 'density': density, 'mass': mass,
                'pressure': pressure, 'velocityX': velocityX,
                'velocityY': velocityY, 'energy': energy,
                'smoothingLength': smoothingLength, 'neighbors': neighbors 
            })
    return Dt, times, particles

# Функция обработки одного файла
def process_single_file(args):
    file_path, output_folder, size_scat = args
    Dt, times, particles = load_particles_from_binary(file_path)
    
    fig, axes = plt.subplots(1, 2, figsize=(19.2, 10.8), dpi=96)
    
    x_coords = [p['x'] for p in particles]
    y_coords = [p['y'] for p in particles]
    
    pressure_values = [p['pressure'] for p in particles]
    sc1 = axes[0].scatter(x_coords, y_coords, c=pressure_values, cmap='magma', s=size_scat, edgecolors="none")
    axes[0].set_title("Давление p", fontsize=20, pad=20)
    axes[0].tick_params(axis='both', labelsize=24)
    axes[0].set_xlabel("x", fontsize=24, labelpad=10)
    axes[0].set_ylabel("y", fontsize=24, labelpad=10)
    axes[0].set_aspect("equal")
    
    divider1 = make_axes_locatable(axes[0])
    cax1 = divider1.append_axes("right", size="5%", pad=0.1)
    cbar1 = plt.colorbar(sc1, cax=cax1)
    cbar1.ax.tick_params(labelsize=16)
    
    density_values = [p['density'] for p in particles]
    sc2 = axes[1].scatter(x_coords, y_coords, c=density_values, cmap='viridis', s=size_scat, edgecolors="none")
    axes[1].set_title(r'Плотность $\rho$', fontsize=20, pad=20)
    axes[1].tick_params(axis='both', labelsize=24)
    axes[1].set_xlabel("x", fontsize=24, labelpad=10)
    axes[1].set_ylabel("y", fontsize=24, labelpad=10)
    axes[1].set_aspect("equal")
    
    divider2 = make_axes_locatable(axes[1])
    cax2 = divider2.append_axes("right", size="5%", pad=0.1)
    cbar2 = plt.colorbar(sc2, cax=cax2)
    cbar2.ax.tick_params(labelsize=16)
    
    sample_particles = np.random.choice(particles, size=max(1, len(particles) // 10000), replace=False)
    for p in sample_particles:
        velocity_magnitude = np.hypot(p['velocityX'], p['velocityY'])
        if velocity_magnitude > 0:
            velocityX_normalized = p['velocityX'] / velocity_magnitude
            velocityY_normalized = p['velocityY'] / velocity_magnitude
        else:
            velocityX_normalized = 0
            velocityY_normalized = 0
        
        axes[1].quiver(p['x'], p['y'], 
                      velocityX_normalized * 0.1, 
                      velocityY_normalized * 0.1, 
                      angles='xy', scale_units='xy', scale=1, 
                      color='r', width=0.003)
    
    output_file = os.path.join(output_folder, f"{os.path.basename(file_path).split('.')[0]}.jpeg")
    plt.tight_layout()
    plt.savefig(output_file, dpi=96, bbox_inches='tight', pad_inches=0.1)
    plt.close()

# Основная функция с распараллеливанием
def process_files_parallel(input_folder, output_folder, size_scat):
    os.makedirs(output_folder, exist_ok=True)
    
    bin_files = [f for f in os.listdir(input_folder) if f.endswith('.bin')]
    bin_files.sort(key=lambda x: int(x.split('.')[0]))
    
    with Pool(processes=cpu_count()) as pool:
        args_list = [
            (os.path.join(input_folder, f), output_folder, size_scat) 
            for f in bin_files
        ]
        pool.map(process_single_file, args_list)

if __name__ == '__main__':
    input_folder = r'C:/Users/SKENDLI/Desktop/diplom/diplom_sph/sph_method/data'
    output_folder = r'C:/Users/SKENDLI/Desktop/diplom/diplom_sph/sph_method/Визуализация/Images'
    size_scat = 3
    
    process_files_parallel(input_folder, output_folder, size_scat)
