
import numpy as np
import pandas as pd

S = 10
R = 28
B = 8/3

x_next = lambda x,y,z,dt:x + dt * (S * (y-x))
y_next = lambda x,y,z,dt:y + dt * (-x*z + R*x-y)
z_next = lambda x,y,z,dt:z + dt * (x*y - B*z)

val_i = np.array([1, 0, 0])
dt = 0.01
t_final = 50

DF = pd.DataFrame({"x": [], "y": [], "z":[], "t":[]})
def save_csv(DF2,m_xyz: np.array, t:float):
    nueva_fila = [m_xyz[0],m_xyz[1],m_xyz[2],t]
    DF2.loc[len(DF2)] = nueva_fila
    return DF2
def euler(i_xyz: np.array, t_final:float, dt:float, t_inicial = 0):
    t = t_inicial
    save_csv(DF,i_xyz, t)
    while t <= t_final:
        xyz = np.zeros(len(i_xyz))
        xyz[0] = x_next(i_xyz[0], i_xyz[1], i_xyz[2], dt) 
        xyz[1] = y_next(i_xyz[0], i_xyz[1], i_xyz[2], dt) 
        xyz[2] = z_next(i_xyz[0], i_xyz[1], i_xyz[2], dt) 
        t += dt
        i_xyz = xyz
        save_csv(DF,xyz, t)

euler(val_i, t_final=t_final, dt=dt)
DF.to_csv("Resultados/Atractor de Lorentz.csv")

DF = pd.read_csv("Resultados/Atractor de Lorentz.csv")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


# Preparar datos
df_sorted = DF.sort_values('t')  # Ordenar por tiempo

# Inicializar figura
fig = plt.figure(figsize=(10, 7), facecolor="w")
ax = fig.add_subplot(111, projection='3d')

# Función de actualización para cada frame
def update(frame):
    ax.clear()
    # Graficar hasta el frame actual
    data = df_sorted[df_sorted['t'] <= frame]
    sc = ax.scatter(data['x'], data['y'], data['z'], c=data['t'], cmap='viridis')
    ax.set_title(f'Tiempo t = {frame:.2f}')
    return sc

# Crear animación
ani = FuncAnimation(
    fig, update,
    frames=np.linspace(DF['t'].min(), DF['t'].max(), 150),  # Pasos de tiempo
    interval=100,  # Tiempo entre frames (ms)
    blit=False
)

ax.set_axis_off()
ax.grid(False)
ani.save('Resultados/animacion_Atractor.mp4', writer='ffmpeg')
print("Animacion lista")