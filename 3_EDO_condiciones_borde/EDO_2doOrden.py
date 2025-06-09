import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import math

m = 1
c = 1
k = 0.1

y0 = 1
dy0 = 0

Z = np.array((y0, dy0))
t = np.array([0.0, 100.0])

def dZdt(t, Z):
    return [Z[1], (-c*Z[1]-k*Z[0])/m]

sol = solve_ivp(dZdt, t, Z, method='RK45', dense_output=True)

t_eval = np.linspace(t[0], t[1], 1000)
y_sol = sol.sol(t_eval)

print(y_sol)
pos = y_sol[0].T
vel = y_sol[1].T

plt.plot(pos, vel, label="estado", color="blue")
plt.xlabel('posicion')
plt.ylabel('velocidad')
plt.title(' usando RK45')
plt.grid()
plt.savefig("Resultados/EDO_2doOrden")

plt.figure(figsize=(8,6))
plt.plot(t_eval, y_sol[0].T, label="posicion", color="red")
plt.xlabel('tiempo')
plt.ylabel('posicion')
plt.title(' usando RK45')
plt.grid()
plt.savefig("Resultados/EDO_2doOrden_1")
t = 100
y_1= 1
Y = np.array([y_1, y0])
dt = 0.1
x1 = lambda x0,x_1,dt : (-k+2*m/(math.pow(dt,2)))/(m/(math.pow(dt,2))+(c)/(2*dt)) *x0 + ((-m)/(math.pow(dt,2))+(c)/(2*dt))/((m/(math.pow(dt,2))+(c)/(2*dt)))*x_1 
for i in range(t-1):
    X = x1(x0=Y[-1], x_1=Y[-2], dt=dt)
    Y = np.append(Y, X)
Y = Y[1:]

print(Y)
plt.figure(figsize=(8,6))
plt.plot(range(t), Y, label="posicion", color="red")
plt.xlabel('tiempo')
plt.ylabel('posicion')
plt.title('Diferencia Central')
plt.grid()
plt.savefig("Resultados/EDO_2doOrden_2")