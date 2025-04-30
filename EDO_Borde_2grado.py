import numpy as np
import math
class EDO_Borde_2grado():
    """Coeficientes de una EDO de segundo grado\n
    a*d²y/dx²+b*dy/dx+c*(y-ext)=0"""
    def __init__(self, particiones:int,intervalo:tuple,borde:tuple, ext:float=0, coef_a:float=0,coef_b:float = 0,coef_c:float=0):
        self.a = coef_a
        self.b = coef_b
        self.c = coef_c
        self.ext = ext
        self.borde = borde
        self.particiones = particiones
        self.inicio = intervalo[0]
        self.final = intervalo[1]
        self.dx = (intervalo[1]-intervalo[0])/particiones
        self.ecu = f"{self.a}*d²y/dx + {self.b}*dy/dx + {self.c}*(y-{self.ext}) = 0"
        
    def generar_matriz_vector(self):
        self.espacio = np.linspace(self.inicio, self.final, self.particiones)
        self.A = np.zeros((len(self.espacio),len(self.espacio)))
        self.A[0,0] = 1
        self.A[-1,-1]=1
        for i in range(1,len(self.espacio)-1):
            self.A[i][i] = self.dx**2 * self.c - 2*self.a
            self.A[i][i-1] = self.a-self.dx*self.b
            self.A[i][i+1] = self.a+self.dx*self.b
        self.V = np.zeros(len(self.espacio))
        self.V[0] = self.borde[0]
        self.V[-1] = self.borde[-1]
        for v in range(1,len(self.espacio)-1):
            self.V[v] = self.c*self.dx ** 2 * self.ext
    def resolver_EDO(self):
        """Resuelve la EDO, con la matriz A y el vector V, retornando un dicionario con las respuestas"""
        X = np.linalg.solve(self.A, self.V)
        self.solucion = {"x":tuple(self.espacio), "y": tuple(X)}
    def presentar_solucion(self):
        import matplotlib.pyplot as plt
        # Graficar los puntos de la solución
        plt.plot(EDO.solucion["x"], EDO.solucion["y"], marker='o', label=f'Solución EDO')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title(f'Solución de la EDO de segundo grado:{self.ecu}')
        plt.legend()
        plt.grid(True)
        plt.savefig("Resultados/Resultado_EDO_Contorno_2grado")
        plt.show()

            
if __name__ == "__main__":
    
    def Vibrar(k,modo = 0, intervalo = (0,1),particiones = 100, borde = (0,0)):
        espacio = np.linspace(intervalo[0], intervalo[1], particiones)
        A = np.zeros((len(espacio)-2,len(espacio)-2))
        for i in range(0,len(A)):
            A[i][i] = -2 * k
            try:
                A[i][i-1] = 1
                A[i][i+1] = 1
            except IndexError:
                pass 
        A[0][-1] = 0
        autovalores, autovectores = np.linalg.eig(A)
        print("autovalor", autovalores[modo])
        Y = [borde[0]]
        v = autovectores
        for i in v:
            Y.append(i[modo])
        Y.append(borde[1])
        return espacio, Y
    
    import matplotlib.pyplot as plt
    
    
    plt.figure(figsize=(8,6))
    
    for i in range(3,6):
        x,y = Vibrar(6, i)
        plt.plot(x,y)
    
    plt.savefig("Resultados/EDO---")
    
        
        