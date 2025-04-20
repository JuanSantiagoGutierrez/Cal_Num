import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt

def generar_indices(n: int, traslacion:int = 0) -> np.array:
    """Genera un array simétrico de n/2 enteros alrededor de cero incluyendo el 0 \n
    ,y si n es par prefiere a la derecha \n
    Traslada el array antes de devolverlo"""
    array = np.arange(-int(n/2), int(n/2)+1)
    if n%2 == 0:
        array = array[1:]
    return array + traslacion

def gen_matrix(derivada: int = 1, truncamiento: int = 2, traslacion: int = 0) -> np.ndarray:
    """Crea una matrix cuadrada para el calculo una n-derivada usando expansiones de Taylor alrededor de un punto"""
    puntos = derivada + truncamiento
    indices = generar_indices(puntos, traslacion)
    
    # Definir una matrix cuadrada vacia
    A = np.ones((puntos,puntos))

    # Los elementos de la matrix se calculan así
    def coef(i:int ,j:int) -> float:
        return math.pow(j,i) / math.factorial(i)

    # Recorrer la matrix y generar cada elemento
    for i in range(puntos):
        for j in indices:
            k = j + sum(1 for i in indices if indices[i] < 0) # Ajusta el indice j, con el indice de la matrix de forma que indices[0] -> i,0 de la matrix
            A[i][k] = round(coef(i=i,j=j), truncamiento + 2)
            
    # Retornar la matrix
    return A
        
def coef_derivada(matrix: np.ndarray[float], derivada:int =1) -> np.ndarray[float]:
    """gen_matrix() -> A\n
    Matrix A para derivadas, aplica las condiciones para la n-derivada\n
    devolviendo los coeficientes X \n 
    X -> evaluar_derivada()"""
    # Define un vector vacio y aplica las condiciones 
    v = np.zeros((len(matrix)))
    v[derivada] = 1
    # Resuelve el sistema y lo retorna
    return np.linalg.solve(a=matrix, b=v)

def evaluar_derivada(X:np.ndarray,f:np.ufunc ,x:float=0.0, dx: float = 0.1, traslacion:int = 0) -> float:
    """coef_derivada() -> X \n
    Evalua la derivada en un punto en base a los coeficientes X y a cual derivada se dirije \n
    -> f'(x)"""
    df = 0
    # Usa los coeficientes de X y los índices de los puntos para calcular df, recorriendo ambos arrays en simultaneo. 
    for i,j in zip(X, generar_indices(len(X), traslacion)): 
        df += i * f(x + j * dx )
    return df / dx

def derivar(x:float,f:np.ufunc,dx:float,derivada: int =1, orden:int = 2, traslacion:int=0)-> float:
    """Deriva una funcion númericamente, usando expansiones de Taylor y orden de Landau"""
    # Calcular la matrix, por defecto da la 1er derivada
    A = gen_matrix(derivada=derivada,truncamiento=orden,traslacion= t)
   
    # Calcular los coeficientes
    X = coef_derivada(matrix=A,derivada=derivada)
    
    # Evaluar la derivada
    df = evaluar_derivada(X=X, f=f, x=x, dx=dx,traslacion=traslacion)
    return df

if __name__ == "__main__":
    # Variables
    t = 1
    derivada = 1
    orden = 2
    dx = 2
    DF = pd.read_csv("TPs/Reporte diario COVID-19.csv")
    def contagios(dia):
        """contagios(dia)""" # Aclarar la funcion para la documentación
        return DF["acumulado"][dia-1]
    L = []
    for dia in DF["dia"]:
        t = 0
        while True:
            try:
                df = derivar(x=dia, f=contagios, dx=dx, derivada=derivada, orden=orden, traslacion=t)
                break
            except:
                if dia < len(DF)//2: t+=1 
                else: t-=1
        L += [df]
    DF = DF.assign(df1_2_Orden_2_p = L)
    DF.to_csv("Reporte_1.csv")
    
    
    import matplotlib.pyplot as plt
    plt.figure(figsize=(8,6))
    plt.plot(DF["dia"], DF["acumulado"])
    plt.savefig("graf_acumulado.png")
    plt.figure(figsize=(8,6))
    plt.plot(DF["dia"], DF["nuevo"])
    plt.savefig("graf_nuevos.png")
    plt.figure(figsize=(8,6))
    plt.plot(DF["dia"], L)
    plt.savefig("graf_L.png")