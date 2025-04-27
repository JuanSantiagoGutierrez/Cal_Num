import numpy as np
import math
import pandas as pd
import scipy as sc

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
    def coef(i: np.ndarray, j: np.ndarray) -> np.ndarray:
        return (j ** i) / np.vectorize(math.factorial)(i)

    # Crear mallas de índices (i, j) para la matriz
    i_grid, j_grid = np.meshgrid(np.arange(puntos), indices, indexing='ij')
    
    # Generar la matriz de una vez
    A = coef(i_grid, j_grid)
    
    # Redondear (opcional, solo si es necesario para precisión visual)
    A = np.round(A, truncamiento + 2)
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

def evaluar_derivada(X:np.ndarray,f:np.ufunc ,x:float=0.0, dx: float = 0.1, traslacion:int = 0, derivada:int = 1) -> float:
    """coef_derivada() -> X \n
    Evalua la derivada en un punto en base a los coeficientes X y a cual derivada se dirije \n
    -> f'(x)"""
    df = 0
    # Usa los coeficientes de X y los índices de los puntos para calcular df, recorriendo ambos arrays en simultaneo. 
    for i,j in zip(X, generar_indices(len(X), traslacion)): 
        df += i * f(x + j * dx )
    return df / math.pow(dx,derivada)

def derivar(x:float,f:np.ufunc,dx:float,derivada: int =1, orden:int = 2, traslacion:int=0)-> float:
    """Deriva una funcion númericamente, usando expansiones de Taylor y orden de Landau"""
    # Calcular la matrix, por defecto da la 1er derivada
    A = gen_matrix(derivada=derivada,truncamiento=orden,traslacion= traslacion)
   
    # Calcular los coeficientes
    X = coef_derivada(matrix=A,derivada=derivada)
    
    # Evaluar la derivada
    df = evaluar_derivada(X=X, f=f, x=x, dx=dx,traslacion=traslacion, derivada=derivada)
    return df

def extrapolacion_Richarson(D0:float, D1:float, dx0:float, dx1:float, n:int)-> tuple[float]:
    """Aproximación de Richarson, El orden aumenta en e^(n+2)"""
    df = np.round(D0 + ((D0 - D1)/(math.pow((dx1/dx0),n)-1)),5) 
    err = abs(np.round((D0-D1)/(math.pow(dx1,n)-math.pow(dx0,n)),5))
    return (df, err)

"""1. Aproxime la derivada primera para cada día con un esquema de orden 2 y un paso de 1 día."""
def ejer():
    # Variables
    traslacion = 5
    derivada = 2
    orden = 3
    pasos = 1 # 
    def contagios(dia):
        return DF["acumulado"][dia-1]
    L = []
    for dia in DF["dia"]:
        traslacion = 0
        while True:
            try:
                df = derivar(x=dia, f=contagios, dx=pasos, derivada=derivada, orden=orden, traslacion=traslacion)
                break
            except:
                if dia < len(DF)//2: traslacion+=1 
                else: traslacion-=1
        L += [df]
    DF = DF.assign(df_ejer1 = L)
    return 
    
if __name__ == "__main__":
    derivada = 1
    A = gen_matrix(derivada, 2, 0)
    X = coef_derivada(A, derivada)
    print("Matriz :\n", A)
    print("Coef.:\n", X)
    
    