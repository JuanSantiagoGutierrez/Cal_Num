class my_funcion():
    def __init__(self, funcion, derivada = None, tolerancia: float = 0.01, valor_inicial = None , valor_final = None):
        self.f   = funcion
        self.df  = derivada
        self.tol = tolerancia
        self.a   = valor_inicial
        self.b   = valor_final
    def buscar_raiz(self, metodo:str, csv = False, redondeo:int = 5, max_iter:int = 100):
        r = redondeo
        if metodo not in ("Newton-Raphson","Regula-Fansi","Biseccion"):
            return (None,'Metodos:"Newton-Raphson","Regula-Fansi","Biseccion"')
        if csv: 
            import pandas as pd
            if metodo == "Newton-Raphson":
                columns = ["iteración","a","f(a)","xr","f(xr)","error"]
                def guardar_inf(self):
                    new_line = [i, a, fa, xr,fxr,error]
                    df.loc[len(df)] = new_line
            else: 
                columns = ["iteracion","a","f(a)","b","f(b)","xr","f(xr)","error"]
                def guardar_inf(self):
                    new_line = [i, a, fa, b, fb, xr,fxr,error]
                    df.loc[len(df)] = new_line
            df = pd.DataFrame(columns=columns)
            df = df.reset_index(drop=True)
            # Función para guardar el csv
            def guardar_csv(self):
                from datetime import date
                date = date.today()
                df.to_csv(f"Resultados/{metodo}_{date}.csv")
                self.dataframe = df

        # Método Newton-Raphson    
        if metodo == "Newton-Raphson":
            # Variable intermedia: Evita que se altere el intervalo [a,b]
            ar  = round(self.a,r)
            for i in range(1, max_iter+1):
                a = round(ar,r)
                fa  = round(self.f(a),r)
                dfa = round(self.df(a),r)
                # Evitar la excepción por div0
                if dfa== 0:
                    if csv:guardar_csv(self)
                    return (None, f"Error, división por cero por que df({round(a,r)})=0")
                # Aplicar el método
                xr  = round(a - (fa)/(dfa), r)
                fxr = round(self.f(xr), r)
                # Calcular el error
                error = abs(fxr)
                # Comparar con la tolerancia
                if error <= self.tol:
                    # Si la tolerancia es mayor guardar y retornar
                    if csv: 
                        guardar_inf(self)
                        guardar_csv(self)
                    return (xr,f"Raíz encontrada en x={xr} con un error = {error}, i = {i}")
                # Si la tolerancia es menor entonces guardar y seguir iterando
                else: ar = xr
                if csv: guardar_inf(self)
            return (None, "Se ha alcanzado el máx. de iteraciones, se recomienda revisar el csv")
        # Método Biseción o Regula-Fansi
        elif metodo in ["Biseccion","Regula-Fansi"]:
            ar = self.a
            br = self.b
            i = 0
            if ar == br:
                return (None, "Intervalo inválido, por que a = b")
            while True:
                # Asignación de variables
                if i == 100: 
                    guardar_csv(self)
                    return (xr, error, b, br, f"Iteraciones máx{i}")
                a = round(ar, r)
                b = round(br, r)
                fa= round(self.f(a),r)
                fb= round(self.f(b),r)
                i = i+1
                # Intervalo valido
                if fa*fb<0:
                    # Los dos métodos difieren aquí
                    if metodo == "Biseccion":
                        xr = round(a+(b-a)/2,r)
                        fxr= round(self.f(xr),r)
                        error= abs(fxr)
                    elif metodo == "Regula-Fansi":
                        try: 
                            m = (fb-fa)/(b-a)
                            xr = round(b-fb/m,r)
                        except ZeroDivisionError:
                            return (None, f"División por cero en a={a} y b={b}")
                        fxr = round(self.f(xr),r)
                        error = abs(fxr)
                # Si es 0, asignar la raíz correspondiente
                elif fa*fb == 0:
                    if fa  == 0:
                        xr  = a
                        fxr = fa
                    else: 
                        xr  = b
                        fxr = fb
                    error = 0.0
                # Sí es posítivo el producto de las imagenes retornar no válido
                else: return (None, f"Intervalo no válido, las dos imágenes son posítivas") 
                # Guardar los datos
                if csv:
                    guardar_inf(self)
                # Comprobar el error, si es menor a la tolerancia guardar y retornar
                if error <= self.tol:
                    if csv: guardar_csv(self)
                    return (xr, f"Raíz encontrada en x = {xr} con error = {error}, i = {i}")
                # Si no, seguir iterando
                else: 
                    if fa*fxr<0:
                        br=xr
                    else: ar=xr
        else: return (None, "Error")

import numpy as np
def f(x):
    cosh = np.cosh(x * 160 / 2)
    y = cosh - 1 - 15 * x
    return y
def df(x):
    return np.cos(x)
iter = [0.0001, 0.01]
tol = 1E-8
func = my_funcion(funcion=f,derivada=df, tolerancia = tol, valor_inicial=iter[0], valor_final=iter[1])
# print(func.buscar_raiz("Newton-Raphson", redondeo=4, csv = True))
print(func.buscar_raiz("Biseccion", redondeo=12, csv=True))
print(func.buscar_raiz("Regula-Fansi", redondeo=12, csv = True))
