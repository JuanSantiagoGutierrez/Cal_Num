import pandas as pd
class EDO():
    def __init__(self, dfdx, x0, f0, dx):
        self.dfdx = dfdx
        self.x0 = x0
        self.f0 = f0
        self.dx = dx
        
    def derivada0(self, x0, f0):
        self.derivada = self.dfdx(x0= x0, dx = self.dx,f0= f0)
        return self.derivada
    
    def EDO_Euler(self, x0, f0):
        self.f1 = f0 + self.dx * self.derivada0(x0, f0)
        return self.f1
    
    def predicion_EDO_Euler(self,pasos, csv_name = "Resultado EDO_Euler.csv"):    
        DF = pd.DataFrame({"x":[],"f(x)":[],"dfdx[x]":[], "f(x+dx)":[]})
        f0 = self.f0
        x0 = self.x0
        dfdx = self.derivada0(x0 = self.x0, f0 = self.f0)
        f1 = self.EDO_Euler(x0=x0, f0=f0)
        for i in range(1, pasos+1):
            L = [x0, f0, dfdx, f1]
            DF.loc[len(DF)] = L
            x0 = self.x0 + i * self.dx
            f0 = f1
            dfdx = self.derivada0(x0=x0, f0=f0)
            f1 = self.EDO_Euler(x0=x0, f0=f0)
        DF.to_csv(f"Resultados/{csv_name}")
        return DF
    def predicion_EDO_Euler_Optimizado(self, pasos, csv_name = "EDO_Euler_optimizado.csv"):
        DF = pd.DataFrame({"x":[],"f(x)":[],"f_ev":[],"dfdx_ev":[],"dfdx[x]":[], "f(x+dx)":[]})
        f0 = self.f0
        x0 = self.x0
        dfdx_ev = self.derivada0(x0 = self.x0, f0 = self.f0)
        f_ev = self.EDO_Euler(x0=x0, f0=f0)
        dfdx = self.derivada0(x0 = self.x0+self.dx, f0 = f_ev)
        f1 = f0 + self.dx*(dfdx+dfdx_ev) / 2
        for i in range(1, pasos+1):
            L = [x0, f0,dfdx_ev, f_ev, dfdx, f1]
            DF.loc[len(DF)] = L
            f0 = f1
            x0 = self.x0 + i * self.dx
            dfdx = self.derivada0(x0 = x0, f0 = f0)
            f_ev = self.EDO_Euler(x0=x0+self.dx, f0=f0)
            dfdx_ev = self.derivada0(x0 = x0+self.dx, f0 = f_ev)
            
            f1 = f0 + self.dx*(dfdx+dfdx_ev) / 2
        DF.to_csv(f"Resultados/{csv_name}")
        return DF
            
            
            
            


def dfdx(x0, dx, f0):
    return -2*(x0**3) + 12*(x0**2) -20*x0 + 8.5

Func = EDO(dfdx=dfdx, x0=0,f0=2 , dx=0.01)
DF = Func.predicion_EDO_Euler(pasos=400)

Func2 = EDO(dfdx=dfdx, x0=1,f0=2 , dx=0.01)
DF2 = Func2.predicion_EDO_Euler_Optimizado(pasos=200,csv_name="Resultado2 EDO_Euler.csv")

#Func3 = EDO(dfdx=dfdx, x0=0,f0=1 , dx=0.00001)
#DF3 = Func3.predicion_EDO_Euler(pasos=100000,csv_name="Resultado3 EDO_Euler.csv")
import matplotlib.pyplot as plt
import numpy as np

X = DF["x"]
Y = DF["f(x)"]
fig = plt.plot(X,Y)
X2 = DF2["x"]
Y2 = DF2["f(x)"]
fig = plt.plot(X2,Y2, color="red")

#DF3 = pd.read_csv("Resultados/Resultado3 EDO_Euler.csv")
#X3 = DF3["x"]
#Y3 = DF3["f(x)"]
#fig = plt.plot(X3,Y3, color="green")
plt.savefig("graf2.png")
