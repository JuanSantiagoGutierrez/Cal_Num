import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

csv = "Resultados/Regula-Fansi_2025-04-15.csv"
df = pd.read_csv(csv)


resolucion = 25
def f(x):
    return x**3 + x*2 +1

long = resolucion*int(df["a"][0])
long2 = resolucion*int(df["b"][0])
x = np.array(range(long,long2))/resolucion
plt.plot(x,f((np.linspace(df["a"][0], df["b"][0], len(x)))))
for i in range(10):
    p_x = [df["a"][i],df["b"][i]]
    p_y = [df["f(a)"][i],df["f(b)"][i]]

    plt.plot(p_x,p_y,"*", color = "red" if i == 0 else "green")
    coeff = np.polyfit(p_x,p_y,1)
    recta = np.poly1d(coeff)
    
    y = np.zeros(len(x))

    for j in range(len(x)):
        y[j] = np.polyval(recta, x[j])
    plt.plot(x,y)
    
    plt.plot(x,np.zeros(len(x)))
    
    
    plt.savefig(f"graf.png")
    from time import sleep
    sleep(4)
    


