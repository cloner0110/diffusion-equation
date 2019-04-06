from scipy.sparse import diags
import numpy as np
from math import pi
import matplotlib.pyplot as plt
from scipy.sparse import diags
def k_calculator(T):
    k=T/((1.829*10**(-4))*T**(n))+0.0245
    return k
def implicitEuler(dt,dx,t_max,x_max,alpha):
    s = alpha*dt/dx**2  
    y = np.arange(0,x_max+dx,dx) 
    t = np.arange(0,t_max+dt,dt)
    dx = 0.001
    k=237 ; h=35 ; Acond=pi*0.005*0.005
    Aconv=pi*2*0.005*0.1
    nt = len(t) 
    ny = len(y) 
    V0 = 120+273 
    Vinf=273+15
    V = np.zeros((ny,)) 
    A = diags([-s, 1+2*s, -s], [-1, 0, 1], shape=(ny-2, ny-2)).toarray() 
    for n in range(nt): 
        Vn = V ;B = Vn[1:-1] 
        B[0] = B[0]+s*Vn[0] ;B[-1] = B[-1]+s*Vn[-1]
        V[1:-1] = np.linalg.solve(A,B) 
        V[0]=120+273
        t_ave=(Vn[-1]+Vn[-2])/2
        #switcher for 4th part of the Exercise , disable if you want to use constant K
        #c=dx*h*Aconv/(k_calculator(t_ave)*Acond)
        c=dx*h*Aconv/(k*Acond)
        V[-1]=((c*Vinf)+V[-2])/(1+c)
    return y,V
dt = 0.001 ;dx=0.001
alpha = 97.1*10**(-6) 
x_max = 0.1 
y = np.arange(0,x_max+dx,dx) 
nx = len(y)
plt.figure(figsize=(7,5))
V = np.zeros((nx,)) # initial condition
V[:] = 273+30
for time in np.linspace(0,100.0,5):
    y,V = implicitEuler(dt,dx,time,x_max,alpha)
    plt.plot(y,V,'k')
plt.show()

