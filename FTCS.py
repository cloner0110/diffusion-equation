import numpy as np 
from math import pi
import matplotlib.pyplot as plt 
alpha=97.1*10**(-6)
uinitial=273+30
k=237 ; h=35 ; Acond=pi*0.005*0.005
Aconv=pi*2*0.005*0.1
Nx=10 
dx=0.1/(Nx-1)
dt=0.25
Vinf=15+273
timestep=1000
c=dx*h*Aconv/(k*Acond)
u_tmp=np.zeros(Nx)
u=np.zeros(Nx)
x=np.linspace(0,0.1,Nx)
u[:]=uinitial
s = alpha*dt/dx**2
for j in range(600):
    u[0]=120+273
    u[-1]=((c*Vinf)+u_tmp[-2])/(1+c)
    for i in range(1,Nx-1):
        u[i]=u_tmp[i]+s*(u_tmp[i-1]-(2*u_tmp[i])+u_tmp[i+1])
    u_tmp=u
    if j%30==0:
        plt.plot(x,u)
plt.show()