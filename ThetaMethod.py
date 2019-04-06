from scipy.sparse import diags
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from math import pi
def k_calculator(T):
    k=T/((1.829*10**(-4))*T**(2))+0.0245
    return k
def theta_Method(dx,h,k,nx,dt,nt,D,V,ntout,theta):
    Vout = [] 
    V0 = 120+273 
    Vinf=273+15
    s = D*dt/dx**2
    A = diags([-theta*s, 1+(2*theta*s), -theta*s], [-1, 0, 1], 
          shape=(nx-2, nx-2)).toarray() 
    B1 =diags([theta*s, 1-(2*theta*s), theta*s],[-1, 0, 1], shape=(nx-2, nx-2)).toarray()
    for n in range(1,nt): 
        Vn = V
        if theta==1:
            Acond=pi*0.005*0.005
            Aconv=pi*2*0.005*0.1
            t_ave=(Vn[-1]+Vn[-2])/2
            #switcher for 4th part of the Exercise , disable if you want to use constant K
            #c=dx*h*Aconv/(k_calculator(t_ave)*Acond)
            c=dx*h*Aconv/(k*Acond)
            B = Vn[1:-1]
            B[0] = B[0]+s*Vn[0]
            B[-1] = B[-1]+s*Vn[-1]
        else:
            Acond=pi*0.005*0.005
            Aconv=pi*2*0.005*0.1
            c=dx*h*Aconv/(k*Acond)
            B = np.dot(Vn[1:-1],B1) 
            B[0] = B[0]+(1-theta)*s*(V0+V0)
            B[-1] = B[-1]+(1-theta)*s*(Vn[-1]+Vn[-1])
        V[1:-1] = np.linalg.solve(A,B)
        V[0]=120+273
        V[-1]=((c*Vinf)+V[-2])/(1+c)
        if n % int(nt/float(ntout)) == 0 or n==nt-1:
            Vout.append(V.copy()) 
    return Vout
dt = 0.001 
dx = 0.001
k=237 ; h=35
theta=1
alpha = 97.1*10**(-6) 
y_max = 0.1 
y = np.arange(0,y_max+dx,dx) 
nx = len(y)
nt = 9000
plt.figure(figsize=(7,5))
V = np.zeros((nx,)) # initial condition
V[:] = 273+30
Vout= theta_Method(dx,h,k,nx,dt,nt,alpha,V,5,theta)
for V in Vout:
    plt.plot(y,V)
    plt.hold
plt.show()
