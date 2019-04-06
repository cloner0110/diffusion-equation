from scipy.sparse import diags
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
def implicitEuler(dt,dy,t_max,y_max,alpha,V0,V1):
    s = alpha*dt/dy**2  
    y = np.arange(0,y_max+dy,dy) 
    t = np.arange(0,t_max+dt,dt)
    nt = len(t) 
    ny = len(y) 
    V = np.zeros((ny,)) 
    V[0] = V0 
    V[-1] = V1 
    A = diags([-s, 1+2*s, -s], [-1, 0, 1], shape=(ny-2, ny-2)).toarray() 
    for n in range(nt): 
        Vn = V ;B = Vn[1:-1] 
        B[0] = B[0]+s*V0 ;B[-1] = B[-1]+s*V1
        V[1:-1] = np.linalg.solve(A,B) 
    return y,t,V,s
dt = 0.01 
dy = 0.0005 
alpha = 2*10**(-4) 
y_max = 0.04 # in m
V0 = 10.0 
V1 = 0.0 
for time in np.linspace(0,1.0,10):
    y,t,V,s = implicitEuler(dt,dy,time,y_max,alpha,V0,V1)
    plt.plot(y,V,'k')
plt.show()
