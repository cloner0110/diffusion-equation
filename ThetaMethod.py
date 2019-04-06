from scipy.sparse import diags
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
def theta_Method(dy,ny,dt,nt,D,V,ntout,theta):
    Vout = [] 
    V0 = 200 
    V1 = 200 
    s = D*dt/dy**2  
    A = diags([-theta*s, 1+(2*theta*s), -theta*s], [-1, 0, 1], 
          shape=(ny-2, ny-2)).toarray() 
    B1 =diags([theta*s, 1-(2*theta*s), theta*s],[-1, 0, 1], shape=(ny-2, ny-2)).toarray()
    for n in range(1,nt): 
        Vn = V
        if theta==1:
            B = Vn[1:-1]
            B[0] = B[0]+s*V0
            B[-1] = B[-1]+s*V1
        else:
            B = np.dot(Vn[1:-1],B1) 
            B[0] = B[0]+(1-theta)*s*(V0+V0)
            B[-1] = B[-1]+(1-theta)*s*(V1+V1)
        V[1:-1] = np.linalg.solve(A,B)
        if n % int(nt/float(ntout)) == 0 or n==nt-1:
            Vout.append(V.copy()) 
    return Vout,s
dt = 0.001 
dy = 0.001 
theta=0.5
alpha = 97.1*10**(-6) 
y_max = 0.04 
y = np.arange(0,y_max+dy,dy) 
ny = len(y)
nt = 3000
plt.figure(figsize=(7,5))
V = np.zeros((ny,)) # initial condition
V[0] = 0
Vout,s = theta_Method(dy,ny,dt,nt,alpha,V,10,theta)
for V in Vout:
    plt.plot(y,V,'k')
plt.show()
