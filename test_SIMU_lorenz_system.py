import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def lorenz(t,x): #defining the lorenz equation to be solved
    Beta = np.array([10, 28, 8 / 3])  # chaotic values
    dx = np.array([Beta[0]*(x[1]-x[0]),
                   x[0]*(Beta[1]-x[2]) - x[1],
                   x[0]*x[1] - Beta[2]*x[2]])
    return dx

x0 = np.array([0, 1, 20]) #initial condition
dt = 0.001 #resolution
timeSpan = np.arange(dt, 50+dt, dt)
#timeSpan = np.expand_dims(timeSpan, axis=0)

sol = solve_ivp(lorenz,[timeSpan[0],timeSpan[-1]],x0,method='RK45',t_eval=timeSpan,rtol=1e-12)
result = sol.y

ax = plt.axes(projection='3d')
ax.plot3D(result[0,:],result[1,:],result[2,:])
