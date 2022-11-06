import numpy as np
import matplotlib.pyplot as plt

def FDM_we(I,V,f,c,L,dt,C,T):
    Nt = int(round(T/dt)) #number of points in time discretization
    t = np.linspace(0,Nt*dt,Nt+1) #temporal discretization
    dx_t = dt*c/float(C) #spatial discretization from Courrant number
    Nx = int(round(L/dx_t)) #number of points in space discretization
    x = np.linspace(0,L,Nx+1) #spatial discretization
    C2 = C**2 #Courrant number squared

    #compaptibility with x and t
    dx = x[1]-x[0]
    dt = t[1] - t[0]

    if f is None or f == 0:
        f = lambda x, t: 0
    if V is None or V == 0:
        V = lambda x: 0

    u = np.zeros(Nx+1) #u[n+1]
    u_n = np.zeros(Nx+1) #u[n]
    u_nm1 = np.zeros(Nx+1) #u[n-1]

    for i in range(0,Nx+1):
        u_n[i] = I(x[i]) #initial condition in u_n
        axi.cla()
        axi.plot(x, u)
        plt.pause(0.001)

    n = 0
    for i2 in range(1,Nx):
        u[i2] = u_n[i2] + dt*V(x[i]) + 0.5*C2*(u_n[i2-1] - 2*u_n[i2] + u_n[i2+1]) + 0.5*dt**2*f(x[i], t[n]) #cálculo de

        axi.cla()
        axi.plot(x, u)
        plt.pause(0.001)
    u[0] = 0
    u[Nx] = 0

    #switch variables before next step
    u_nm1[:] = u_n
    u_n[:] = u

    for n in range(1,Nt):
        for i3 in range(1,Nx):
            u[i3] = -u_nm1[i3] + 2*u_n[i3] + C2*(u_n[i3-1] - 2*u_n[i3] + u_n[i3+1]) + dt**2*f(x[i], t[n])
            print(dt**2*f(x[i], t[n]))

        u[0] = 0
        u[Nx] = 0

        axi.cla()
        axi.plot(x, u)
        plt.pause(0.001)

        #switch variables before nect step
        u_nm1[:] = u_n
        u_n[:] = u

    return u,x,t

def plot_u(I,V,f,c,L,dt,C,T):
    u, x, t = FDM_we(I,V,f,c, L, dt, C, T)
    #plt.plot(x,u,'r-')



#Running a case: vibration of a string
def guitar(C):
    L = 0.75
    x0 = 0.8*L
    a = 0.005
    freq = 440
    wavelength = 2*L
    c = freq*wavelength
    omega = 2*np.pi*freq
    #num_periods = 1
    T = 2*np.pi/omega
    dt = L/50./c
    V=0
    f=0

    def I(x):
        if x<x0:
            return a*x/x0
        else:
            return a/(L-x0)*(L-x)

    return plot_u(I,V,f,c,L,dt,C,T)

fig, axi = plt.subplots()
C = 0.75
guitar(C)








