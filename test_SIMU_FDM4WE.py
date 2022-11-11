#Implementation of a FDM simulator based on the 1D wave equation model
#Based on the book:
#Finite difference methods for wave equations
#by Hans Petter Langtangen and Svein Linge

import numpy as np
import matplotlib.pyplot as plt

def FDM_we(I,V,f,c,L,dt,C,T,new_velo,veloc=None): #constant velocity by default
    Nt = int(round(T/dt)) #number of points in time discretization
    t = np.linspace(0,Nt*dt,Nt+1) #temporal discretization
    dx_t = dt*c/float(C) #spatial discretization from Courrant number
    dt2 = dt * dt #dt squared
    Nx = int(round(L/dx_t)) #number of points in space discretization
    x = np.linspace(0,L,Nx+1) #spatial discretization

    if veloc is None:
        print(f'Velocity in the medium = {c} m/s')
        C2 = C**2 #Courrant number squared
    else:
        print(f'Velocity in the first medium = {c} m/s')
        print(f'Velocity in the second medium = {new_velo} m/s')
        C2 = (dt/dx_t)**2 #Courrant number squared considering variable velocity
        xx,zz = np.meshgrid(np.arange(Nx),0)
        q = np.zeros(Nx)
        idx_newvelo = int(np.round(2*Nx/3))
        q[:xx[0,idx_newvelo]] = c
        q[xx[0,idx_newvelo]:] = new_velo
        print(q)

    #compaptibility with x and t
    dx = x[1] - x[0]
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
        axi.set_ylim([-0.005,0.005])
        axi.set_xlim([0,L])
        plt.pause(0.001)

    # Border conditions
    # n = 0
    # for i2 in range(1,Nx):
    #     u[i2] = u_n[i2] + dt*V(x[i]) + 0.5*C2*(u_n[i2-1] - 2*u_n[i2] + u_n[i2+1]) + 0.5*dt**2*f(x[i], t[n]) #cálculo de
    #
    #     axi.cla()
    #     axi.plot(x, u)
    #     plt.pause(0.001)
    # u[0] = 0
    # u[Nx] = 0
    #
    # #switch variables before next step
    # u_nm1[:] = u_n
    # u_n[:] = u

    #Neumann border conditions
    for i4 in range(0, Nx + 1):
        if i4<Nx:
            ip1 = i4+1
        else:
            ip1 = i4-1
        if i4>0:
            im1 = i4-1
        else:
            im1 = i4+1

        u[i4] = u_n[i4] + dt2 * V(x[i4]) + 0.5 * C2 * (u_n[im1] - 2 * u_n[i4] + u_n[ip1]) + 0.5 * dt2 * f(x[i4], t[0])

        axi.cla()
        axi.plot(x, u)
        axi.set_ylim([-0.005, 0.005])
        axi.set_xlim([0, L])
        plt.pause(0.001)

    # Update data structures for next step
    u_nm1[:] = u_n
    u_n[:] = u  # safe, but slower

    for n in range(1,Nt):
        if veloc is None:
            for i3 in range(1, Nx+1):
                if i3 < Nx:
                    ip2 = i3 + 1
                else:
                    ip2 = i3 - 1
                if i3 > 0:
                    im2 = i3 - 1
                else:
                    im2 = i3 + 1

                u[i3] = -u_nm1[i3] + 2*u_n[i3] + C2*(u_n[im2] - 2*u_n[i3] + u_n[ip2]) + dt2*f(x[i], t[n])
        else:
            for i3 in range(1, Nx - 1):
                if i3 < Nx:
                    ip2 = i3 + 1
                else:
                    ip2 = i3 - 1
                if i3 > 0:
                    im2 = i3 - 1
                else:
                    im2 = i3 + 1

                u[i3] = -u_nm1[i3] + 2 * u_n[i3] + C2 * (0.5*(q[i3]+q[ip2])*(u_n[ip2] - u_n[i3]) - \
                0.5*(q[i3] + q[im2])*(u_n[i3] - u_n[im2]))+ dt2 * f(x[i], t[n])

        u[0] = 0
        u[Nx] = 0

        axi.cla()
        axi.plot(x, u)
        axi.set_ylim([-0.005, 0.005])
        axi.set_xlim([0, L])
        plt.pause(0.001)

        #switch variables before nect step
        u_nm1[:] = u_n
        u_n[:] = u

    return u,x,t

def run(I,V,f,c,L,dt,C,T):
    new_velo = 1000  # new velocity in m/s
    u, x, t = FDM_we(I,V,f,c, L, dt, C, T,new_velo)
    #veloc=
#Running a case: vibration of a string
fig, axi = plt.subplots()
C = 0.75
L = 0.75
x0 = 0.8*L
a = 0.005
freq = 440
wavelength = 2*L
c = freq*wavelength
omega = 2*np.pi*freq
num_periods = 10
T = 2*np.pi/omega*num_periods
dt = L/50./c
V=0
f=0

def I(x):
    if x<x0:
        return a*x/x0
    else:
        return a/(L-x0)*(L-x)

run(I,V,f,c,L,dt,C,T)









