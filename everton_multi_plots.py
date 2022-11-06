import numpy as np
import matplotlib.pyplot as plt

fig, axi = plt.subplots()
for i in range(0,50):
    t = np.linspace(0,(i*np.pi),1000)
    axi.cla()
    axi.plot(t, np.cos(t))
    axi.set_xlabel('t')
    axi.set_ylabel('cos(t)')
    axi.set_title(f'iteração = {i}')
    plt.pause(0.001)

