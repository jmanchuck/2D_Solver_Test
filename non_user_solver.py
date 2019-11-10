import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd

c = 2.99 * (10 ** 8)
mu = 1.26 * (10 ** (-6))
epsilon = 8.85 * (10 ** (-12))
size = 7
h = np.zeros((size, size))
ex = np.zeros((size + 1, size))
ey = np.zeros((size, size + 1))

dt = 1 / (10 * c)
dx = 1
steps = 400
time = 2 * dt
end_time = dt * steps  # run for 100 timesteps

sigma = dt * 6
step = 0

# Gaussian properties
sigma = dt * 30
pulseStart = 0
pulseEnd = 6 * sigma
pulseMid = (pulseStart + pulseEnd) / 2
maxima = (1 / (sigma * math.sqrt(2*math.pi)))

while time < end_time:

    h_prev = np.copy(h)
    ex_prev = np.copy(ex)
    ey_prev = np.copy(ey)

    # update h
    h = h_prev + (dt / (dx * mu)) * (ex[1:] - ex[:-1] - ey[:, 1:] + ey[1:, :-1])

    # override h bottom left with Gaussian
    if pulseStart < time < pulseEnd:
        h[2][2] = (1 / (sigma * math.sqrt(2 * math.pi))) * (math.exp(-((time - pulseMid) ** 2) / (2 * (sigma ** 2))))

    ex[1:-1, :] = ex_prev[1:-1, :] + (dt / (dx * epsilon)) * (h[:-1] - h[1:])
    ey[:, 1:-1] = ey_prev[:, 1:-1] - (dt / (dx * epsilon)) * (h[:, :-1] - h[:, 1:])

    time += dt
    step += 1
    print(step)

    # dfx = pd.DataFrame(ex)
    # print(dfx)
    df = pd.DataFrame(h)
    print(df)

    # plt.pcolormesh(h)
    # plt.pause(0.001)
    # plt.draw()


