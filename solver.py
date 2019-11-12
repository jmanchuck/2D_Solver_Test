import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd
from time import process_time

c = 2.99 * (10 ** 8)
mu = 1.26 * (10 ** (-6))
epsilon = 8.85 * (10 ** (-12))

sigma_w = 1
sigma_t = 1 / sigma_w
omega_0 = 0
n_max = 1
s = 10
stability = 0.2

omega_max = omega_0 + 3 * sigma_w
lambda_min = math.pi * 2 * c / (n_max * omega_max)
l_x, l_y = 20 * lambda_min, 20 * lambda_min

ds = lambda_min / s

dt = ds * stability / c
steps = 5000
time = 0
end_time = dt * steps  # run for 200 timesteps

step = 0

# Gaussian derivative properties
pulseMid = 3 * sigma_t
def gaussian_derivative_maxima():
    t = 0
    current = 0
    prev = 0
    while t < pulseMid:
        current = (-t + pulseMid) * (1 / (sigma_t * math.sqrt(2 * math.pi))) * (
        math.exp(-((t - pulseMid) ** 2) / (2 * (sigma_t ** 2))))
        if current < prev:
            return prev
        t += dt
        prev = current
    return None
maxima = gaussian_derivative_maxima()
print(maxima)

time = 0
size = int(l_x / ds)
h = np.zeros((size, size))
ex = np.zeros((size + 1, size))
ey = np.zeros((size, size + 1))

print(end_time)
frames = np.arange(0, end_time, dt)

fig, ax = plt.subplots()
# mesh = plt.pcolormesh(h)

print('maximum', maxima)

constant_e = (dt / (ds * epsilon))
constant_mu = (dt / (ds * mu))


def update(time):

    global h
    global step
    h_prev = h
    ex_prev = ex
    ey_prev = ey

    # update h
    h = h_prev + constant_mu * ((ex_prev[1:] - ex_prev[:-1]) - (ey_prev[:, 1:] - ey_prev[:, :-1]))

    # override h bottom left with Gaussian
    if 0 < time < pulseMid * 3:
        h[50][50] = (-time + pulseMid) * (1 / (sigma_t * math.sqrt(2 * math.pi))) * (
            math.exp(-((time - pulseMid) ** 2) / (2 * (sigma_t ** 2))))
    ex[1:-1, :] = ex_prev[1:-1, :] + constant_e * (h[1:, :] - h[:-1, :])
    ey[:, 1:-1] = ey_prev[:, 1:-1] - constant_e * (h[:, 1:] - h[:, :-1])

    ax.set_title(step)

    return h


# plot properties
im = plt.imshow(h, animated=True, vmax=maxima, vmin=-maxima, aspect='auto')
plt.set_cmap('seismic')
ax.xaxis.set_ticks_position('top') # the rest is the same
ax.xaxis.set_label_position('top')
fig.colorbar(im)


def animate(time):
    start = process_time()
    im.set_array(update(time))
    print(process_time() - start)
    return im,


anim = animation.FuncAnimation(fig, animate, frames=np.arange(0, end_time, dt), interval=1, blit=True, repeat=False)

plt.show()
exit()
