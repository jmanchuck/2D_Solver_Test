import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from time import process_time

c = 2.99 * (10 ** 8)
mu = 1.26 * (10 ** (-6))
epsilon = 8.85 * (10 ** (-12))

sigma_w = 1
sigma_t = 1 / sigma_w
omega_0 = 0
n_max = 1
s = 10
stability = 0.1

omega_max = omega_0 + 3 * sigma_w
lambda_min = math.pi * 2 * c / (n_max * omega_max)
l_x, l_y = 10 * lambda_min, 10 * lambda_min


ds = lambda_min / s

dt = ds * stability / (c)
steps = 1000
time = 0
end_time = dt * steps  # run for 200 timesteps

step = 0

# Gaussian properties
maxima = (1 / (sigma_t * math.sqrt(2*math.pi)))
pulseMid = 3 * sigma_t

time = 0
size = int(l_x / ds)
h = np.zeros((size, size))
ex = np.zeros((size + 1, size))
ey = np.zeros((size, size + 1))

fig, ax = plt.subplots()

mesh = ax.pcolormesh(h)

constant_e = (dt / (ds * epsilon))
constant_mu = (dt / (ds * mu))

def update(time):

    start = process_time()

    global h
    global step
    h_prev = h
    ex_prev = ex
    ey_prev = ey


    # update h
    h = h_prev + constant_mu * ((ex_prev[1:] - ex_prev[:-1]) - (ey_prev[:, 1:] - ey_prev[:, :-1]))

    # override h bottom left with Gaussian
    if 0 < time < 300 * dt:
        h[50][50] = (-time + pulseMid) * (1 / (sigma_t * math.sqrt(2 * math.pi))) * (
            math.exp(-((time - pulseMid) ** 2) / (2 * (sigma_t ** 2))))

    ex[1:-1, :] = ex_prev[1:-1, :] + constant_e * (h[1:, :] - h[:-1, :])
    ey[:, 1:-1] = ey_prev[:, 1:-1] - constant_e * (h[:, 1:] - h[:, :-1])

    # run = process_time() - start
    # print(process_time() - run)
    # print(time / dt)

    ax.set_title(time / dt)

    mesh = ax.pcolormesh(h)
    step += 1

    t = process_time()
    print(t - start)

    return mesh


def animate(time):

    return update(time)


anim = animation.FuncAnimation(fig, func=animate, frames=np.arange(0, end_time, dt), interval=1, repeat=False, cache_frame_data=False)

plt.show()
exit()