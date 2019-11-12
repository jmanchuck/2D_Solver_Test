import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd
from time import process_time

# Physical constants
c = 2.99 * (10 ** 8)
mu = 1.26 * (10 ** (-6))
epsilon = 8.85 * (10 ** (-12))

h_np_arrs = []

# ALL USER INPUTS HERE
sigma_w = 1
sigma_t = 1 / sigma_w
omega_0 = 0
n_max = 1
s = 10
stability = 0.2

# improvements:
# - take into account epsilon max to determine how we will sample a wavelength (to get dx or ds)
# - user inputs real time simulation length in seconds, use dt to calculate how many steps


# Simulation constants, calculate from user input but stays constant throughout simulation
omega_max = omega_0 + 3 * sigma_w
lambda_min = math.pi * 2 * c / (n_max * omega_max)
l_x, l_y = 20 * lambda_min, 20 * lambda_min
ds = lambda_min / s
dt = ds * stability / c
steps = 5000
end_time = dt * steps
constant_eps = (dt / (ds * epsilon))
constant_mu = (dt / (ds * mu))

# Gaussian derivative properties
pulseMid = 3 * sigma_t  # the middle of the gaussian derivative (approximate)
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

# Simulation variables, these change throughout the simulation
time = 0
step = 0
size = int(l_x / ds)
h = np.zeros((size, size))
ex = np.zeros((size + 1, size))
ey = np.zeros((size, size + 1))

# animation plotting properties
frames = np.arange(0, end_time, dt)
fig, ax = plt.subplots()
# mesh = plt.pcolormesh(h)

# material properties
eps_arr = (dt / ds) * (np.ones((size, size)) / epsilon)
mu_arr = (dt / ds) * (np.ones((size, size)) / mu)
epsilon_relative = 10

# adding square material in top right
eps_arr[:(size // 3), 2 * (size // 3):] *= epsilon_relative
# mu_arr[:(size // 3), 2 * (size // 3):] *= 1

# reflecting box in the bottom left
row_upper = 2 * (size // 3)
row_lower = size - 1
col_upper = 0
col_lower = (size // 3)

def update(time):

    global h
    global step
    h_prev = h
    ex_prev = ex
    ey_prev = ey

    # update h
    h = h_prev + mu_arr * ((ex_prev[1:] - ex_prev[:-1]) - (ey_prev[:, 1:] - ey_prev[:, :-1]))

    # override h bottom left with Gaussian
    if 0 < time < pulseMid * 4:
        magnitude = (-time + pulseMid) * (1 / (sigma_t * math.sqrt(2 * math.pi))) * (
            math.exp(-((time - pulseMid) ** 2) / (2 * (sigma_t ** 2))))
        h[size//2][size//2] = magnitude
        # h[-2][-2] = magnitude

    h[row_upper:row_lower, col_upper:col_lower] = 0
    ex[1:-1, :] = ex_prev[1:-1, :] + eps_arr[1:, :] * (h[1:, :] - h[:-1, :])
    ey[:, 1:-1] = ey_prev[:, 1:-1] - eps_arr[:, 1:] * (h[:, 1:] - h[:, :-1])

    ex[row_upper:row_lower, col_upper:col_lower] = 0
    ey[row_upper:row_lower, col_upper:col_lower] = 0
    step += 1

    h_np_arrs.append(h)

    return h


# plot properties
im = plt.imshow(h, animated=True, vmax=maxima, vmin=-maxima, aspect='auto')
plt.set_cmap('seismic')
ax.xaxis.set_ticks_position('top') # the rest is the same
ax.xaxis.set_label_position('top')
fig.colorbar(im)


def animate(time):
    if step % 10 == 0:
        ax.set_title("Time Step = {}".format(step))
        # plt.savefig(str(step) + '.png')

    # if step % 100 == 0:
    #     print('.', end=' ')
    im.set_array(update(time))
    return im,


anim = animation.FuncAnimation(fig, animate, frames=np.arange(0, end_time, dt), interval=1, blit=False, repeat=False)
# anim.save(filename="movie.gif", fps=20)

# Set up formatting for the movie files
# Writer = animation.writers['pillow']
# writer = Writer(metadata=dict(artist='Me'), bitrate=1800)
# anim.save('movie.gif', writer=writer)
plt.show()
exit()