import numpy as np
import math

c = 299792458
mu = 4 * math.pi * (10 ** (-7))
epsilon = math.sqrt(mu * (c ** 2))

h = np.zeros((50, 50))
ex = np.zeros((50, 50))
ey = np.zeros((50, 50))

h_prev = h
ex_prev = ex
ey_prev = ey

dt = 1 / (c + 1)
dx = 1
steps = 100
step = 0
end_time = dt * steps  # run for 100 timesteps

sigma = dt * 2
start_pulse_step = 0
end_pulse_step = 6


def gaussian_excitation(start, end, sigma, time):
    t_mean = (start + end) // 2
    excitation = (1 / (sigma * math.sqrt(2 * math.pi))) * math.exp(-((time - t_mean) ** 2) / (2 * (sigma ** 2)))
    return excitation

for i in range(10):
    print(gaussian_excitation(0, 8, 4, i))


for i in range(2, steps):
    if step <= end_pulse_step:
        h[24][24] = gaussian_excitation(start_pulse_step * dt, end_pulse_step * dt, sigma, i * dt)

    h_prev = h
    h = np.zeros((50, 50))
    for j in range(50):
        for k in range(50):
            if k < 49:
                ex_down = ex[j][k+1]
            else:
                ex_down = 0
            if j < 49:
                ey_right = ey[j+1][k]
            else:
                ey_right = 0
            h[j][k] = h_prev[j][k] + (dt / (mu * dx)) * (ex_down - ex[j][k] - ey_right + ey[j][k])

    print(h)