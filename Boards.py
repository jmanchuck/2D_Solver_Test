import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate

global c
c = 2.99 * (10 ** 8)

global epsilon
epsilon = 8.85 * (10 ** (-12))

global mu
mu = 1.26 * (10 ** (-6))


class Matrix:

    def __init__(self, dimensions, half=False):
        self.x_size = dimensions[1]
        self.y_size = dimensions[0]
        self.half = half

        self.current = np.zeros([self.y_size, self.x_size])
        self.previous = np.zeros([self.y_size, self.x_size])

    def updatePrev(self):
        """

        Returns: deep copy of the current array to prevent assigning by reference

        """

        self.previous = np.copy(self.current)

class Solver:

    def __init__(self, size, time, dt, dx):

        # matrices
        self.ex = Matrix((size, size))
        self.ey = Matrix((size, size))
        self.h = Matrix((size, size))

        # size
        self.size = size

        # mesh
        self.dt = dt
        self.dx = dx

        # length of simulation
        self.endTime = time

        # current time
        self.time = 0
        self.step = 0

        # point source variables
        self.pulse = False
        self.sigma = None
        self.pulseStart = None
        self.pulseMax = None
        self.pulseMiddle = None
        self.gaussianMax = None
        self.pulseLocation = size // 2  # set the pulse location to be single point at middle of array

        # init the above variables for a gaussian point source with sigma, pulse start time, pulse end time
        self.gaussianSource(sigma=40 * self.dt, pulseStart=0, pulseMax=240 * self.dt)

        self.cmap = plt.get_cmap('viridis')

        # catch user error for mesh input
        if c * dt > dx:
            raise ValueError("Time step ", dt, "is too large compared to distance step ", dx)

    def gaussianSource(self, sigma, pulseStart, pulseMax):
        self.sigma = sigma
        self.pulseStart = pulseStart
        self.pulseMax = pulseMax
        self.pulseMiddle = (self.pulseStart + self.pulseMax) / 2
        self.gaussianMax = (1 / (self.sigma * math.sqrt(2*math.pi)))

        # user defines sigma_time or sigma_freq
        # ft gaussian gives 1/sigma, def f_max = 3(1/sigma)

    def gaussian(self):
        magnitude = (1 / (self.sigma * math.sqrt(2 * math.pi))) * math.exp(
            -((self.time - self.pulseMiddle) ** 2) / (2 * (self.sigma ** 2)))
        print(magnitude)
        return magnitude

    def update(self):
        """
        Applies all 3 update equations on H, Ex and Ey in order
        """
        self.h.updatePrev()
        self.ex.updatePrev()
        self.ey.updatePrev()

        # if time passed pulse start, pulse boolean is true and update method will use override method
        # otherwise turn off override method
        if self.time >= self.pulseStart and not self.pulse:
            self.pulse = True
        if self.time > self.pulseMax and self.pulse:
            self.pulse = False

        self.h.current[1:-1, 1:-1] = self.h.previous[1:-1, 1:-1] + (self.dt / (self.dx * mu)) * (
                    self.ex.current[2:, 1:-1] - self.ex.current[:-2, 1:-1] - self.ey.current[1:-1, 2:] + self.ey.current[1:-1, :-2])

        if self.pulse:
            self.h.current[self.pulseLocation][self.pulseLocation] = self.gaussian()

        self.ex.current[1:-1, 1:-1] = self.ex.previous[1:-1, 1:-1] + (self.dt / (self.dx * epsilon)) * (self.h.current[1:-1, 2:] - self.h.current[1:-1, :-2])
        self.ey.current[1:-1, 1:-1] = self.ey.previous[1:-1, 1:-1] - (self.dt / (self.dx * epsilon)) * (self.h.current[2:, 1:-1] - self.h.current[:-2, 1:-1])

        self.time += self.dt
        self.step += 1

        # df = pd.DataFrame(self.h.current)
        # print(df)
        # df1 = pd.DataFrame(self.ex.current)
        # print(df1)
        #
        # df2 = pd.DataFrame(self.ey.current)
        # print(df2)

        # plt.pcolormesh(test, vmin=1, vmax=self.gaussianMax, shading="gouraud", cmap=self.cmap)
        # plt.pause(0.001)
        # plt.draw()
        print(self.step)


def test():
    # testing
    dt = 1 / (10 * c)
    dx = 1
    time = 500 * dt  # simulate for _ time steps
    size = 50  # size of square array
    solver = Solver(size=size, time=time, dt=dt, dx=dx)
    # print(solver.h.current)

    while solver.time < solver.endTime:
        solver.update()


def make_existing_gif():
    import os
    import imageio

    directory = "Figs/"
    images = []
    file_names = sorted(os.listdir(directory))
    for file in file_names:
        if file.endswith(".png"):
            file_path = os.path.join(directory, file)
            images.append(imageio.imread(file_path))

    imageio.mimsave("movie.gif", images, fps=10)


if __name__ == "__main__":
    test()
    # make_existing_gif()
    print("finished")
