import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd

# constants
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

    def place(self, i, j, magnitude):
        """
        Args:
            i: index of row
            j: index of col
            magnitude: magnitude of the excitation
        """

        # if out of range
        if i < 0 or i > self.x_size - 1:
            raise ValueError("Tried to place at row " + i)
        if j < 0 or j > self.y_size - 1:
            raise ValueError("Tried to place at col " + j)

        # update the current array at row col
        self.current[i][j] = magnitude

    def currentCopy(self):
        """

        Returns: deep copy of the current array to prevent assigning by reference

        """

        return np.copy(self.current)

    def getCurrent(self, i, j):
        """
        Args:
            i: index of row
            j: index of column

        Returns: the field magnitude at that point, if outside of edges it returns 0 (reflect condition)
        """
        if i < 0 or i >= self.y_size or j < 0 or j >= self.x_size:
            return 0
        else:
            return self.current[i][j]


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

        # point source variables
        self.pulse = False
        self.sigma = None
        self.pulseStart = None
        self.pulseMax = None
        self.pulseMiddle = None
        self.pulseLocation = size // 2  # set the pulse location to be single point at middle of array

        # init the above variables for a gaussian point source with sigma, pulse start time, pulse end time
        self.gaussianSource(sigma=2 * self.dt, pulseStart=0, pulseMax=5 * self.dt)

        # catch user error for mesh input
        if c * dt > dx:
            raise ValueError("Time step ", dt, "is too large compared to distance step ", dx)

    def gaussianSource(self, sigma, pulseStart, pulseMax):
        self.sigma = sigma
        self.pulseStart = pulseStart
        self.pulseMax = pulseMax
        self.pulseMiddle = (self.pulseStart + self.pulseMax) / 2

    def gaussian(self):
        magnitude = (1 / (self.sigma * math.sqrt(2 * math.pi))) * math.exp(
            -((self.time - self.pulseMiddle) ** 2) / (2 * (self.sigma ** 2)))
        # print(magnitude)
        return magnitude

    def update(self):
        """
        Applies all 3 update equations on H, Ex and Ey in order
        """
        self.h.previous = self.h.currentCopy()
        self.ex.previous = self.ex.currentCopy()
        self.ey.previous = self.ey.currentCopy()

        # if time passed pulse start, pulse boolean is true and update method will use override method
        # otherwise turn off override method
        if self.time >= self.pulseStart and not self.pulse:
            self.pulse = True
        if self.time > self.pulseMax and self.pulse:
            self.pulse = False

        for i in range(self.h.y_size):
            for j in range(self.h.x_size):
                self.update_h(i, j)

        for i in range(self.ex.y_size):
            for j in range(self.ex.x_size):
                self.update_x(i, j)
                self.update_y(i, j)

        self.time += self.dt

        # saving the image with name
        if self.time // self.dt < 10:
            plt.imsave("Figs/00{}.png".format(self.time // self.dt), self.h.current)
        elif self.time // self.dt < 100:
            plt.imsave("Figs/0{}.png".format(self.time // self.dt), self.h.current)
        elif self.time // self.dt >= 100:
            plt.imsave("Figs/{}.png".format(self.time // self.dt), self.h.current)

        df = pd.DataFrame(self.h.current)
        print(df)

        df1 = pd.DataFrame(self.ex.current)
        print(df1)

        df2 = pd.DataFrame(self.ey.current)
        print(df2)

        # if self.time == 0:
        #     mode = 'w'
        # else:
        #     mode = 'a'
        # df.to_csv('h_' + str(self.size) + '.csv', mode=mode, index=False, header=False)

    def update_h(self, i, j):
        """
        Args:
            i: row index
            j: column index
        Applies update equation for a specific point
        """

        # if within time frame where we want to have a pulse, h_next value gets override by pulse value
        if self.pulse and i == self.pulseLocation and j == self.pulseLocation:
            h_next = self.gaussian()

        else:
            h_next = self.h.previous[i][j] + (self.dt / (self.dx * mu)) * (
                        (self.ex.getCurrent(i, j + 1) - self.ex.getCurrent(i, j - 1)) -
                        (self.ey.getCurrent(i + 1, j) - self.ey.getCurrent(i - 1, j)))
        self.h.place(j, i, h_next)

    def update_x(self, i, j):
        e_next = self.ex.previous[i][j] + (self.dt / (self.dx * epsilon)) * (
                    self.h.getCurrent(i, j + 1) - self.h.getCurrent(i, j - 1))
        self.ex.place(j, i, e_next)

    def update_y(self, i, j):
        e_next = self.ey.previous[i][j] - (self.dt / (self.dx * epsilon)) * (
                    self.h.getCurrent(i + 1, j) - self.h.getCurrent(i - 1, j))
        self.ey.place(j, i, e_next)


def test():
    # testing
    dt = 1 / (2 * c)
    dx = 1
    time = 4 * dt  # simulate for _ time steps
    size = 5  # size of square array
    solver = Solver(size=size, time=time, dt=dt, dx=dx)
    # print(solver.h.current)

    while solver.time < solver.endTime:
        solver.update()
        # print(solver.h.current)
        # print(solver.ex.current)
        # print(solver.ey.current)

    print("making gif...")


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

    imageio.mimsave("movie.gif", images, fps=5)


if __name__ == "__main__":
    test()
    # make_existing_gif()
    print("finished")
