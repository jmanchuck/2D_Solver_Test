import numpy as np

# ADDDED STUFF

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
        self.current = np.zeros([self.y_size, self.x_size])
        self.previous = np.zeros([self.y_size, self.x_size])
        self.time = 0
        self.half = half
        # TODO
        ## add delta T delta X

    def place(self, row, col, magnitude):
        """
        Args:
            x: x coordinate in plane
            y: y coordinate in plane
            magnitude: magnitude of the excitation
        """

        # if out of range
        if row < 0 or row > self.x_size - 1:
            raise Exception
        if col < 0 or col > self.y_size - 1:
            raise Exception
        self.current[row][col] = magnitude

    def currentCopy(self):
        """

        Returns: deep copy of the current array to prevent assigning by reference

        """

        # make into np array
        new = []
        for i in self.current:
            new.append(i)
        return new

    def increment(self):
        # increment by delta T
        self.time += 1

    def getCurrent(self, i, j):
        """

        Args:
            i: index
            j: index

        Returns: the field magnitude at that point, if outside of edges it returns 0 (reflect condition)

        """
        if i < 0 or i >= self.y_size or j < 0 or j >= self.x_size:
            return 0
        else:
            return self.current[i][j]


class Solver:

    def __init__(self, size, time, deltaT, deltaX):
        self.ex = Matrix((size, size))
        self.ey = Matrix((size, size))
        self.h = Matrix((size, size))
        self.deltaT = deltaT
        self.deltaX = deltaX
        self.endTime = time
        self.time = 0
        self.const = 1

    def update(self):
        """
        Applies all 3 update equations on H, Ex and Ey in order
        """
        # temp_h = self.h.currentCopy()
        # temp_x = self.ex.currentCopy()
        # temp_y = self.ey.currentCopy()

        self.h.previous = self.h.currentCopy()
        self.ex.previous = self.ex.currentCopy()
        self.ey.previous = self.ey.currentCopy()

        for i in range(self.h.y_size):
            for j in range(self.h.x_size):
                self.update_h(i, j)

        for i in range(self.mesh):
            for j in range(self.mesh):
                self.update_x(i, j)
                self.update_y(i, j)

        self.time += 1

    def update_h(self, i, j):
        """
        Args:
            i: row index
            j: column index
        Applies update equation for a specific point
        """
        h_next = self.h.previous[i][j] + (1 / mu) * ((self.ex.getCurrent(i, j + 1) - self.ex.getCurrent(i, j - 1)) -
                                                     (self.ey.getCurrent(i + 1, j) - self.ey.getCurrent(i - 1, j)))
        self.h.place(j, i, h_next)
        self.h.increment()

    def update_x(self, i, j):
        e_next = self.ex.previous[i][j] + (1 / epsilon) * (self.h.getCurrent(i, j + 1) - self.h.getCurrent(i, j - 1))
        self.ex.place(j, i, e_next)
        self.ex.increment()

    def update_y(self, i, j):
        e_next = self.ey.previous[i][j] - (1 / epsilon) * (self.h.getCurrent(i + 1, j) - self.h.getCurrent(i - 1, j))
        self.ey.place(j, i, e_next)
        self.ey.increment()


def main():
    solver = Solver(5, 10)
    solver.h.place(6, 6, 0.01)
    print(solver.h.current)
    solver.ex.place(6, 6, 0.1)
    solver.ey.place(6, 6, 0.1)

    while solver.time < solver.endTime:
        solver.update()
        print(solver.h.current)
        # print(solver.ex.current)
        # print(solver.ey.current)


if __name__ == "__main__":
    main()
