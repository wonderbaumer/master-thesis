import numpy as np

class multi_parameter():
    def __init__(self , t , particle):
        self.time = t
        self.B = particle.B
        self.epsilon = particle.epsilon
        self.delta = particle.delta

    def r(self):
        r0 = 1.0
        r10 = self.B / (3 * (1 - self.B)) * (self.time - np.sin(self.time))
        r01 = -self.B / (3 * (1 - self.B)) * (np.sin(self.time) + 5 * self.time)

        rtot = r0 + self.epsilon * r10 + self.delta * r01

        return rtot