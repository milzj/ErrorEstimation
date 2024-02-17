import numpy as np

class Experiment(object):

    def __init__(self):

        self._Alpha = [0.0, 1e-3, 1.0]
        self._N = [2**n for n in range(5, 10)]

    @property
    def Nref(self):
       return 4*self._N[-1]

    @property
    def Alpha(self):
        return self._Alpha

    @property
    def N(self):
        return self._N
