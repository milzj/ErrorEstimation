import numpy as np

class Experiment(object):

    def __init__(self, name):

        if "LinearProblem"  in name:

            self._N = [2**n for n in range(8, 16)]
            self._Nref = 16*self._N[-1]

        else:

            self._N = [2**n for n in range(5, 10)]
            self._Nref = 4*self._N[-1]

        self._Alpha = [0.0, 1e-3, 1.0]


    @property
    def Nref(self):
       return self._Nref

    @property
    def Alpha(self):
        return self._Alpha

    @property
    def N(self):
        return self._N
