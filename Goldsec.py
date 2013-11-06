# -*- coding: utf-8 -*-
from numpy import *
from math  import *
from Bracket import bracket

class Goldsec:
    TINY = 1.0e-12
    def __init__(self, maxit, eps, func):
        self.MAXIT = maxit
        self.EPS   = eps
        self.func  = func
        self.df    = lambda x : self.dfunc(x)

    def bracket(self, a, b):
        (self.a, self.b) = bracket(a, b, self.df)

    def dfunc(self, x):
        return (self.func(x+Goldsec.TINY) - self.func(x)) / Goldsec.TINY

    def minimize(self):
        x1 = self.a
        x3 = self.b
        GOLD = (1.0 + sqrt(5.0)) / 2.0
        x2 = (1.0 * x1 + GOLD * x3) / (1.0 + GOLD)
        for it in range(self.MAXIT):
            l = [x1, x2, x3]
            l.sort()
            (x1, x2, x3) = l
            if(abs(x1 - x3) < Goldsec.TINY):
                break

            x4  = x1 + (x3 - x2)
            print('[%03d] %f %f %f %f' % (it+1, x1, x2, x3, x4))
            f2 = self.func(x2)
            f4 = self.func(x4)
            if x2 < x4:
                if f2 < f4:
                    x3 = x4
                else:
                    x1 = x2
                    x2 = x4
            else:
                if f2 < f4:
                    x1 = x4
                else:
                    x3 = x2
                    x2 = x4

        self.fmin = self.func(x1)
        self.xmin = x1
        return self.xmin

if __name__=='__main__':
    f = lambda x : (x - 2) ** 2
    gc = Goldsec(100, 1.0e-8, f)
    gc.bracket(00.0, 1.0)
    gc.minimize()