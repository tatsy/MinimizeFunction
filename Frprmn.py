# -*- coding: utf-8 -*-
from math import *
from numpy import *
from Linemethod import Linemethod

class Frprmn(Linemethod):
    def __init__(self, maxit, eps, p, xi, func):
        Linemethod.__init__(self, maxit, eps, p, xi, func)

    def minimize(self):
        g  = zeros(self.n)
        h  = zeros(self.n)
        fp = self.func(self.p)
        self.xi = self.grad(p)
        for j in range(self.n):
            g[j] = -self.xi[j]
            self.xi[j] = h[j] = g[j]

        fret = 0.0
        for it in range(self.MAXIT):
            try:
                fret = self.linmin()
            except Exception, msg:
                print(msg)

            if (2.0 * abs(fret - fp)) <= (self.EPS * (abs(fret) + abs(fp) + Frprmn.TINY)):
                return self.p

            fp = fret
            self.xi = self.grad(self.p)
            test = 0.0
            den = max(abs(fp), 1.0)
            for j in range(self.n):
                temp = abs(self.xi[j]) * max(abs(self.p[j]), 1.0) / den
                if temp > test:
                    test = temp

            if test < Frprmn.TINY:
                return self.p

            dgg = gg = 0.0
            for j in range(self.n):
                gg += g[j] * g[j]
                # dgg += self.xi[j] * self.xi[j]
                dgg += (self.xi[j] + g[j]) * self.xi[j]

            if gg == 0.0:
                return self.p

            gam = dgg / gg
            for j in range(self.n):
                g[j]= -self.xi[j]
                self.xi[j] = h[j] = g[j] + gam * h[j]

    def grad(self, p):
        pt = p.copy()
        deriv = zeros(self.n)
        for i in range(self.n):
            pt[i] += Frprmn.TINY
            deriv[i] = (self.func(pt) - self.func(p)) / Frprmn.TINY
            pt[i] -= Frprmn.TINY
        return deriv

def func(x):
    return (x[0] - 1.5) * (x[0] - 1.5) + x[0] * x[0] + 2.0 * x[0] * x[1] + x[1] * x[1] + (x[2]-1) * (x[2]-1)

if __name__=='__main__':
    f  = lambda x : func(x)
    p  = array([2.0, 3.0, 4.0])
    xi = array([0.0, 0.0, 0.0])
    fr = Frprmn(200, 1.0e-8, p, xi, f)
    pmin = fr.minimize()
    print(pmin)
