# -*- coding: utf-8 -*-
from math import *
from numpy import *
from Brent import Brent

class Dlinemethod:
    def __init__(self, maxit, eps, p, xi, func):
        self.MAXIT = maxit
        self.EPS   = eps
        self.p     = p
        self.xi    = xi
        self.n     = len(p)
        self.func  = func

    def linmin(self):
        df1dim = Df1dim(self.p, self.xi, self.func)
        a = 0.0
        b = 1.0
        br = Brent(self.MAXIT, self.EPS, self.func)
        br.bracket(a, b)
        xmin = br.minimize()
        for j in range(self.n):
            self.xi[j] *= xmin
            self.p[j]  += self.xi[j]
        return br.fmin

class Df1dim:
    TINY = 1.0e-8
    def __init__(self, p, xi, func):
        self.p    = p
        self.xi   = xi
        self.func = func

    def eval(self, x):
        xt = zeros(self.n)
        for j in range(self.n):
            xt[j] = self.p[j] + x * self.xi[j]
        return self.func(xt)

    def df(self, x):
        df1 = 0.0
        dft = self.grad(self.xt)

    def grad(self, p):
        pt = p.copy()
        deriv = zeros(self.n)
        for i in range(self.n):
            pt[i] += Df1dim.TINY
            deriv[i] = (self.func(pt) - self.func(p)) / Df1dim.TINY
            pt[i] -= Df1dim.TINY
        return deriv
