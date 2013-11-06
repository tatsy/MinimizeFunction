# -*- coding: utf-8 -*-
from numpy import *
from math import *

class Qnewton:
    TINY   = 1.0e-8
    STPMAX = 100.0

    def __init__(self, maxit, eps, p, xi, func):
        self.MAXIT = maxit
        self.EPS   = eps
        self.p     = p
        self.xi    = xi
        self.func  = func

    def dfpmin(self):
        fp  = self.func(self.p)
        g   = self.grad(self.p)
        hes = zeros(self.n, self.n)
        sm  = 0.0
        for i in range(self.n):
            for j in range(self.n):
                hes[i,j] = 0.0
            hes[i,i] = 1.0
            self.xi[i] = -g[i]
            sm += self.p[i] * self.p[i]

        stpmax = Qnewton.STPMAX * max(sqrt(sm), self.n)

    def grad(self, p):
        pt = p.copy()
        deriv = zeros(self.n)
        for i in range(self.n):
            pt[i] += Qnewton.TINY
            deriv[i] = (self.func(pt) - self.func(p)) / Qnewton.TINY
            pt[i] -= Qnewton.TINY
        return deriv
