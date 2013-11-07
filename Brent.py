#-*- coding: utf-8 -*-
from math import *
from numpy import *
from Bracket import Bracketmethod

class Brent(Bracketmethod):
    def __init__(self, maxit, tol, func):
        self.maxit = maxit
        self.tol   = tol
        self.func  = func

    def minimize(self):
        CGOLD = 0.3819660
        ZEPS  = 1.0e-12
        d = e = 0.0

        a = self.ax if self.ax < self.cx else self.cx
        b = self.ax if self.ax > self.cx else self.cx
        x = w= v = self.bx
        fw = fv = fx = self.func(x)

        for it in range(self.maxit):
            xm = 0.5 * (a + b)
            tol1 = self.tol * abs(x) + ZEPS
            tol2 = 2.0 * tol1
            if(abs(x - xm) <= (tol2 - 0.5 * (b - a))):
                self.fmin = fx
                self.xmin = x
                return self.xmin

            if abs(e) > tol1:
                r = (x - w) * (fx - fv)
                q = (x - v) * (fx - fw)
                p = (x - v) * q - (x - w) * r
                q = 2.0 * (q - r)
                p = -p if q < 0.0 else p
                q = abs(q)
                etemp = e
                e = d
                if abs(p) >= abs(0.5 * q * etemp) or p <= q * (a - x) or p >= q * (b - x):
                    e = a - x if x >= xm else b - x
                    d = CGOLD * e
                else:
                    d = p / q
                    u = x + d
                    if u - a < tol2 or b - u < tol2:
                        d = tol1 * self.sign(xm - x)
            else:
                e = a - x if x >= xm else b - x
                d = CGOLD * e

            u = x + d if abs(d) >= tol1 else x + tol1 * self.sign(d)
            fu = self.func(u)

            if fu <= fx:
                if u >= x:
                    a = x
                else:
                    b = x

                (v, w, x) = (w, x, u)
                (fv, fw, fx) = (fw, fx, fu)
            else:
                if u < x:
                    a = u
                else:
                    b = u
                if fu <= fw or w == x:
                    (v, w, fv, fw) = (w, u, fw, fu)
                elif fu <= fv or v == x or v == w:
                    (v, fv) = (u, fu)

        print('not converge')
        self.xmin = 0.0
        self.fmin = 0.0

if __name__ == '__main__':
    f = lambda x : (x - 1.0) * (x - 11.0) * (x - 2) * (x - 4)
    brent = Brent(100, 3.0e-8, f)
    brent.bracket(0.0, 1.0)
    print('%f, %f' % (brent.ax, brent.bx))
    brent.minimize()
    print('xmin = %f, fmin = %f' % (brent.xmin, brent.fmin))
