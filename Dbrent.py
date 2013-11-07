# -*- coding: utf-8 -*-
from numpy import *
from Decker import Decker

class Brent(Decker):
    # コンストラクタ
    def __init__(self, maxit, eps, func):
        Decker.__init__(self, maxit, eps, func)

    # 関数の最小化 (bracketingされている前提)
    def minimize(self):
        if self.df(self.a) * self.df(self.b) > Brent.TINY:
            raise Exception, 'Please set initial two points with different sign'

        c = self.a
        d = 0.0
        mflag = True
        for it in range(self.MAXIT):
            fa = self.df(self.a)
            fb = self.df(self.b)
            fc = self.df(c)
            if abs(fa - fc) > Brent.TINY and abs(fb - fc) > Brent.TINY:
                s = self.inv_quad(self.a, self.b, c, self.df)
            else:
                s = self.b - (self.b - self.a) / (fb - fa + Brent.TINY) * fb

            cond1 = (s < (3.0 * self.a + self.b) / 4.0 or self.b < s)
            cond2 = (mflag and abs(s - self.b) >= abs(self.b - c) / 2.0)
            cond3 = (not mflag and abs(s - self.b) >= abs(c - d) / 2.0)
            cond4 = (mflag and abs(self.b - c) < Brent.TINY)
            cond5 = (not mflag and abs(c - d) < Brent.TINY)
            if cond1 or cond2 or cond3 or cond4 or cond5:
                s = (self.a + self.b) / 2.0
                mflag = True
            else:
                mflag = False

            fs = self.df(s)
            d = c
            c = self.b
            if fa * fs < 0.0:
                self.b = s
            else:
                self.a = s

            if abs(fa) < abs(fb):
                (self.a, self.b) = (self.b, self.a)

            if abs(self.a - self.b) < Brent.TINY:
                break

        self.fmin = self.func(self.a)
        self.xmin = self.a
        return self.xmin

    # 逆二次補間
    def inv_quad(self, a, b, c, func):
        fa = func(a)
        fb = func(b)
        fc = func(c)
        s = 0.0
        s += (a * fb * fc) / ((fa - fb) * (fa - fc))
        s += (fa * b * fc) / ((fb - fa) * (fb - fc))
        s += (fa * fb * c) / ((fc - fa) * (fc - fb))
        return s

if __name__=='__main__':
    func = lambda x : (x - 1) * (x - 4)
    br = Brent(100, 1.0e-8, func)
    br.bracket(0.0, 1.0)
    print('a = %f, b = %f' % (br.a, br.b))
    br.minimize()
    print('xmin = %f, fmin = %f' % (br.xmin, br.fmin))

