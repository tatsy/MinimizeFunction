# -*- coding: utf-8 -*-
from math  import *
from numpy import *
from Goldsec import Goldsec

# Decker法を実行するクラス
class Decker(Goldsec):

    # コンストラクタ
    def __init__(self, maxit, eps, func):
        Goldsec.__init__(self, maxit, eps, func)

    # 関数の最小化 (bracketingされている前提)
    def minimize(self):
        TINY = 1.0e-8
        if self.dfunc(self.a) * self.dfunc(self.b) >= TINY:
            raise Exception, 'Please set initial two points with different sign'

        c = self.a
        for it in range(self.MAXIT):
            fb = self.dfunc(self.b)
            fc = self.dfunc(c)
            s = self.b - (self.b - c) / (fb - fc + Decker.TINY) * fb
            m = (self.a + self.b) / 2.0
            if (self.b < s < m) or (m < s < self.b):
                newb = s
            else:
                newb = m

            if self.dfunc(self.a) * self.dfunc(newb) >= 0:
                self.a = self.b
            c = self.b
            self.b = newb

            print('[%03d] a = %.8f, b = %.8f' % (it, self.a, self.b))
            if abs(self.a - self.b) < TINY:
                break

        self.fmin = func(self.a)
        self.xmin = self.a

    def dfunc(self, x):
        return (self.func(x+Decker.TINY) - self.func(x)) / Decker.TINY

if __name__=='__main__':
    func = lambda x : (x - 3) * (x - 1)
    deck = Decker(100, 1.0e-12, func)
    deck.bracket(-10.0, 10.0)
    print('a = %f, b = %f' % (deck.a, deck.b))
    deck.minimize()
    print('xmin = %f, fmin = %f' % (deck.xmin, deck.fmin))

