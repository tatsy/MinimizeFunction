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
        if self.dfunc(self.ax) * self.dfunc(self.bx) >= TINY:
            raise Exception, 'Please set initial two points with different sign'

        cx = self.ax
        for it in range(self.MAXIT):
            fb = self.dfunc(self.bx)
            fc = self.dfunc(cx)
            s = self.bx - (self.bx - cx) / (fb - fc + Decker.TINY) * fb
            m = (self.ax + self.bx) / 2.0
            if (self.bx < s < m) or (m < s < self.bx):
                newbx = s
            else:
                newbx = m

            if self.dfunc(self.ax) * self.dfunc(newbx) >= 0:
                self.ax = self.bx
            cx = self.bx
            self.bx = newbx

            print('[%03d] a = %.8f, b = %.8f' % (it, self.ax, self.bx))
            if abs(self.ax - self.bx) < TINY:
                break

        self.fmin = func(self.ax)
        self.xmin = self.ax

    def dfunc(self, x):
        return (self.func(x+Decker.TINY) - self.func(x)) / Decker.TINY

if __name__=='__main__':
    func = lambda x : (x - 3) * (x - 1)
    deck = Decker(100, 1.0e-20, func)
    deck.bracket(-10.0, 10.0)
    print('a = %f, b = %f' % (deck.ax, deck.bx))
    deck.minimize()
    print('xmin = %f, fmin = %f' % (deck.xmin, deck.fmin))

