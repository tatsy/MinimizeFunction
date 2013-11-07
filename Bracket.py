#-*- coding: utf-8 -*-
from math import *
from numpy import *

class Bracketmethod:
    # コンストラクタ
    def __init__(self, func):
        self.func = func

    def sign(self, x):
        if x > 0.0:
            return 1.0
        elif x < 0.0:
            return -1.0
        else:
            return 0.0

    # bracketingメソッド
    def bracket(self, a, b):
        # 定数
        GOLD   = (1.0 + sqrt(5.0)) / 2.0
        GLIMIT = 100.0
        TINY   = 1.0e-20

        self.ax = a
        self.bx = b
        fa = self.func(self.ax)
        fb = self.func(self.bx)
        if fb > fa:
            (self.ax, self.bx) = (self.bx, self.ax)
            (fa, fb) = (fb, fa)

        self.cx = self.bx + GOLD * (self.bx - self.ax)
        fc = self.func(self.cx)
        while fb > fc:
            print('%f %f %f' % (self.ax, self.bx, self.cx))
            r = (self.bx - self.ax) * (fb - fc)
            q = (self.bx - self.cx) * (fb - fa)
            u = self.bx - ((self.bx - self.cx) * q - (self.bx - self.ax) * r) / (2.0 * max(abs(q-r), TINY) * sign(q - r))
            ulim = self.bx + GLIMIT * (self.cx - self.bx)

            if (self.bx - u) * (u - self.cx) > 0.0:
                fu = self.func(u)
                if fu < fc:
                    self.ax = self.bx
                    self.bx = u
                    fa = fb
                    fb = fu
                    return
                elif fu > fb:
                    self.cx = u
                    fc = fu
                    return
                u = self.cx + GOLD * (self.cx - self.bx)
                fu = self.func(u)
            elif (self.cx - u) * (u - ulim) > 0.0:
                fu = self.func(u)
                if fu < fc:
                    (self.bx, self.cx, u) = (self.cx, u, u + GOLD * (u -self.cx))
                    (fb, fc, fu) = (fc, fu, self.func(u))
            elif (u - ulim) * (ulim - self.cx) >= 0.0:
                u = ulim
                fu = self.func(u)
            else:
                u = self.cx + self.cx + GOLD * (self.cx - self.bx)
                fu = self.func(u)

            (self.ax,self.bx,self.cx) = (self.bx,self.cx,u)
            (fa,fb,fc) = (fb,fc,fu)

if __name__=='__main__':
    f  = lambda x : (x - 2.0) * (x - 4.0)
    br = Bracketmethod(f)
    br.bracket(0.0, 1.0)
    print('%f %f' % (br.ax, br.bx))

