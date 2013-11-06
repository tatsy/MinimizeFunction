# -*- coding: utf-8 -*-
from math import *
from numpy import *
from Linemethod import Linemethod

# Powell法で関数を最小化するクラス
class Powell(Linemethod):
    TINY = 1.0e-8

    # コンストラクタ
    def __init__(self, maxit, eps, ftol, p, xi, func):
        Linemethod.__init__(self, maxit, eps, p, xi, func)
        self.ftol  = ftol
        self.func  = func
        self.ximat = zeros((self.n, self.n))
        for i in range(self.n):
            self.ximat[i,i]= 1.0

    # 関数の最小化
    def minimize(self):
        pt   = zeros(self.n)
        ptt  = zeros(self.n)
        fret = self.func(self.p)
        for j in range(self.n):
            pt[j] = self.p[j]

        for it in range(self.MAXIT):
            # どの方向に進むのが一番関数値が減るかを求める
            fp   = fret
            ibig = 0
            dl   = 0.0
            for i in range(self.n):
                for j in range(self.n):
                    self.xi[j] = self.ximat[j,i]
                fptt = fret
                try:
                    fret = self.linmin()
                except Exception, msg:
                    print(msg)
                # print('fptt = %f, fret = %f' % (fptt, fret))
                if fptt - fret > dl:
                    dl = fptt - fret
                    ibig = i

            # 終了条件判定 (勾配が一定以下になったら終了)
            if 2.0 * (fp - fret) <= self.ftol * (abs(fp) + abs(fret) + Powell.TINY):
                return self.p

            # 次の繰り返し計算の準備
            for j in range(self.n):
                ptt[j] = 2.0 * self.p[j] - pt[j]
                self.xi[j] = self.p[j] - pt[j]
                pt[j] = self.p[j]

            fptt = self.func(ptt)
            if fptt < fp:
                t = 2.0 * (fp - 2.0 * fret + fptt) * (fp - fret - dl) ** 2 - dl * (fp - fptt) ** 2
                if t < 0.0:
                    try:
                        fret = self.linmin()
                    except Exception, msg:
                        print(msg)
                    for j in range(self.n):
                        self.ximat[j, ibig]   = self.ximat[j, self.n-1]
                        self.ximat[j, self.n-1] = self.xi[j]

def func(x):
    return (x[0] - 1.5) * (x[0] - 1.5) + x[0] * x[0] + 2.0 * x[0] * x[1] + x[1] * x[1] + (x[2]-1) * (x[2]-1)

if __name__=='__main__':
    f  = lambda x : func(x)
    p  = zeros(3)
    xi = array([0.0, 0.0, 0.0])
    pw = Powell(200, 1.0e-8, 1.0e-8, p, xi, f)
    pmin = pw.minimize()
    print(pmin)

