# -*- coding:utf-8 -*-
from numpy import *
from Brent import Brent

# 特定方向で関数を最小化するクラス
class Linemethod:
    TINY = 1.0e-8

    # コンストラクタ
    def __init__(self, maxit, eps, p, xi, func):
        self.MAXIT = maxit
        self.EPS   = eps
        self.p     = p
        self.xi    = xi
        self.n     = len(p)
        self.func  = func

    # ある次元方向へ進んだ時の関数の最小値を見つける
    def linmin(self):
        a = 0.0
        b = 1.0
        f1dim  = F1dim(self.p, self.xi, self.func)
        f1d    = lambda x : f1dim.eval(x)
        br = Brent(self.MAXIT, self.EPS, f1d)
        br.bracket(a, b)
        br.minimize()
        for j in range(self.n):
            self.xi[j] *= br.xmin
            self.p[j]  += self.xi[j]

        self.xmin = self.xi
        self.fmin = br.fmin
        return self.fmin

# 実数値ベクトル関数を方向xiに限定した関数を作るクラス
class F1dim:
    def __init__(self, p, xi, func):
        self.p    = p
        self.xi   = xi
        self.func = func
        self.n    = len(p)
        self.xt   = [0] * self.n

    def eval(self, x):
        for j in range(self.n):
            self.xt[j] = self.p[j] + x * self.xi[j]
        return self.func(self.xt)

# 最小化する関数
def func(x):
    return (x[0] - 1.5) * (x[0] - 1.5) + x[0] * x[0] + 2.0 * x[0] * x[1] + x[1] * x[1] + (x[2]-1) * (x[2]-1) - 10.0

if __name__=='__main__':
    f  = lambda x : func(x)
    p  = zeros(3)
    xi = array([1.0, 0.0, 0.0])
    lm = Linemethod(100, 1.0e-8, p, xi, f)
    lm.linmin()
    print(lm.xmin)
    print(lm.fmin)
