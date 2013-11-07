# -*- coding: utf-8 -*-
from math import *
from numpy import *

# 値の入れ替え
def swap(a, b):
    return (b, a)

# Amoebaクラス
class Amoeba:
    # コンストラクタ
    def __init__(self, ftol):
        self.ftol  = ftol   # 精度
        self.nfunc = 0      # 関数の評価回数
        self.mpts  = 0      # simplexの点の数
        self.ndim  = 0      # 次元
        self.fmin  = 0.0    # 結果の格納先
        self.y     = []
        self.p     = []
        self.psum  = []

    # 最小化メソッド
    def minimize(self, point, delt, func):
        # 初期点の設定
        dels = [delt] * len(point)
        self.ndim = len(point)
        pp = zeros((self.ndim + 1, self.ndim))
        for i in range(self.ndim + 1):
            for j in range(self.ndim):
                pp[i, j] = point[j]
            if i != 0:
                pp[i, i - 1] += dels[i - 1]

        # 繰り返し計算の用意
        NMAX = 2000
        TINY = 1.0e-8
        self.mpts = len(pp)
        self.ndim = len(pp[0])

        # 関数の評価
        self.p    = pp
        self.y    = zeros(self.mpts)
        self.psum = zeros(self.ndim)
        pmin = zeros(self.ndim)
        x = zeros(self.ndim)
        for i in range(self.mpts):
            for j in range(self.ndim):
                x[j] = self.p[i,j]
            self.y[i] = func(x)

        nfunc = 0
        self.get_psum()
        while True:
            # 最小値と最大値を与える要素の次元を調べる
            ilo = 0
            if self.y[0] > self.y[1]:
                ihi  = 1
                inhi = 0
            else:
                ihi  = 0
                inhi = 1

            for i in range(self.mpts):
                if self.y[i] <= self.y[ilo]:
                    ilo = i
                if self.y[i] > self.y[ihi]:
                    inhi = ihi
                    ihi = i
                elif (self.y[i] > self.y[inhi]) and (i != ihi):
                    inhi = i

            rtol = 2.0 * abs(self.y[ihi] - self.y[ilo]) / (abs(self.y[ihi]) + abs(self.y[ihi]) + TINY)
            # 精度が一定以下になったら計算を終了
            if rtol < self.ftol:
                (self.y[0], self.y[ilo]) = swap(self.y[0], self.y[ilo])
                for i in range(self.ndim):
                    (self.p[0, i], self.p[ilo, i]) = swap(self.p[0, i], self.p[ilo, i])
                    pmin[i] = self.p[0, i]
                self.fmin = self.y[0]
                return pmin

            if nfunc >= NMAX:
                raise Exception, 'Too much iteration'

            nfunc += 2

            # より評価値が下がるようにyを変更
            ytry = self.amotry(ihi, -1.0, func)
            if ytry <= self.y[ilo]:
                ytry = self.amotry(ihi, 2.0, func)
            elif ytry >= self.y[inhi]:
                ysave = self.y[ihi]
                ytry = self.amotry(ihi, 0.5, func)
                if ytry >= ysave:
                    for i in range(self.mpts):
                        if i != ilo:
                            for j in range(self.ndim):
                                self.p[i, j] = self.psum[j] = 0.5 * (self.p[i, j] + self.p[ilo, j])
                            self.y[i] = func(self.psum)

                    nfunc += self.ndim
                    self.get_psum()
            else:
                nfunc -= 1

    # 各次元の数の合計
    def get_psum(self):
        for j in range(self.ndim):
            sm = 0.0
            for i in range(self.mpts):
                sm += self.p[i, j]
            self.psum[j] = sm

    # 各次元を圧縮した時の関数の評価値を計算する
    def amotry(self, ihi, fac, func):
        ptry = zeros(self.ndim)
        fac1 = (1.0 - fac) / self.ndim
        fac2 = fac1 - fac
        for j in range(self.ndim):
            ptry[j] = self.psum[j] * fac1 - self.p[ihi][j] * fac2

        ytry = func(ptry)
        if ytry < self.y[ihi]:
            self.y[ihi] = ytry
            for j in range(self.ndim):
                self.psum[j] += ptry[j] - self.p[ihi, j]
                self.p[ihi, j] = ptry[j]
        return ytry

# 評価する関数
def f(x):
    ret = 0.0
    ndim = len(x)
    for i in range(ndim):
        d = (x[i] - (i+1))
        ret += d * d
    return ret

def g(x):
    return (x[0] - 1.5) * (x[0] - 1.5) + x[0] * x[0] + 2.0 * x[0] * x[1] + x[1] * x[1] + (x[2]-1) * (x[2]-1)

if __name__ == '__main__':
    a = Amoeba(1.0e-8)
    point = array([0] * 3)
    func = lambda x : g(x)
    r = a.minimize(point, 0.5, func)
    print(r)
