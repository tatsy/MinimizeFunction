# -*- coding: utf-8 -*-
from math import *
from numpy import *

# 主双対内点法による線形計画の解法
class Intpt:
    BIG   = 1.0e12
    DELTA = 0.02
    SIGMA = 0.9

    # コンストラクタ
    def __init__(self, A, b, c, x):
        self.A = A
        self.b = b
        self.c = c
        self.m = len(A)
        self.n = len(A[0])
        self.x = zeros(self.n)
        self.y = zeros(self.m)
        self.z = zeros(self.n)

    # 目的関数の最小化
    def minimize(self, maxit, eps):
        At = self.A.transpose().copy()

        rpfact = 1.0 + sqrt(dot(self.b, self.b))
        rdfact = 1.0 + sqrt(dot(self.c, self.c))
        for j in range(self.n):
            self.x[j] = 1000.0
            self.z[j] = 1000.0

        for i in range(self.m):
            self.y[i] = 1000.0

        normrp_old = Intpt.BIG
        normrd_old = Intpt.BIG
        rp = zeros(self.m)
        rd = zeros(self.n)
        d  = zeros(self.n)

        tempn = zeros(self.n)
        rhs   = zeros(self.m)
        dx    = zeros(self.n)
        dy    = zeros(self.m)
        dz    = zeros(self.n)

        status = 0
        for it in range(maxit):
            # 収束判定
            ax = dot(self.A, self.x)
            for i in range(self.m):
                rp[i] = ax[i] - self.b[i]
            normrp = sqrt(dot(rp, rp)) / rpfact
            aty = dot(At, self.y)
            for j in range(self.n):
                rd[j] = aty[j] + self.z[j] - self.c[j]
            normrd = sqrt(dot(rd, rd))/ rdfact
            gamma  = dot(self.x, self.z)
            mu = Intpt.DELTA * gamma / self.n
            pobj = dot(self.c, self.x)
            dobj = dot(self.b, self.y)
            print('pobj = %e, dobj = %e' % (pobj, dobj))
            gamma_norm = gamma / (1.0 + abs(pobj))

            if (normrp < eps) and (normrd < eps) and (gamma_norm < eps):
                # 収束した
                status = 0
                return status
            if (normrp > 1000.0 * normrp_old) and (normrp > eps):
                # 主問題が実行不能
                status = 1
                return status
            if (normrd > 1000.0 * normrd_old) and (normrd > eps):
                # 双対問題が実行不能
                status = 2
                return status

            # 繰り返しの計算
            for j in range(self.n):
                d[j] = self.x[j] / self.z[j]

            D = diag(d)
            adat = dot(dot(self.A, D), At)
            for j in range(self.n):
                tempn[j] = self.x[j] - mu / self.z[j] - d[j] * rd[j]
            tempm = dot(self.A, tempn)
            for i in range(self.m):
                rhs[i] = -rp[i] + tempm[i]

            dy = dot(linalg.pinv(adat), rhs)
            tempn = dot(At, dy)
            for j in range(self.n):
                dz[j] = -tempn[j] - rd[j]
            for j in range(self.n):
                dx[j] = -d[j] * dz[j] + mu / self.z[j] - self.x[j]

            # ステップ幅を計算
            alpha_p = 1.0
            for j in range(self.n):
                if self.x[j] + alpha_p * dx[j] < 0.0:
                    alpha_p = -self.x[j] / dx[j]

            alpha_d = 1.0
            for j in range(self.n):
                if self.z[j] + alpha_d * dz[j] < 0.0:
                    alpha_d = -self.z[j] / dz[j]

            alpha_p = min(alpha_p * Intpt.SIGMA, 1.0)
            alpha_d = min(alpha_d * Intpt.SIGMA, 1.0)
            for j in range(self.n):
                self.x[j] += alpha_p * dx[j]
                self.z[j] += alpha_d * dz[j]
            for i in range(self.m):
                self.y[i] += alpha_d * dy[i]
            normrp_old = normrp
            normrd_old = normrd

        # 繰り返し回数がmaxitを超えた
        status = 3
        return status

if __name__=='__main__':
    x  = zeros(2)
    c  = -array([400.0, 300.0])
    b  = -array([3800.0, 2100.0, 1200.0])
    A  = -array([[60.0, 40.0], [20.0, 30.0], [20.0, 10.0]])
    ip = Intpt(A, b, c, x)
    s  = ip.minimize(30, 1.0e-8)
    print(s)
    print(ip.x)

