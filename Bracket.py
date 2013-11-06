from math import *
from numpy import *

def bracket(a, b, func):
    fa = func(a)
    fb = func(b)
    if fa > 0.0 and fb > 0.0 and fa < fb:
        (a, b) = (b, a)

    if fa < 0.0 and fb < 0.0 and fa > fb:
        (a, b) = (b, a)

    it = 0
    MAXIT = 200
    while True:
        fa = func(a)
        fb = func(b)
        # print('[%03d] %f, %f' % (it+1, fa, fb))
        if fa * fb <= 0.0:
            break
        b = b + (b - a)

        it += 1
        if it > MAXIT:
            raise Exception, '[Bracket] Too much iteration'

    return (a, b)
