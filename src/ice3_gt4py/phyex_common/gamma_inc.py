# -*- coding: utf-8 -*-
from math import gamma
from ifs_physics_common.utils.f2py import ported_method
import logging

from numpy import exp, log, finfo


@ported_method(from_file="PHYEX/src/common/aux/gamma_inc.F90")
def gamma_inc(a: float, x: float):
    """Compute the genrnalized gamma function

     The purpose of this function is to compute the generalized
    incomplete Gamma function of its argument.

                              /X
                        1     |
     GAMMA_INC(A,X)= -------- | Z**(A-1) EXP(-Z) dZ
                     GAMMA(A) |
                              /0

     Reference : Press, Teukolsky, Vetterling and Flannery: Numerical Recipes, 209-213

     Args:
         a (float): _description_
         x (float): _description_
    """

    zeps = 3e-7
    itmax = 100
    zfpmin = 1e-30

    if x < 0 or a <= 0:
        raise ValueError(f"Invalid arguments: x < 0 or a <= 0")

    if x < a + 1:
        ap = a
        zsum = 1 / a
        zdel = zsum
        jn = 1

        # LOOP_SERIES: DO
        while not (abs(zdel) < abs(zsum) * zeps):
            ap += 1
            zdel *= x / ap
            zsum += zdel
            jn += 1

            if jn > itmax:
                raise ValueError(
                    "a argument is too large or ITMAX is too small the incomplete GAMMA_INC "
                    "function cannot be evaluated correctly by the series method"
                )

        return zsum * exp(-x + a * log(x) - log(gamma(a)))

    else:
        b = x + 1 - a
        c = 1 / finfo.tiny
        d = 1 / b
        h = d
        jn = 1

        while not (abs(zdel - 1) < zeps):

            an = -jn * (jn - a)
            b += 2
            d = an * d + b

            if abs(d) < finfo.tiny:
                d = zfpmin
            if abs(c) < finfo.tiny:
                c = zfpmin

            d = 1 / d
            zdel = d * c
            h = h * zdel
            jn += 1

            if jn > itmax:
                raise ValueError(
                    "a argument is too large or ITMAX is too small the incomplete GAMMA_INC "
                    "function cannot be evaluated correctly by the series method"
                )

        return 1 - h * exp(-x + a * log(x) - log(gamma(a)))