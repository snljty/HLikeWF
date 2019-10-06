#! /usr/bin/env python3
# -*- Coding: UTF-8 -*-

r"""
Generate the associated Legendre functions.

$$
P_l^{\left\vert m\right\vert}\left(w\right)\equiv
\frac{1}{2^ll!}\left(1-w^2\right)^{\frac{\left\vert m\right\vert}2}
\frac{\mathrm d^{l+\left\vert m\right\vert}}
{\mathrm dw^{l+\left\vert m\right\vert}}
\left(w^2-1\right)^l,l=0,1,2,\cdots
$$
"""

from __future__ import division
import sympy as sym

def P(l: int, m: int, w = sym.Symbol('w')):
    ans = (1 / (2 ** l * sym.factorial(l)))
    ans *= (1 - w ** 2) ** (sym.functions.elementary.complexes.Abs(m) / 2)
    ans *= sym.diff(((w ** 2 - 1) ** l), w, (l + sym.functions.elementary.complexes.Abs(m)))
    return ans

test = sym.Symbol('w')

