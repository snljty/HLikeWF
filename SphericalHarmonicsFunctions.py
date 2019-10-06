#! /usr/bin/env python3
# -*- Coding: UTF-8 -*-

r"""
Generate the spherical harmonics functions

$$
Y_{l,m}\left(\theta,\phi\right)=S_{l,m}\left(\theta\right)T_m\left(\phi\right)
$$
"""

from __future__ import division
import sympy as sym

w = sym.Symbol('w')
theta = sym.Symbol('theta')
phi = sym.Symbol('phi')

def P(l: int, m: int, w = sym.Symbol('w')):
    r'''
Generate the associated Legendre functions.

$$
P_l^{\left\vert m\right\vert}\left(w\right)\equiv
\frac{1}{2^ll!}\left(1-w^2\right)^{\frac{\left\vert m\right\vert}2}
\frac{\mathrm d^{l+\left\vert m\right\vert}}
{\mathrm dw^{l+\left\vert m\right\vert}}
\left(w^2-1\right)^l,l=0,1,2,\cdots
$$
'''
    ans = (1 / (2 ** l * sym.factorial(l)))
    ans *= (1 - w ** 2) ** (sym.functions.elementary.complexes.Abs(m) / 2)
    ans *= sym.diff(((w ** 2 - 1) ** l), w, (l + sym.functions.elementary.complexes.Abs(m)))
    return ans

def S(l: int, m: int, theta = sym.Symbol('theta')):
    r'''
$$
S_{l,m}\left(\theta\right)=
\left(
\frac{2l+1}{2}\frac{\left(l-\left\vert m\right\vert\right)!}
{\left(l+\left\vert m\right\vert\right)!}
\right)^{\frac 1 2}
P^{\left\vert m\right\vert}_l
\left(\mathrm{cos}\left(\theta\right)\right)
$$
'''

    ans = (2 * l + 1) / 2
    ans *= sym.factorial(l - sym.functions.elementary.complexes.Abs(m))
    ans /= sym.factorial(l + sym.functions.elementary.complexes.Abs(m))
    ans **= 1 / 2
    ans *= P(l, m, sym.cos(theta))
    return ans

def T(m: int, phi = sym.Symbol('phi')):
    r'''
$$
T_m\left(\phi\right)=\frac 1{\sqrt{2\pi}}\mathrm e^{\mathrm im\phi}
$$
'''
    ans = 1 / sym.sqrt(2 * sym.pi)
    ans *= sym.E ** (sym.I * m * phi)
    return ans

def Y(l: int, m: int, theta = sym.Symbol('theta'), phi = sym.Symbol('phi')):
    r'''
$$
Y_{l,m}\left(\theta,\phi\right)=S_{l,m}\left(\theta\right)T_m\left(\phi\right)
$$
'''
    return S(l, m, theta) * T(m, phi)

