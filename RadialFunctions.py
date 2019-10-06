#! /usr/bin/env python3
# -*- Coding: UTF-8 -*-

r"""
b_0 is calculated through the normalized integral.
$$
\int\limits_0^{+\infty}
{\left\vert R_{nl}\left(r\right)\right\vert}^2r^2\mathrm dr=1
$$

$$
b_{0\left(n,l\right)}=\frac {\left(\frac{2Z}{na_0}\right)^{\frac {2l+3} 2}}
{\sqrt{\sum\limits_{j=0}^{n-l-1}
\sum\limits_{k=0}^{n-l-1}\left(
\left(\prod\limits_{u=0}^{j-1}\frac{u+l+1-n}{\left(u+1\right)
\left(u+2l+2\right)}\right)
\left(\prod\limits_{v=0}^{k-1}\frac{v+l+1-n}{\left(v+1\right)
\left(v+2l+2\right)}\right)
\left(\left(j+k+2l+2\right)!\right)\right)}}
$$

$$
b_{j+1}=\frac{2Z}{na_0}\frac{j+l+1-n}{\left(j+1\right)\left(j+2l+2\right)}b_j
$$

$$
R_{nl}\left(r\right)=
r^l\mathrm e^{-\frac{Zr}{na_0}}\sum\limits_{j=0}^{n-l-1}{b_{j\left(n,l\right)}r^j}
$$
"""

from __future__ import division
import sympy as sym

# Z = sym.Symbol('Z')
# a0 = sym.Symbol('a0')
Z = 1
a0 = 1.
r = sym.Symbol('r')

def Calc_b_0(n: int, l: int):
    '''
Calculates b0
'''
    global Z, a0
    tmp = [0., 0., 0., 0.]
    tmp[0] = ((2 * Z) / (n * a0)) ** ((2 * l + 3) / 2)
    tmp[1] = 0.
    for j in range(n - l):
        tmp[2] = 0.
        for k in range(n - l):
            tmp[3] = sym.factorial(j + k + 2 * l + 2)
            for u in range(j):
                tmp[3] *= (u + l + 1 - n) / ((u + 1) * (u + 2 * l + 2))
            for v in range(k):
                tmp[3] *= (v + l + 1 - n) / ((v + 1) * (v + 2 * l + 2))
            tmp[2] += tmp[3]
        tmp[1] += tmp[2]
    ans = tmp[0] / sym.sqrt(tmp[1])
    return ans
               
def R(n: int, l: int, r = sym.Symbol('r')):
    r'''
Using b_0 and the recursion relationship to compute radius functions
'''
    global Z, a0
    # Calculate b_{0(n,l)} here.
    b = [Calc_b_0(n, l),]
    for index in range(n - l - 1):
        tmp = (2 * Z) / (n * a0)
        tmp *= (index + l + 1 - n)/((index + 1) * (index + 2 * l + 2))
        tmp *= b[index]
        b.append(tmp)
    ans = 0.
    for index in range(n - l): ans += b[index] * r ** index
    ans *= r ** l * sym.E ** (- (Z * r) / (n * a0))
    return ans


        
