#! /usr/bin/env python3
# -*- Coding: UTF-8 -*-

r"""
Time-independant normalized wave function for hydrogen-like atoms.

The comments below quoted with "$$ $$" are LaTeX codes, 
which are mathematical formulas used.
If you want to see the graphic, copy it to a LaTeX template below:

\documentclass[11pt,a4paper]{article}
\pagestyle{empty}
\begin{document}
\begin{center}
$$
%% Paste the LaTeX expression here.
$$
\end{center}
\end{document}

and save it as a xxx.ltx file.
Then compile it with command:
latex xxx.ltx
Render it with:
dvipng --png -T tight -D 1024 -bg Transparent -o xxx.png xxx.dvi
Or if you want a white-background version, render with:
dvipng --png -T tight -D 1024 -o xxx.png xxx.dvi
If you like a pdf file, just compile with:
pdflatex xxx.ltx
"""

from __future__ import division
import sympy as sym

Z = sym.Symbol('Z')
a0 = sym.Symbol('a0')
# Z is the charge of the nucleus, for hydrogen it is 1
# Z = 1
# a0 is the Bohr radius for hydrogen atom. In atomic unit you can set it to 1.
# a0 = 1.
r = sym.Symbol('r')
w = sym.Symbol('w')
theta = sym.Symbol('theta')
phi = sym.Symbol('phi')

def AssoLegendre(l: int, m: int, w = sym.Symbol('w')):
    r'''
Generate the associated Legendre functions.

$$
AssociatedLegendre_l^{\left\vert m\right\vert}\left(w\right)\equiv
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

def P(l: int, m: int, w = sym.Symbol('w')):
    r'''
Same as associated Legendre functions, but expanded the derivative
expression.
sym.ff is the permutation number.
sym.binomial is the combination number.

$$
\frac{\mathrm d^{l+\left\vert m\right\vert}}
{\mathrm dw^{l+\left\vert m\right\vert}}
=\sum\limits_{j=\frac{l+\left\vert m\right\vert
+\mathrm{Mod}\left(l+\left\vert m\right\vert,2\right)}2}
^l\mathrm C_l^jP_{2j}^{l+\left\vert m\right\vert}
w^{2j-\left(l+\left\vert m\right\vert\right)}
$$
'''
    ans = 0.
    tmp = l + sym.functions.elementary.complexes.Abs(m)
    for j in range((tmp + tmp % 2) // 2, l + 1):
        ans += (sym.binomial(l, j) * sym.ff(2 * j, tmp) * 
                w ** (2 * j - tmp))
    ans *= (1 - w ** 2) ** (sym.functions.elementary.complexes.Abs(m) / 2)
    ans /= 2 ** l * sym.factorial(l)
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

def TReal(m: int, phi = sym.Symbol('phi')):
    r'''
$$
\mathrm{Real}\left(T_m\left(\phi\right)\right)=
\frac 1{\sqrt{2\pi}}\mathrm{cos}\left(m\phi\right)
$$
'''
    ans = 1 / sym.sqrt(2 * sym.pi)
    ans *= sym.cos(m * phi)
    return ans

def TImag(m: int, phi = sym.Symbol('phi')):
    r'''
$$
\mathrm{Imag}\left(T_m\left(\phi\right)\right)=
\frac 1{\sqrt{2\pi}}\mathrm{sin}\left(m\phi\right)
$$
'''
    ans = 1 / sym.sqrt(2 * sym.pi)
    ans *= sym.sin(m * phi)
    return ans

def Y(l: int, m: int, theta = sym.Symbol('theta'), phi = sym.Symbol('phi')):
    r'''
$$
Y_{l,m}\left(\theta,\phi\right)=S_{l,m}\left(\theta\right)T_m\left(\phi\right)
$$
'''
    return S(l, m, theta) * T(m, phi)

def YReal(l: int, m: int, theta = sym.Symbol('theta'), phi = sym.Symbol('phi')):
    r'''
$$
\mathrm{Real}\left(Y_{l,m}\left(\theta,\phi\right)\right)=
S_{l,m}\left(\theta\right)\mathrm{Real}\left(T_m\left(\phi\right)\right)
$$
'''
    return S(l, m, theta) * TReal(m, phi)

def YImag(l: int, m: int, theta = sym.Symbol('theta'), phi = sym.Symbol('phi')):
    r'''
$$
\mathrm{Imag}\left(Y_{l,m}\left(\theta,\phi\right)\right)=
S_{l,m}\left(\theta\right)\mathrm{Imag}\left(T_m\left(\phi\right)\right)
$$
'''
    return S(l, m, theta) * TImag(m, phi)

def Calc_b_0(n: int, l: int):
    r'''
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
$$
R_{nl}\left(r\right)=
r^l\mathrm e^{-\frac{Zr}{na_0}}\sum\limits_{j=0}^{n-l-1}{b_{j\left(n,l\right)}r^j}
$$
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

def psi(n: int, l: int, m: int, r = sym.Symbol('r'),
        theta = sym.Symbol('theta'), phi = sym.Symbol('phi')):
    r'''
$$
\psi\left(r,\theta,\phi\right)
=R_{n,l}\left(r\right)Y_{l,m}\left(\theta,\phi\right)
$$
'''
    return R(n, l, r) * Y(l, m, theta, phi)

def psiReal(n: int, l: int, m: int, r = sym.Symbol('r'),
            theta = sym.Symbol('theta'), phi = sym.Symbol('phi')):
    r'''
$$
\mathrm{Real}\left(\psi\left(r,\theta,\phi\right)\right)
=R_{n,l}\left(r\right)\mathrm{Real}
\left(Y_{l,m}\left(\theta,\phi\right)\right)
$$
'''
    return R(n, l, r) * YReal(l, m, theta, phi)

def psiImag(n: int, l: int, m: int, r = sym.Symbol('r'),
            theta = sym.Symbol('theta'), phi = sym.Symbol('phi')):
    r'''
$$
\mathrm{Imag}\left(\psi\left(r,\theta,\phi\right)\right)
=R_{n,l}\left(r\right)\mathrm{Imag}
\left(Y_{l,m}\left(\theta,\phi\right)\right)
$$
'''
    return R(n, l, r) * YImag(l, m, theta, phi)

