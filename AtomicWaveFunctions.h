# ifndef __Hydrogen_Like_AO_H__
# define __Hydrogen_Like_AO_H__

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

/*
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
*/

/* define double precision complex number structure */
typedef struct _dcomplex {double real; double imag;} dcomplex;

/* Z is the charge of the nucleus, for hydrogen it is 1 */
# define Z 1
/* a0 is the Bohr radius for hydrogen atom. In atomic unit you can set it to 1. */
# define a0 1.

/* Factorial number */
int Factorial(int n);

/* Permutation number */
int Perm(int n, int m);

/* Combination number */
int Comb(int n, int m);

/*
Generate the associated Legendre functions.
$$
P_l^{\left\vert m\right\vert}\left(w\right)\equiv
\frac{1}{2^ll!}\left(1-w^2\right)^{\frac{\left\vert m\right\vert}2}
\frac{\mathrm d^{l+\left\vert m\right\vert}}
{\mathrm dw^{l+\left\vert m\right\vert}}
\left(w^2-1\right)^l,l=0,1,2,\cdots
$$
$$
\frac{\mathrm d^{l+\left\vert m\right\vert}}
{\mathrm dw^{l+\left\vert m\right\vert}}
=\sum\limits_{j=\frac{l+\left\vert m\right\vert
+\mathrm{Mod}\left(l+\left\vert m\right\vert,2\right)}2}
^l\mathrm C_l^jP_{2j}^{l+\left\vert m\right\vert}
w^{2j-\left(l+\left\vert m\right\vert\right)}
$$
*/
double P(int l, int m, double w);

/*
$$
S_{l,m}\left(\theta\right)=
\left(
\frac{2l+1}{2}\frac{\left(l-\left\vert m\right\vert\right)!}
{\left(l+\left\vert m\right\vert\right)!}
\right)^{\frac 1 2}
P^{\left\vert m\right\vert}_l
\left(\mathrm{cos}\left(\theta\right)\right)
$$
*/
double S(int l, int m, double theta);

/*
$$
\mathrm{Real}\left(T_m\left(\phi\right)\right)=
\frac 1{\sqrt{2\pi}}\mathrm{cos}\left(m\phi\right)
$$
*/
double TReal(int m, double phi);

/*
$$
\mathrm{Imag}\left(T_m\left(\phi\right)\right)=
\frac 1{\sqrt{2\pi}}\mathrm{sin}\left(m\phi\right)
$$
*/
double TImag(int m, double phi);

/*
$$
T_m\left(\phi\right)=\frac 1{\sqrt{2\pi}}\mathrm e^{\mathrm im\phi}
$$
*/
dcomplex T(int m, double phi);

/*
$$
\mathrm{Real}\left(Y_{l,m}\left(\theta,\phi\right)\right)=
S_{l,m}\left(\theta\right)\mathrm{Real}\left(T_m\left(\phi\right)\right)
$$
*/
double YReal(int l, int m, double theta, double phi);

/*
$$
\mathrm{Imag}\left(Y_{l,m}\left(\theta,\phi\right)\right)=
S_{l,m}\left(\theta\right)\mathrm{Imag}\left(T_m\left(\phi\right)\right)
$$
*/
double YImag(int l, int m, double theta, double phi);

/*
$$
Y_{l,m}\left(\theta,\phi\right)=S_{l,m}\left(\theta\right)T_m\left(\phi\right)
$$
*/
dcomplex Y(int l, int m, double theta, double phi);

/*
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
*/
double Calc_b_0(int n, int l);

/*
$$
R_{nl}\left(r\right)=
r^l\mathrm e^{-\frac{Zr}{na_0}}\sum\limits_{j=0}^{n-l-1}{b_{j\left(n,l\right)}r^j}
$$
*/
double R(int n, int l, double r);

/*
$$
\mathrm{Real}\left(\psi\left(r,\theta,\phi\right)\right)
=R_{n,l}\left(r\right)\mathrm{Real}
\left(Y_{l,m}\left(\theta,\phi\right)\right)
$$
*/
double psiReal(int n,int l, int m, double r, double theta, double phi);

/*
$$
\mathrm{Imag}\left(\psi\left(r,\theta,\phi\right)\right)
=R_{n,l}\left(r\right)\mathrm{Imag}
\left(Y_{l,m}\left(\theta,\phi\right)\right)
$$
*/
double psiImag(int n,int l, int m, double r, double theta, double phi);

/*
$$
\psi\left(r,\theta,\phi\right)
=R_{n,l}\left(r\right)Y_{l,m}\left(\theta,\phi\right)
$$
*/
dcomplex psi(int n,int l, int m, double r, double theta, double phi);

# endif

