# ifndef __Hydrogen_Like_AO_H__
# include "AtomicWaveFunctions.h"
# endif

int Factorial(int n)
{
  int index;
  int ans;

  ans = 1;
  for(index = n; index > 0; index --)
    ans *= index;
  return ans;
}

int Perm(int n, int m)
{
  int index;
  int ans;

  ans = 1;
  for(index = n; index > n - m; index --)
    ans *= index;
  return ans;
}

int Comb(int n, int m)
{
  if(2 * m > n)
    return Comb(n, n - m);
  return Perm(n, m) / Factorial(m);
}

double P(int l, int m, double w)
{
  double ans;
  int tmp;
  int j;

  ans = 0.;
  tmp = l + abs(m);
  for(j = (tmp + tmp % 2) / 2; j <= l; j ++)
    ans += Comb(l, j) * Perm(2 * j, tmp) * pow(w, 2 * j - tmp);
  ans *= pow((1 - pow(w, 2)), (abs(m) / 2.));
  ans /= pow(2, l) * Factorial(l);
  return ans;
}

double S(int l, int m, double theta)
{
  double ans;

  ans = (2 * l + 1) / 2.;
  ans *= Factorial(l - abs(m));
  ans /= Factorial(l + abs(m));
  ans = pow(ans, 1 / 2.);
  ans *= P(l, m, cos(theta));
  return ans;
}

double TReal(int m, double phi)
{
  double ans;

  ans = 1 / sqrt(2 * M_PI);
  ans *= cos(m * phi);
  return ans;
}

double TImag(int m, double phi)
{
  double ans;

  ans = 1 / sqrt(2 * M_PI);
  ans *= sin(m * phi);
  return ans;
}

dcomplex T(int m, double phi)
{
  dcomplex ans;

  ans.real = TReal(m, phi);
  ans.imag = TImag(m, phi);
  return ans;
}

double YReal(int l, int m, double theta, double phi)
{
  return S(l, m, theta) * TReal(m, phi);
}

double YImag(int l, int m, double theta, double phi)
{
  return S(l, m, theta) * TImag(m, phi);
}

dcomplex Y(int l, int m, double theta, double phi)
{
  dcomplex ans;

  ans.real = YReal(l, m, theta, phi);
  ans.imag = YImag(l, m, theta, phi);
  return ans;
}

double Calc_b_0(int n, int l)
{
  double ans;
  double tmp[4];
  int j, k, u, v;

  tmp[0] = pow((double)(2 * Z) / (n * a0), (2 * l + 3) / 2.);
  tmp[1] = 0.;
  for(j = 0; j < n - l; j ++)
  {
    tmp[2] = 0.;
    for(k = 0; k < n - l; k ++)
    {
      tmp[3] = Factorial(j + k + 2 * l + 2);
      for(u = 0; u < j; u ++)
        tmp[3] *= (double)(u + l + 1 - n) / ((u + 1) * (u + 2 * l + 2));
      for(v = 0; v < k; v ++)
        tmp[3] *= (double)(v + l + 1 - n) / ((v + 1) * (v + 2 * l + 2));
      tmp[2] += tmp[3];
    }
    tmp[1] += tmp[2];
  }
  ans = tmp[0] / sqrt(tmp[1]);
  return ans;
}

double R(int n, int l, double r)
{
  double* b = NULL;
  double tmp;
  double ans;
  int index;

  b = (double*)malloc((n - l) * sizeof(double));
  /* Calculate b_{0(n,l)} here. */
  b[0] = Calc_b_0(n, l);
  for(index = 0; index < n - l - 1; index ++)
  {
    tmp = (double)(2 * Z) / (n * a0);
    tmp *= (double)(index + l + 1 - n) / ((index + 1) * (index + 2 * l + 2));
    tmp *= b[index];
    b[index + 1] = tmp;
  }
  ans = 0.;
  for(index = 0; index < n - l; index ++)
    ans += b[index] * pow(r, index);
  ans *= pow(r, l) * exp(- (double)(Z * r) / (double)(n * a0));
  free(b);
  b = NULL;
  return ans;
}

double psiReal(int n,int l, int m, double r, double theta, double phi)
{
  return R(n, l, r) * YReal(l, m, theta, phi);
}

double psiImag(int n,int l, int m, double r, double theta, double phi)
{
  return R(n, l, r) * YImag(l, m, theta, phi);
}

dcomplex psi(int n,int l, int m, double r, double theta, double phi)
{
  dcomplex ans;

  ans.real = psiReal(n, l, m, r, theta, phi);
  ans.imag = psiImag(n, l, m, r, theta, phi);
  return ans;
}

