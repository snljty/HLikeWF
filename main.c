/* Uncomment these lines if the head file is not included in the compile options.
# ifndef __Hydrogen_Like_AO_H__
# include "AtomicWaveFunctions.h"
# endif
*/

int main()
{
  /* Apply the functions here. */
  # define Z 1
  # define a0 1.
  printf("T(0, M_PI_4) = %.6lf + %.6lf i\n", T(0, M_PI_4).real, 
         T(0, M_PI_4).imag);
  printf("TReal(1, M_PI_4) = %.6lf\n", TReal(1, M_PI_4));
  printf("TImag(1, M_PI_4) = %.6lf\n", TImag(1, M_PI_4));
  printf("S(2, 1, M_PI_4) = %.6lf\n", S(2, 1, M_PI_4));
  printf("R(3, 2, 1.5) = %.6lf\n", R(3, 2, 1.5));
  printf("YImag(3, 1, M_PI_4, M_PI / 6) = %.6lf\n", 
         YImag(3, 1, M_PI_4, M_PI / 6));
  printf("psiReal(3, 2, -2, 1.5, M_PI / 5, M_PI / 7) = %.6lf\n",
         psiReal(3, 2, -2, 1.5, M_PI / 5, M_PI / 7));
  puts("\nShould be:");
  puts("T(0, M_PI_4) = 0.398942 + 0.000000 i");
  puts("TReal(1, M_PI_4) = 0.282095");
  puts("TImag(1, M_PI_4) = 0.282095");
  puts("S(2, 1, M_PI_4) = 0.968246");
  puts("R(3, 2, 1.5) = 0.012304");
  puts("YImag(3, 1, M_PI_4, M_PI / 6) = 0.399915");
  puts("psiReal(3, 2, -2, 1.5, M_PI / 5, M_PI / 7) = 0.001024");

  system("1>NUL PAUSE");

  return 0;
}

