@echo off

rem compile method:

gcc -o AtomicWaveFunctions.o -c AtomicWaveFunctions.c

gcc -shared -fPIC -o libAtomicWaveFunctions.dll AtomicWaveFunctions.o
ar rcs libAtomicWaveFunctions.lib AtomicWaveFunctions.o

gcc -o main.exe main.c -L . -l AtomicWaveFunctions ^
-I . -include "AtomicWaveFunctions.h"

