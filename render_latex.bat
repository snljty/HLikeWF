@echo off

setlocal enableextensions
for /r . %%i in (*.ltx) do (
latex %%~nxi
del %%~ni.aux
del %%~ni.log
dvipng --png -T tight -D 1024 -bg Transparent -o %%~ni.png %%~ni.dvi
del %%~ni.dvi
)

