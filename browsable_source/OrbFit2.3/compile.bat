echo off
echo. This batch program compiles all OrbFit sources 
echo. using the Digital Visual FORTRAN compiler. You must 
echo. have the compiler in order to execute this batch file.
echo. See the file README.windows
nmake /nologo winstall
echo.
echo. Compilation complete!
echo.
pause
