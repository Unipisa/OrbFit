echo off
echo.
echo. ==== FITOBS.exe: Running Automated Test ==== 
call cleanfit
..\..\bin\fitobs.exe < 1220T-2.inp
echo. ==== FITOBS.exe: Automated Test Completed ==== 
pause
