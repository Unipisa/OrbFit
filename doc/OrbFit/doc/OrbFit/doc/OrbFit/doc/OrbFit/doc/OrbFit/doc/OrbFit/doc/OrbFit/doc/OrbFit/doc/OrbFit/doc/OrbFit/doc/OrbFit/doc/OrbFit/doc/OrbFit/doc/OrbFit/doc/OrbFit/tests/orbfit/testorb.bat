echo off
echo.
echo. ==== Orbfit.exe: Running Automated Test ==== 

call cleanorb

echo 1998PB1 | ..\..\bin\orbfit.exe
echo 1996RO13_gauss | ..\..\bin\orbfit.exe
echo 1996RO13_diffcor | ..\..\bin\orbfit.exe
echo 1996RO13_ephem1 | ..\..\bin\orbfit.exe
echo 1996RO13_ephem2 | ..\..\bin\orbfit.exe
echo 1996RO13_ephem3 | ..\..\bin\orbfit.exe
del "1996RO13=1220T-2.rwo"
echo 1996RO13_ident1 | ..\..\bin\orbfit.exe
ren "1996RO13=1220T-2.rwo" "1996RO13=1220T-2.ident1.rwo"
copy "1996RO13=1220T-2.rwo.sav" "1996RO13=1220T-2.rwo"
copy "1996RO13=1220T-2.rwo" "1996RO13=1220T-2.in.rwo"
echo 1996RO13_ident2 | ..\..\bin\orbfit.exe
ren "1996RO13=1220T-2.rwo" "1996RO13=1220T-2.ident2.rwo"

echo. ==== Orbfit.exe: Automated Test Completed ==== 
pause
