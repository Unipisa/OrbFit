echo off
rem This batch removes all object code, test lfiles anf the JPL Ephemerides.
rem It is intended solely to prepare the directory for distribution.
nmake /nologo windist
pause
