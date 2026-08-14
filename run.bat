@echo off
REM Convenience launcher: pins the interpreter that has the dependencies installed.
REM The system default "python" on this machine is a broken msys64 build whose pip
REM cannot import, so calling it directly fails with ModuleNotFoundError: cv2.
REM
REM Usage:
REM   run.bat analyze "*.tif" -o results --verify
REM   run.bat wheel results -o wheel_30deg --diameter 100 --width 10 --sector 30 ^
REM           --rim-depth 2 --areal-density 40 --grain-element R3D3 --verify
REM   run.bat --help

setlocal
set "PY=%LOCALAPPDATA%\Programs\Python\Python312\python.exe"

if not exist "%PY%" (
    echo ERROR: Python 3.12 not found at "%PY%"
    echo Install it, or edit PY in this file to point at an interpreter that has
    echo the packages in requirements.txt installed.
    exit /b 1
)

"%PY%" -W ignore -m semgrit %*
