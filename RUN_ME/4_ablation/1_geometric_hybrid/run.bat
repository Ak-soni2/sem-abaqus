@echo off
cd /d "%~dp0"
abaqus verify -user_exp
abaqus job=geometric_hybrid input=geometric_hybrid.inp user=vumat_grind.for double=both cpus=1 datacheck
if errorlevel 1 exit /b 1
abaqus job=geometric_hybrid input=geometric_hybrid.inp user=vumat_grind.for double=both cpus=8 interactive
