@echo off
cd /d "%~dp0"
call abaqus verify -user_exp
call abaqus job=geometric_hybrid input=geometric_hybrid.inp user=vumat_grind.for double=both cpus=1 datacheck
if errorlevel 1 exit /b 1
call abaqus job=geometric_hybrid input=geometric_hybrid.inp user=vumat_grind.for double=both cpus=8 interactive
