@echo off
cd /d "%~dp0"
call abaqus verify -user_exp
call abaqus job=forced_brittle_check input=forced_brittle.inp user=vumat_grind.for double=both cpus=1 datacheck
if errorlevel 1 exit /b 1
call abaqus job=forced_brittle input=forced_brittle.inp user=vumat_grind.for double=both cpus=8 interactive
