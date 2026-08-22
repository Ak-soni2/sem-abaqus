@echo off
cd /d "%~dp0"
abaqus verify -user_explicit
abaqus job=forced_brittle input=forced_brittle.inp user=vumat_grind.for double=both cpus=1 datacheck
if errorlevel 1 exit /b 1
abaqus job=forced_brittle input=forced_brittle.inp user=vumat_grind.for double=both cpus=8 interactive
