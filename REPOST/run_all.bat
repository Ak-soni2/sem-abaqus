@echo off
rem =====================================================================
rem  Re-post-process all six completed jobs. NO Abaqus solve -- this only
rem  reads .odb files you already have. Minutes, not hours.
rem
rem  USAGE -- run it FROM the folder that holds the .odb files. Every
rem  output is written next to its .odb, so the working directory has to
rem  be that folder or the comparison table at the end finds nothing.
rem       cd /d D:\temp
rem       REPOST\run_all.bat
rem
rem  WHY RE-RUN AT ALL. The first pass used an older copy of the script
rem  sitting in D:\temp. Two things were wrong with it:
rem
rem    1. It divided by the MEAN element volume, 2.700e-11 mm3, on a graded
rem       mesh where the elements that actually get removed are the surface
rem       layer at 1.3846e-11 mm3. Every specific energy came out 1.95x low.
rem
rem    2. It never read the SDVs, so the ductile/brittle branch map -- the
rem       one figure this whole model exists to produce -- was never
rem       extracted.
rem
rem  AND ONE BUG IN HOW IT WAS CALLED. find_report() falls back to "first
rem  *_report.json in the folder", so with energy_criterion_report.json
rem  sitting in D:\temp all six jobs were resolved against it. Five of six
rem  therefore used the wrong dc (87.75 nm vs 52.92 nm), the wrong layer
rem  geometry and the wrong predicted split. The forces and energies come
rem  straight off the .odb and were unaffected, but the removal and branch
rem  numbers were not. Every line below passes its own report explicitly.
rem =====================================================================
setlocal
set PP=%~dp0postprocess_odb.py
set HS=%~dp0hotspot.py
set R=%~dp0reports

if not exist "single_abrasive1.odb" if not exist "single_abrasive2_sic.odb" (
  echo ERROR: no .odb files in the current directory.
  echo Run this FROM the folder holding the .odb files:
  echo     cd /d D:\temp
  echo     REPOST\run_all.bat
  exit /b 1
)

echo.
echo ############ sandstone ############
call :job single_abrasive1
call :job multi_abrasive1
call :job energy_abrasive1
echo.
echo ############ silicon carbide ############
call :job single_abrasive2_sic
call :job MULTI_abrasive1_sic
call :job ENERGY_abrasive1_sic

echo.
echo ############ comparison ############
abaqus python "%HS%" --table
echo.
echo Done. Send back: the console text above, plus *_summary.json,
echo *_hotspot.json and *_sdv.csv.
goto :eof

:job
if not exist "%~1.odb" (
  echo   SKIP %~1 -- no %~1.odb in this folder
  goto :eof
)
echo.
echo --- %~1 ---
abaqus python "%PP%" "%~1.odb" "%R%\%~1_report.json"
abaqus python "%HS%" "%~1.odb"
goto :eof
