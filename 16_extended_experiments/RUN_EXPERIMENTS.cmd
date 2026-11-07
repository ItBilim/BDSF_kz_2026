@echo off
setlocal
set OPENBLAS_NUM_THREADS=1
set OMP_NUM_THREADS=1
"%~dp0..\runtime\python\python.exe" -X utf8 "%~dp0run_experiments.py" %*
if errorlevel 1 exit /b %errorlevel%
"%~dp0..\runtime\python\python.exe" -X utf8 "%~dp0verify_extension.py"
exit /b %errorlevel%
