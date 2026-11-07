@echo off
setlocal
set OPENBLAS_NUM_THREADS=1
set OMP_NUM_THREADS=1
"%~dp0..\runtime\python\python.exe" -X utf8 "%~dp0run_all.py" %*
exit /b %errorlevel%
