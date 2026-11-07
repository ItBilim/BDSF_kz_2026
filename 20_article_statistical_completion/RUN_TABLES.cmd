@echo off
setlocal
set OPENBLAS_NUM_THREADS=1
set OMP_NUM_THREADS=1
"%~dp0..\runtime\python\python.exe" -X utf8 "%~dp0build_article_completion.py"
exit /b %errorlevel%
