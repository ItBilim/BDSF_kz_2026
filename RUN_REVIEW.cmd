@echo off
setlocal
"%~dp0runtime\python\python.exe" -I -X utf8 "%~dp013_reproducibility_package\13_pipeline.py" --project-root "%~dp0."
exit /b %errorlevel%
