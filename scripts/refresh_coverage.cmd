@echo off
rem Brings every heavy test's verdict and coverage up to date, starting at once (hours; see
rem scripts\refresh_coverage.py). Double-click, or run with arguments for it.
cd /d "%~dp0.."
python scripts\refresh_coverage.py %*
pause
