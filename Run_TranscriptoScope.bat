@echo off
setlocal
cd /d "%~dp0"

wscript "%~dp0Launch_TranscriptoScope.vbs"
if errorlevel 1 (
  echo.
  echo App launch failed. Review the message above.
  pause
  exit /b 1
)
