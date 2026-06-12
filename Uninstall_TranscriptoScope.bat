@echo off
setlocal
cd /d "%~dp0"

powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0scripts\windows_uninstall.ps1"
if errorlevel 1 (
  echo.
  echo Uninstall failed. Review the message above.
  pause
  exit /b 1
)

echo.
pause
