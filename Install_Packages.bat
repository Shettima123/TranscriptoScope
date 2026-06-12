@echo off
setlocal
cd /d "%~dp0"

powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0scripts\windows_install_packages.ps1"
if errorlevel 1 (
  echo.
  echo Package setup failed. Review the message above.
  pause
  exit /b 1
)

echo.
echo Package setup complete.
pause
