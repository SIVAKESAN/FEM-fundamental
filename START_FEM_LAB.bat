@echo off
setlocal
cd /d "%~dp0"
title FEM Learning Laboratory

echo ============================================================
echo                  FEM LEARNING LABORATORY

echo         Python FEM solver + Black and White Website

echo ============================================================
echo.

if not exist ".venv\Scripts\python.exe" (
  echo [1/3] Creating Python virtual environment...
  python -m venv .venv
  if errorlevel 1 (
    echo.
    echo ERROR: Python is not available or virtual environment failed.
    echo Install Python 3.10+ and enable Add Python to PATH.
    pause
    exit /b 1
  )
)

echo [2/3] Checking/installing Python dependencies...
".venv\Scripts\python.exe" -m pip install -r requirements.txt
if errorlevel 1 (
  echo.
  echo ERROR: Could not install Python packages.
  echo Check your internet connection and try again.
  pause
  exit /b 1
)

echo.
echo [3/3] Starting Python FEM server on port 8010...
echo.
echo Open this URL in your browser:
echo.
echo               http://127.0.0.1:8010
echo.
echo Do NOT close this window while using the website.
echo Press CTRL+C to stop the server.
echo.
start "" /min cmd /c "timeout /t 6 /nobreak >nul & start http://127.0.0.1:8010"
".venv\Scripts\python.exe" -m uvicorn backend.app:app --host 127.0.0.1 --port 8010
if errorlevel 1 (
  echo.
  echo Cannot start port 8010. Make sure no other FEM Lab is running.
)
pause
