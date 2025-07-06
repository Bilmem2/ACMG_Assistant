@echo off
REM ================================================================
REM ACMG Variant Classification Assistant - Build Script
REM ================================================================
REM This script builds a standalone executable for Windows
REM No Python installation required for end users
REM ================================================================

echo.
echo 🧬 ACMG Variant Classification Assistant
echo ========================================
echo Building standalone executable (test mode disabled)...
echo.

REM Check if Python is available
python --version >nul 2>&1
if %errorlevel% neq 0 (
    echo ❌ Python is not installed or not in PATH
    echo Please install Python 3.8+ and try again
    pause
    exit /b 1
)

echo ✅ Python is available

REM Install build requirements
echo.
echo 📦 Installing build requirements...
pip install -r requirements_build.txt

if %errorlevel% neq 0 (
    echo ❌ Failed to install build requirements
    pause
    exit /b 1
)

echo ✅ Build requirements installed

REM Run the build script
echo.
echo 🔨 Building executable...
python build_executable.py

if %errorlevel% neq 0 (
    echo ❌ Build failed
    pause
    exit /b 1
)

echo.
echo 🎉 Build completed successfully!
echo.
echo 📂 Distribution folder: ACMG_Assistant_Portable
echo 💾 Executable size: ~50-100 MB
echo.
echo 📋 Next steps:
echo 1. Test the executable in ACMG_Assistant_Portable folder
echo 2. Zip the folder for distribution
echo 3. Share with users who don't have Python
echo.
echo Press any key to exit...
pause >nul
