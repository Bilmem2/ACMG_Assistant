#!/usr/bin/env python3
"""
PyInstaller build script for ACMG Variant Classification Assistant
================================================================

This script creates a standalone executable that can run without Python installation.
The executable supports normal mode only (test mode requires Python).
"""

import os
import sys
import shutil
import subprocess
from pathlib import Path

def build_executable():
    """Build the standalone executable using PyInstaller"""
    
    print("🔧 Building ACMG Assistant Executable...")
    print("=" * 50)
    
    # Check if PyInstaller is installed
    try:
        import PyInstaller
        print(f"✅ PyInstaller version: {PyInstaller.__version__}")
    except ImportError:
        print("❌ PyInstaller not found. Installing...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", "pyinstaller"])
        print("✅ PyInstaller installed successfully")
    
    # Build configuration
    script_name = "acmg_assistant.py"
    exe_name = "ACMG_Assistant"
    dist_folder = "ACMG_Assistant_Portable"
    
    # PyInstaller command - Console app for better interaction
    cmd = [
        "pyinstaller",
        "--onefile",                    # Single executable file
        "--console",                    # Keep console window
        "--name", exe_name,             # Executable name
        "--add-data", "config;config",  # Include config folder
        "--add-data", "core;core",      # Include core folder
        "--add-data", "utils;utils",    # Include utils folder
        "--hidden-import", "colorama",  # Ensure colorama is included
        "--hidden-import", "requests",  # Ensure requests is included
        "--hidden-import", "scipy",     # Ensure scipy is included
        "--hidden-import", "numpy",     # Ensure numpy is included
        "--exclude-module", "test_*",   # Exclude test modules
        "--exclude-module", "matplotlib", # Exclude matplotlib
        "--exclude-module", "PIL",      # Exclude PIL
        "--exclude-module", "tkinter",  # Exclude tkinter
        script_name
    ]
    
    print(f"🔨 Building console executable (test mode disabled)...")
    print(f"Command: {' '.join(cmd)}")
    
    try:
        # Run PyInstaller
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print("✅ Build completed successfully!")
        
        # Create distribution folder
        create_distribution_folder(exe_name, dist_folder)
        
        print(f"\n🎉 Executable created successfully!")
        print(f"📁 Location: {dist_folder}")
        print(f"📝 Instructions:")
        print(f"   1. Test the executable in the {dist_folder} folder")
        print(f"   2. Zip the entire {dist_folder} folder for distribution")
        print(f"   3. Share with users who don't have Python installed")
        
    except subprocess.CalledProcessError as e:
        print(f"❌ Build failed: {e}")
        if e.stdout:
            print(f"stdout: {e.stdout}")
        if e.stderr:
            print(f"stderr: {e.stderr}")
        return False
    
    return True

def create_distribution_folder(exe_name, dist_folder):
    """Create a clean distribution folder with necessary files"""
    
    print(f"\n📦 Creating distribution folder: {dist_folder}")
    
    # Create distribution folder
    if os.path.exists(dist_folder):
        shutil.rmtree(dist_folder)
    os.makedirs(dist_folder)
    
    # Copy executable
    exe_path = f"dist/{exe_name}.exe"
    if os.path.exists(exe_path):
        shutil.copy2(exe_path, dist_folder)
        print(f"✅ Copied {exe_name}.exe")
    else:
        print(f"❌ Executable not found: {exe_path}")
        return False
    
    # Create batch files for easy use
    create_batch_files(dist_folder, exe_name)
    
    # Copy documentation
    copy_documentation(dist_folder)
    
    # Create README for executable
    create_executable_readme(dist_folder)
    
    print("✅ Distribution folder created successfully!")
    return True

def create_batch_files(dist_folder, exe_name):
    """Create batch files for easy usage"""
    
    # Normal mode batch file
    normal_batch = f"""@echo off
title ACMG Variant Classification Assistant
echo.
echo ======================================
echo  ACMG Variant Classification Assistant
echo ======================================
echo.
echo Starting normal mode...
echo Note: Test mode is only available in Python version
echo.
{exe_name}.exe
echo.
echo Press any key to exit...
pause >nul
"""
    
    with open(f"{dist_folder}/Run_ACMG_Assistant.bat", "w") as f:
        f.write(normal_batch)
    
    # ACMG 2023 mode batch file
    acmg2023_batch = f"""@echo off
title ACMG Assistant - 2023 Guidelines
echo.
echo =======================================
echo  ACMG Assistant - 2023 Guidelines
echo =======================================
echo.
echo Starting with ACMG 2023 guidelines...
echo.
{exe_name}.exe --acmg-2023
echo.
echo Press any key to exit...
pause >nul
"""
    
    with open(f"{dist_folder}/Run_ACMG_2023.bat", "w") as f:
        f.write(acmg2023_batch)
    
    print("✅ Created batch files for easy usage")

def copy_documentation(dist_folder):
    """Copy essential documentation files"""
    
    docs_to_copy = [
        "LICENSE",
        "README.md"
    ]
    
    for doc in docs_to_copy:
        if os.path.exists(doc):
            shutil.copy2(doc, dist_folder)
            print(f"✅ Copied {doc}")

def create_executable_readme(dist_folder):
    """Create a README file specifically for the executable"""
    
    readme_content = """# ACMG Variant Classification Assistant - Portable Version

## Quick Start

1. **Double-click** `Run_ACMG_Assistant.bat` to start the application
2. **Follow the interactive prompts** to analyze your variant
3. **Results** will be saved in the same folder

## Usage Options

- **Run_ACMG_Assistant.bat**: Standard mode with ACMG 2015 guidelines
- **Run_ACMG_2023.bat**: Enhanced mode with ACMG 2023 guidelines
- **ACMG_Assistant.exe**: Direct executable (advanced users)

## Test Mode

⚠️ **Test mode is not available in the executable version**
- Test mode requires Python installation and is intended for developers
- If you need test mode, please use the Python version of the tool

## Before You Start

### Required Information
Before running the tool, gather these scores from external resources:

1. **Varsome** (https://varsome.com) - Primary resource for in silico scores
2. **ClinVar** (https://www.ncbi.nlm.nih.gov/clinvar/) - Variant annotations
3. **gnomAD** (https://gnomad.broadinstitute.org/) - Population frequencies

### Supported Predictors
- **REVEL, CADD, AlphaMissense** (primary predictors)
- **VEST4, PrimateAI, ESM1b** (high priority predictors)
- **SIFT, PolyPhen-2, MutationTaster** (individual predictors)
- **SpliceAI, MMSplice** (splice prediction)
- **PhyloP, GERP++** (conservation scores)

## Important Notes

⚠️ **This tool is for research and educational purposes only**
- Not intended for direct clinical use without professional oversight
- Results must be validated by certified professionals
- Always cross-check with primary sources

## Support

- **GitHub**: https://github.com/Bilmem2/acmg-assessor
- **Author**: Can Sevilmiş
- **LinkedIn**: https://linkedin.com/in/cansevilmiss

For technical support, please create an issue on GitHub.

## License

This project is licensed under the MIT License - see LICENSE file for details.
"""
    
    with open(f"{dist_folder}/README_EXECUTABLE.txt", "w") as f:
        f.write(readme_content)
    
    print("✅ Created executable README")

def main():
    """Main build function"""
    
    if not os.path.exists("acmg_assistant.py"):
        print("❌ Error: acmg_assistant.py not found in current directory")
        sys.exit(1)
    
    success = build_executable()
    
    if success:
        print(f"\n🎉 Build completed successfully!")
        print(f"📁 Executable ready for distribution in 'ACMG_Assistant_Portable' folder")
        print(f"💡 Note: Test mode is disabled in executable - only available in Python version")
    else:
        print(f"\n❌ Build failed. Please check the error messages above.")
        sys.exit(1)

if __name__ == "__main__":
    main()
