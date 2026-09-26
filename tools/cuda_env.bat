@echo off
rem tools\cuda_env.bat -- MSVC environment for the CUDA build (fiber/CUDA_PLAN.md, stage G0).
rem On Windows nvcc compiles the host part of .cu files only through MSVC cl, and cl is
rem not in PATH outside the vcvars environment. The CPU build (preset clang-windows,
rem directory build) does not depend on this script.
rem
rem Usage -- wraps any command:
rem     tools\cuda_env.bat cmake --preset cuda
rem     tools\cuda_env.bat cmake --build build-cuda
rem     tools\cuda_env.bat ctest --test-dir build-cuda -L cuda
rem
rem vcvars64.bat is located via vswhere; fallback is the VS 2026 BuildTools path;
rem manual override -- environment variable VCVARS64.
rem The file is kept ASCII with CRLF: cmd.exe reads .bat in the OEM code page and
rem mis-parses UTF-8 text and LF line endings.
setlocal enabledelayedexpansion

if "%~1"=="" (
    echo [cuda_env] usage: tools\cuda_env.bat ^<command^> [args...]
    exit /b 1
)

set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
set "VCVARS="

if defined VCVARS64 set "VCVARS=%VCVARS64%"

if not defined VCVARS if exist "%VSWHERE%" (
    for /f "usebackq tokens=*" %%I in (`"%VSWHERE%" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do (
        if exist "%%I\VC\Auxiliary\Build\vcvars64.bat" set "VCVARS=%%I\VC\Auxiliary\Build\vcvars64.bat"
    )
)

if not defined VCVARS (
    set "VCVARS=%ProgramFiles(x86)%\Microsoft Visual Studio\18\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
)

if not exist "%VCVARS%" (
    echo [cuda_env] vcvars64.bat not found ^(looked for: "%VCVARS%"^).
    echo [cuda_env] Install MSVC Build Tools ^(nvcc needs cl^) or set VCVARS64.
    exit /b 1
)

rem vcvars64.bat may print a harmless "'vswhere.exe' is not recognized..." from its own
rem internals; the resulting environment (cl, INCLUDE, LIB) is still correct.
call "%VCVARS%" >nul
if errorlevel 1 (
    echo [cuda_env] vcvars64.bat failed.
    exit /b 1
)

%*
exit /b %errorlevel%
