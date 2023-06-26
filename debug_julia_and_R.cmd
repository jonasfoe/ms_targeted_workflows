@echo off
setlocal
setlocal enabledelayedexpansion

ECHO =======================================================
ECHO For RCall to work properly, R_HOME environment variable needs to be set properly
ECHO =======================================================
:ask_julia
ECHO This will recover the julia package library for this project.
ECHO Continue? ^(Y/N^)
set /P INPUT=Type input: 
If /I "!INPUT!"=="y" goto yes_julia
If /I "!INPUT!"=="n" goto no
ECHO Incorrect input & goto ask_julia
:yes_julia
ECHO =======================================================
ECHO Testing RCall from julia:
ECHO =======================================================
ECHO.
julia --project=. ^
  -e "using Pkg" ^
  -e "Pkg.status()" ^
  -e "Pkg.instantiate()" ^
  -e "Pkg.build(\"RCall\")" ^
  -e "using RCall" ^
  -e "println(R\"R.version.string\")"
set evalstatus=%ERRORLEVEL%

ECHO =======================================================
ECHO For juliacall to work properly, julia.exe needs to be in PATH
ECHO =======================================================
:ask_r
ECHO This will recover the R package library for this project.
ECHO Continue? ^(Y/N^)
set /P INPUT=Type input: 
If /I "!INPUT!"=="y" goto yes_r
If /I "!INPUT!"=="n" goto no
ECHO Incorrect input & goto ask_r
:yes_r
ECHO =======================================================
ECHO Testing JuliaCall from R:
ECHO =======================================================
ECHO.
RScript ^
  -e "renv::project()" ^
  -e "renv::restore()" ^
  -e "renv::rebuild('JuliaCall', recursive = FALSE, prompt = FALSE);" ^
  -e "library(JuliaCall); Sys.setenv(JULIA_PROJECT = normalizePath('.')); julia_setup();" ^
  -e "julia_library('Pkg'); julia_eval('Pkg.status()');"
set evalstatus=%ERRORLEVEL%

:no

"Rscript.exe" -e "stopifnot(MsRawAccess::check_license_accepted())"
IF %ERRORLEVEL% NEQ 0 (
	echo The license for the Thermo RawfileReader included with MsRawAccess needs to be check_license_accepted.
	"Rscript.exe" -e "MsRawAccess::print_license()"
	:ask_license
	echo Accept the license? ^(Y/N^)
	set /P INPUT=Type input: 
	If /I "!INPUT!"=="y" goto yes_license
	If /I "!INPUT!"=="n" goto no_license
	echo Incorrect input & goto ask_license
	:yes_license
	echo Registering that the license was accepted on this system.
	"Rscript.exe" -e "MsRawAccess::accept_license()"
)
:no_license

echo.
echo Finished.
IF "%1" NEQ "--close_on_finish" (
	echo Press any key to close this window.
	pause >nul
)
endlocal
exit /b %evalstatus%
