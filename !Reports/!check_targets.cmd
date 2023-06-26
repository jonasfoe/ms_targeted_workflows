@echo off
setlocal
REM This file checks whether all PRM targets in skyline match those in the rawfiles compound annotation
REM Generates / requires skyline id_all_precursors.csv export
REM Author: Jonas D. Förster

set wd=%CD%
echo Working in %wd%
:RProj_search
    if exist ".\<.Rproj" (
        echo Found project root ^(.RProj^) at %CD%
    ) else (
		if "%CD%" == "%CD:~0,3%" (
			echo Can't find the project root ^(.RProj^) in directory tree.
			set evalstatus=1
			goto :STOP
		)
		cd ..
        goto :RProj_search
	)

echo.
"!Skyline Templates\SkylineRunner.exe" ^
	--in="%wd%\id_all.sky" ^
	--report-add="!Skyline Templates\id_all_precursors.skyr" ^
	--report-conflict-resolution=overwrite ^
	--report-format=csv ^
	--report-invariant ^
	--report-name="id_all_precursors" ^
	--report-file="%wd%\id_all_precursors.csv"

echo.
set scriptfile="!Reports/check_targets.R"
echo Running %scriptfile%
"Rscript.exe" -e "root <- getwd(); setwd(commandArgs(TRUE)[1]); source(file.path(root, '%scriptfile%'));" "%wd%"
set evalstatus=%ERRORLEVEL%

echo.
echo Finished.
:STOP
IF "%1" NEQ "--close_on_finish" (
	echo Press any key to close this window.
	pause >nul
)
endlocal
exit /b %evalstatus%
