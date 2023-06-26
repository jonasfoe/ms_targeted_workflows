@echo off
setlocal
REM This file generate sequence coverage plots for all precursors and replicates
REM Generates / requires skyline id_all_precursors.csv and id_all_transitions exports
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
"!Skyline Templates\SkylineRunner.exe" ^
	--in="%wd%\id_all.sky" ^
	--report-add="!Skyline Templates\id_all_transitions.skyr" ^
	--report-conflict-resolution=overwrite ^
	--report-format=csv ^
	--report-invariant ^
	--report-name="id_all_transitions" ^
	--report-file="%wd%\id_all_transitions.csv"

echo.
set scriptfile="!Reports/gen_sequence_coverage_plots.R"
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
