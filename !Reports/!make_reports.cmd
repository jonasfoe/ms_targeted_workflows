@echo off
setlocal
REM This file generates reports based on the *_Report_<i>.yaml
REM Requires skyline id_all.csv (!Skyline Templates/id_all.skyr) and chromatogram_raw_id.csv (!Skyline Templates/chromatogram_raw_id.skyr) export
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

"!Skyline Templates\SkylineRunner.exe" ^
	--in="%wd%\id_all.sky" ^
	--report-add="!Skyline Templates\id_all_chromatograms.skyr" ^
	--report-conflict-resolution=overwrite ^
	--report-format=csv ^
	--report-invariant ^
	--report-name="id_all_chromatograms" ^
	--report-file="%wd%\id_all_chromatograms.csv"

echo.
set scriptfile="!Reports/gen_reports.R"
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
