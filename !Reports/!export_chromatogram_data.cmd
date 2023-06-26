@echo off
setlocal
REM This file exports id_all_precursors.csv, id_all_transitions.csv, id_all_chromatograms.csv from id_all.sky
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
