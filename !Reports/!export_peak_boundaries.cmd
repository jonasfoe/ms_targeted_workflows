@echo off
setlocal
REM This file reapplies all peak boundaries to fix issues where after a reimport of rawdata, transition / peak intensity is not recalculated
REM Requires a default skyline "Peak Boundaries" report template to be available in skyline
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
	--report-format=csv ^
	--report-invariant ^
	--report-name="Peak Boundaries" ^
	--report-file="%wd%\id_all_peak_boundaries.csv"
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
