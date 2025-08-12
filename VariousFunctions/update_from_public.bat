@echo off



cd ..

set UPSTREAM_FOUND="false"
set IS_REGULAR_FORK="false";

FOR /F "tokens=1,2,3" %%i IN ('git remote -v') DO (
	if "%%i %%j"=="origin https://github.com/KULeuvenNeuromechanics/PredSim.git" (
		set IS_REGULAR_FORK="true";
	)
	if "%%i %%j %%k"=="upstream https://github.com/KULeuvenNeuromechanics/PredSim.git (fetch)" (
		set UPSTREAM_FOUND="true"
	)
)

if not(%UPSTREAM_FOUND%=="true") and not(%IS_REGULAR_FORK%=="true") (
	echo setting https://github.com/KULeuvenNeuromechanics/PredSim.git as upstream repository
	echo " git remote add upstream https://github.com/KULeuvenNeuromechanics/PredSim.git"
	echo.
)

if %IS_REGULAR_FORK%=="true" (
	git fetch origin
	git checkout master
	git pull origin
) else (
	git fetch upstream
	git checkout master_public
	git merge upstream/master -m "update from master"
	git push origin

	start https://github.com/KULeuvenNeuromechanics/PredSim_private/compare/master...master_public
)

