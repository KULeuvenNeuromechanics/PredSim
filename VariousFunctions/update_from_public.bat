@echo off

rem +-------------------------------------------------------------------------------
rem | update_from_public.bat
rem |	This script updates the master branch of a (private) fork of the PredSim
rem |	repository based on the master branch of the public PredSim repository. 
rem |	To run this script, go to its location in File Explorer and double-click it
rem |
rem |	More information about...
rem |	1. creating a private fork of the public PredSim repository
rem |	2. (manually) updating the private repo with changes made in the public repo
rem |	3. contributing code from the private repo back to the public repo
rem |	...is available in PredSim/Documentation/PrivateForkPredSim.md
rem |
rem |
rem | Original author: Lars D'Hondt
rem | Original date: 13/August/2025
rem +-------------------------------------------------------------------------------

echo Running update_from_public.bat
echo.
echo This script updates the master branch of a (private) fork of the PredSim repository 
echo based on the master branch of the public PredSim repository. 
echo For more information, see PredSim/Documentation/PrivateForkPredSim.md.
echo.


rem % Go to PredSim/ and run `echo git remote -v` on the command line to list the online  
rem % repos linked to the local repo. For each line returned by this command, get the 
rem % first 3 words (separated by spaces or tabs). These are then used to determine whether 
rem % the current repo is a fork of the public PredSim repo, and whether the public repo
rem % is linked as upstream.

cd ..

rem % Default name and url of the remote
set PUBLIC_REPO_REMOTE_NAME=upstream
set PUBLIC_REPO_REMOTE_URL=https://github.com/KULeuvenNeuromechanics/PredSim.git


set ORIGIN_URL=%PUBLIC_REPO_REMOTE_URL%
set FOUND_PUBLIC_REPO_REMOTE="false"
set IS_FORK="false"

for /F "tokens=1,2,3" %%i in ('git remote -v') do (
	
	if "%%i"=="origin" (
		set ORIGIN_URL=%%j
	)
	
	if "%%i %%j %%k"=="%PUBLIC_REPO_REMOTE_NAME% %PUBLIC_REPO_REMOTE_URL% (fetch)" (
		set FOUND_PUBLIC_REPO_REMOTE="true"
	)
)


rem % Forked repos have a different url than the main PredSim repo.
if not("%PUBLIC_REPO_REMOTE_URL%"=="%ORIGIN_URL%") (set IS_FORK="true")

rem % Finish setting up the remote, if that wasn't done yet.
if %IS_FORK%=="true" if not(%FOUND_PUBLIC_REPO_REMOTE%=="true") (
	echo setting "%PUBLIC_REPO_REMOTE_URL%" as %PUBLIC_REPO_REMOTE_NAME% repository
	echo git remote add %PUBLIC_REPO_REMOTE_NAME% %PUBLIC_REPO_REMOTE_URL%
	echo.
)


rem % Do what we actually wanted to do: update the local master branch.
if %IS_FORK%=="true" (
	rem % Intended use: this script is run from a local clone of a (private) fork of the public repo.

	rem % Perform the steps described in PredSim/Documentation/PrivateForkPredSim.md#manual
	rem % 1.
	echo git fetch %PUBLIC_REPO_REMOTE_NAME%
	rem % 2.
	echo git checkout -B master_public
	echo git push -u origin master_public
	rem % 3.
	echo git merge %PUBLIC_REPO_REMOTE_NAME%/master -m "update from public master"
	rem % 4.
	echo git push origin
	rem % 5.
	start %ORIGIN_URL:~0,-4%/compare/master...master_public
	
) else (

	rem % Unintended use: this script is run from a local clone of the public repo.

	rem % Rather than throwing an error or doing nothing, do what the (confused) user expects to happen:
	rem % Update the local master branch baed on the master branch of the public repo.
	echo git fetch origin
	echo git checkout master
	echo git pull origin
)

rem % The End
echo Done
echo.

rem % Keep the command window open until user closes it
pause
