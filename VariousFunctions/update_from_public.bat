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

rem +-------------------------------------------------------------------------------
rem | Defaults
rem | These do not need to be changed if your (private) fork is set up according to 
rem | the instructions in PredSim/Documentation/PrivateForkPredSim.md

rem | Name and url of the remote.
set PUBLIC_REPO_REMOTE_NAME=upstream
set PUBLIC_REPO_REMOTE_URL=https://github.com/KULeuvenNeuromechanics/PredSim.git

rem | Branch in the private repository that you want to update.
set PRIVATE_REPO_MASTER_BRANCH=master



rem +-------------------------------------------------------------------------------
rem | Show info in command window to inform user what is going on.
echo Running update_from_public.bat
echo.
echo This script updates the %PRIVATE_REPO_MASTER_BRANCH% branch of a (private) fork
echo of the PredSim repository based on the master branch of the public PredSim repo.
echo For more information, see PredSim/Documentation/PrivateForkPredSim.md.
echo.


rem +-------------------------------------------------------------------------------
rem | Detect how the repositories are set up.

rem | Go to PredSim/ and run `echo git remote -v` on the command line to list the
rem | online repos linked to the local repo. For each line returned by this command,
rem | get the  first 3 words (separated by spaces or tabs). These are then used to 
rem | determine whether the current repo is a fork of the public PredSim repo, and
rem | whether the public repo is linked as upstream.

set ORIGIN_URL=%PUBLIC_REPO_REMOTE_URL%
set FOUND_PUBLIC_REPO_REMOTE="false"

cd ..
for /F "tokens=1,2,3" %%i in ('git remote -v') do (
	
	if "%%i"=="origin" (
		set ORIGIN_URL=%%j
	)	
	
	if "%%i %%j %%k"=="%PUBLIC_REPO_REMOTE_NAME% %PUBLIC_REPO_REMOTE_URL% (fetch)" (
		set FOUND_PUBLIC_REPO_REMOTE="true"
	)
)


rem | Check whether the origin is a fork of the remote by comparing their url.
if %PUBLIC_REPO_REMOTE_URL%==%ORIGIN_URL% (
	set IS_FORK="false"
) else (
	set IS_FORK="true"
)


rem | Finish setting up the remote, if that wasn't done yet.
if %IS_FORK%=="true" if not(%FOUND_PUBLIC_REPO_REMOTE%=="true") (
	echo Setting "%PUBLIC_REPO_REMOTE_URL%" as %PUBLIC_REPO_REMOTE_NAME% repository
	git remote add %PUBLIC_REPO_REMOTE_NAME% %PUBLIC_REPO_REMOTE_URL%
	echo.
)

rem +-------------------------------------------------------------------------------
rem | Do what we actually wanted to do: update the local master branch.

if %IS_FORK%=="true" (
	rem | Intended use: this script is run from a local clone of a (private) fork of the public repo.
	
	rem | Perform the steps described in PredSim/Documentation/PrivateForkPredSim.md#manual
	rem | 1.
	git fetch %PUBLIC_REPO_REMOTE_NAME%
	rem | 2.
	git checkout -B master_public %PRIVATE_REPO_MASTER_BRANCH%
	git push -u origin master_public
	rem | 3.
	git merge %PUBLIC_REPO_REMOTE_NAME%/master -m "update from public master"
	rem | 4.
	git push origin
	rem | 5.
	echo.
	echo.
	echo Create a pull request to finish updating the local %PRIVATE_REPO_MASTER_BRANCH% branch.
	start %ORIGIN_URL:~0,-4%/compare/%PRIVATE_REPO_MASTER_BRANCH%...master_public
	
) else (
	rem | Unintended use: this script is run from a local clone of the public repo.

	rem | Rather than throwing an error or doing nothing, do what the (confused) user expects to happen:
	rem | Update the local master branch based on the master branch of the public repo.
	git fetch origin
	git checkout master
	git pull origin
	echo.
	echo Done
)


rem | Keep the command window open until user closes it
pause
