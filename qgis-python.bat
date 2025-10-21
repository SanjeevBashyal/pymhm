REM @ECHO OFF
REM ###########################################################################
REM #
REM #  Script to launch Visual Studio Code with the QGIS Python Environment
REM #
REM ###########################################################################

REM -- Set the QGIS version and root directory --
SET QGIS_VERSION=QGIS 3.40.6
SET QGIS_ROOT="C:\Program Files\%QGIS_VERSION%"

REM -- The folder you want to open in VS Code --
SET VSC_PROJECT_FOLDER="E:\0 Python\pymhm"

REM -- Call the QGIS environment setup batch file --
CALL %QGIS_ROOT%\bin\o4w_env.bat

REM -- (Optional) You can verify the paths are set correctly by un-commenting the next line --
REM path

REM -- Launch Cursor in the specified project folder --
ECHO "Starting Cursor in QGIS environment..."
START "" "%LOCALAPPDATA%\Programs\cursor\Cursor.exe" %VSC_PROJECT_FOLDER%

pause