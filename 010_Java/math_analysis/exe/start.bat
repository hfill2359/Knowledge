@echo off

set CURRENT_DIR=%~dp0

set CLASSPATH=%CLASSPATH%;^
%CURRENT_DIR%/test.jar


java -Xms1024M -Xmx1024M jp.Test

pause
