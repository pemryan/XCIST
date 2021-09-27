ECHO OFF
REM Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE
ECHO ON

set MINGW_PATH="C:\mingw64\bin"

set PATH=%MINGW_PATH%\bin;%PATH%
CD src

%MINGW_PATH%\bin\mingw32-make -f ..\MakeWindows64
move /Y libcatsim64.dll ..\..\catsim\lib
@PAUSE

%MINGW_PATH%\bin\mingw32-make -f ..\MakeWindows64 clean

CD ..

@ECHO Windows build complete
@PAUSE
