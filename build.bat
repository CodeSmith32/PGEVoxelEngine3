@echo off
REM		Build Script

REM Set Compiler Settings Here

set CPP=c++
set GPP=g++
set WINRES=windres
set ARCH=64
set OUTPUT=program.exe
set DEBUGMODE=1


del %OUTPUT% 2>nul
del .objs64\application.o 2>nul

setlocal enabledelayedexpansion

if not exist .objs64 (
	echo Creating Object Directory Structure...
	mkdir .objs64
)

echo Building Game Files...
for %%F in (*.cpp) do (
	if not exist .objs64\%%~nF.o (
		echo Building %%~nF.o
		start /B "%%~nF.o" %CPP% -std=c++20 -c %%F -o .objs64\%%~nF.o
	)
)

REM Wait for building process to finish
:loop
for /f %%G in ('tasklist ^| find /c "%CPP%"') do ( set count=%%G )
if %count%==0 (
	goto linker
) else (
	timeout /t 2 /nobreak>nul
	goto loop
)

:linker

set "files="
for /f "delims=" %%A in ('dir /b /a-d ".objs64\%*" ') do set "files=!files! .objs64\%%A"

echo Linking Executable...
if %ARCH%==64 (
	goto link
)
if %ARCH%==32 (
	goto link
)
echo ARCH Must be 32 or 64! Make sure ARCH matches the compiler's architecture!
goto finish
:link

if %DEBUGMODE% GTR 0 (
	set MWINDOWS=
) else (
	set MWINDOWS=-mwindows
)

%GPP% -L.\library -o %OUTPUT% %files% -static-libstdc++ -lpthread -static -lsetupapi -lwinmm -luser32 -lgdi32 -lopengl32 -lgdiplus -lShlwapi -ldwmapi -lstdc++fs %MWINDOWS%

:finish
if exist .\%OUTPUT% (
	echo Build Success!
) else (
	echo Build Failed!
)

@REM pause>nul

start "" %OUTPUT%