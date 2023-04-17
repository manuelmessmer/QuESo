@echo off
set arg1=%1

IF "%arg1%" == "all" ( call py tibra/tests/test_tibra.py )
IF "%arg1%" == "py" ( call py tibra/tests/test_tibra.py )

IF "%arg1%" == "all" ( call "TestExecutables/Release/test_tibra.exe" )
IF "%arg1%" == "cpp" ( call "TestExecutables/Release/test_tibra.exe" )

IF "%arg1%" == "thingi10k" ( call "TestExecutables/Release/test_tibra_thingi10k.exe -- all %arg2% 5 50" )