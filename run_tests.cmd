@echo off
set arg1=%1
set arg2=%2

IF "%arg1%" == "py" ( call py queso/tests/test_queso.py )

IF "%arg1%" == "cpp" ( call "TestExecutables/Release/test_queso.exe" )

IF "%arg1%" == "thingi10k" (
    call "TestExecutables/Release/test_queso_thingi10k.exe" --run_test=Thingi10KTestSuite/ElementClassificationTest -- small_set %arg2% 5 50
    call "TestExecutables/Release/test_queso_thingi10k.exe" --run_test=Thingi10KTestSuite/ElementClassificationTest -- small_set %arg2% 2 10
    call "TestExecutables/Release/test_queso_thingi10k.exe" --run_test=Thingi10KTestSuite/STLEmbeddingTest -- large_set %arg2% 5 50
    call "TestExecutables/Release/test_queso_thingi10k.exe" --run_test=Thingi10KTestSuite/STLEmbeddingTest -- large_set %arg2% 2 10 )
