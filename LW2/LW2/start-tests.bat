@echo off
echo Running tests...

set PROGRAM="%~1"

call :RunTest 1 "tests\test1.txt" "tests\output1.txt" "tests\result1.txt" "Complete graph"
call :RunTest 2 "tests\test2.txt" "tests\output2.txt" "tests\result2.txt" "Complete bipartite graph"
call :RunTest 3 "tests\test3.txt" "tests\output3.txt" "tests\result3.txt" "Tree"
call :RunTest 4 "tests\test4.txt" "tests\output4.txt" "tests\result4.txt" "Wheel graph"
call :RunTest 5 "tests\test5.txt" "tests\output5.txt" "tests\result5.txt" "From lecture"

echo Testing completed
exit /B 0

:RunTest
set TEST_NUM=%1
set INPUT=%2
set OUTPUT=%3
set EXPECTED=%4
set DESCRIPTION=%5

echo.
echo Test %TEST_NUM%: %DESCRIPTION%
%PROGRAM% %INPUT% %OUTPUT%

fc.exe %OUTPUT% %EXPECTED% >nul
if errorlevel 1 (
    echo [ERROR] Test %TEST_NUM% failed
) else (
    echo [SUCCESS] Test %TEST_NUM% passed
)

exit /B 0
