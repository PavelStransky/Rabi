@echo off

for /L %%i in (1,30,990) do (
    echo Launching Schnellbruder.jl with parameters %%i %1
    julia Schnellbruder.jl %%i %1
)

pause
