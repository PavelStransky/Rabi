@echo off

for /L %%i in (%1,1,%2) do (
    echo Launching Schnellbruder.jl with parameters %%i %3
    julia Schnellbruder.jl %%i %3
)

pause
