@echo off
rem Guo 2019 1T base, law 2 (Swegle) on both fiber contacts, mu = 0.1, V0 = 70, fiber s25 (no RATE), E = 78 GPa.
set OMP_NUM_THREADS=8
set D=C:\Users\Maria\AppData\Local\Temp\claude\c--dip-sph-SPH\9c93d159-f8d5-4f88-9269-5c9cdc59e2d2\scratchpad\local1T
set B=C:\dip\sph\SPH\.claude\worktrees\contact-damped\build-elog\guo2019_shot.exe
set F=C:\dip\sph\SPH\.claude\worktrees\contact-damped\experiments\2026-10-17-guo2019-1t-pinball-vbl\twaron_plastic_s25.fib
cd /d %D%
%B% --config=1T --level=base --frames=off --save-t=2e-05 --threads=8 --n-mon=10 --n-diag=100 --repeat=1 --no-bracket --t-max=0.0015 --fiber-inelastic=%F% --eps-break=0.05 --fiber-E=7.8e+10 --energy-log=100 --v0=70 --mu=0.1 --fs-contact=2 --ff-contact=2 --out=%D%\swegle_mu0.1_v70_E78 > %D%\swegle_mu0.1_v70_E78.log 2>&1
echo rc=%ERRORLEVEL% >> %D%\swegle_mu0.1_v70_E78.log
