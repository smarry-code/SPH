@echo off
rem Guo 2019 1T base, law 1 (pinball) on both fiber contacts, mu = 0.1, V0 = 70, E = 78 GPa, s25 + RATE JC.
set OMP_NUM_THREADS=8
set D=C:\Users\Maria\AppData\Local\Temp\claude\c--dip-sph-SPH\9c93d159-f8d5-4f88-9269-5c9cdc59e2d2\scratchpad\local1T
cd /d %D%
C:\dip\sph\SPH\.claude\worktrees\contact-damped\build-elog\guo2019_shot.exe --config=1T --level=base --v0=70 --mu=0.1 --out=%D%\pinball_rate_mu0.1_v70_E78 --frames=off --save-t=2e-05 --threads=8 --n-mon=10 --n-diag=100 --repeat=1 --no-bracket --t-max=0.0015 --fiber-inelastic=%D%\twaron_plastic_s25_rate.fib --eps-break=0.05 --fiber-E=7.8e+10 --fs-contact=1 --ff-contact=1 --energy-log=100 > %D%\pinball_rate_mu0.1_v70_E78.log 2>&1
echo rc=%ERRORLEVEL% >> %D%\pinball_rate_mu0.1_v70_E78.log
