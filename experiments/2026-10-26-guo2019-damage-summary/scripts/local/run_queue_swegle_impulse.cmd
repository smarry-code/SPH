@echo off
rem Local queue, Guo 2019 1T base, fiber s25 (PLASTIC SY=2.04e9 H=100e9, FAILURE SIGMA_U=2.5e9, no RATE/VISCO),
rem E = 78 GPa, eps_break = 0.05, 8 threads, t_max = 1.5 ms, frames every 20 us, energy log every 100 steps.
rem  1) law 2 (Swegle) on both fiber contacts, mu = 0.3, V0 = 70
rem  2) law 3 (impulse, gamma = 1) on both fiber contacts, mu = 0.1, V0 = 73
set OMP_NUM_THREADS=8
set D=C:\Users\Maria\AppData\Local\Temp\claude\c--dip-sph-SPH\9c93d159-f8d5-4f88-9269-5c9cdc59e2d2\scratchpad\local1T
set B=C:\dip\sph\SPH\.claude\worktrees\contact-damped\build-elog\guo2019_shot.exe
set F=C:\dip\sph\SPH\.claude\worktrees\contact-damped\experiments\2026-10-17-guo2019-1t-pinball-vbl\twaron_plastic_s25.fib
set COMMON=--config=1T --level=base --frames=off --save-t=2e-05 --threads=8 --n-mon=10 --n-diag=100 --repeat=1 --no-bracket --t-max=0.0015 --fiber-inelastic=%F% --eps-break=0.05 --fiber-E=7.8e+10 --energy-log=100
cd /d %D%
echo start swegle %DATE% %TIME% > %D%\queue_swegle_impulse.log
%B% %COMMON% --v0=70 --mu=0.3 --fs-contact=2 --ff-contact=2 --out=%D%\swegle_mu0.3_v70_E78 > %D%\swegle_mu0.3_v70_E78.log 2>&1
echo swegle rc=%ERRORLEVEL% %DATE% %TIME% >> %D%\queue_swegle_impulse.log
%B% %COMMON% --v0=73 --mu=0.1 --fs-contact=3 --ff-contact=3 --out=%D%\impulse_mu0.1_v73_E78 > %D%\impulse_mu0.1_v73_E78.log 2>&1
echo impulse rc=%ERRORLEVEL% %DATE% %TIME% >> %D%\queue_swegle_impulse.log
