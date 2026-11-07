@echo off
setlocal enabledelayedexpansion
rem Owner 2026-10-24: same as run_law3_dmg_mu03_v71.cmd but V0 = 79, starts at once. Guo 1T base, law 3 WITHOUT push-out, dt x0.5,
rem E = 80 GPa, mu = 0.3, fiber s25 + DAMAGE (full failure at 4 %%, DEPS_F = 0.00415); frames (with D) and energy log.
set D=C:\Users\Maria\AppData\Local\Temp\claude\c--dip-sph-SPH\9c93d159-f8d5-4f88-9269-5c9cdc59e2d2\scratchpad\local1T
set B=C:\dip\sph\SPH\.claude\worktrees\contact-damped\build-dmg\guo2019_shot.exe
set OMP_NUM_THREADS=8
cd /d %D%
echo wait !DATE! !TIME! > %D%\law3_dmg_mu03_v79.queue.log
echo start !DATE! !TIME! >> %D%\law3_dmg_mu03_v79.queue.log
%B% --config=1T --level=base --v0=79 --mu=0.3 --fs-contact=3 --ff-contact=3 --dt-scale=0.5 --fiber-inelastic=%D%\twaron_s25_dmg0.00415.fib --eps-break=0.05 --fiber-E=8e+10 --frames=off --save-t=2e-05 --threads=8 --n-mon=10 --n-diag=100 --repeat=1 --no-bracket --t-max=0.0015 --energy-log=100 --out=%D%\law3_dmg_mu03_v79 > %D%\law3_dmg_mu03_v79.log 2>&1
echo done rc=!ERRORLEVEL! !DATE! !TIME! >> %D%\law3_dmg_mu03_v79.queue.log
