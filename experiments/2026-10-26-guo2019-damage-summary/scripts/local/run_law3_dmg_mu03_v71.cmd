@echo off
setlocal enabledelayedexpansion
rem Owner 2026-10-23: Guo 1T base, law 3 WITHOUT push-out on fiber-fiber and fiber-solid, dt x0.5, E = 80 GPa,
rem mu = 0.3, V0 = 71, fiber s25 + DAMAGE (full failure at 4 %%, DEPS_F = 0.00415); frames (with D column) and energy log.
rem Waits until the 2T 150 m/s shot (guo2019_shot.exe) finishes, then runs on 8 threads.
set D=C:\Users\Maria\AppData\Local\Temp\claude\c--dip-sph-SPH\9c93d159-f8d5-4f88-9269-5c9cdc59e2d2\scratchpad\local1T
set B=C:\dip\sph\SPH\.claude\worktrees\contact-damped\build-dmg\guo2019_shot.exe
set OMP_NUM_THREADS=8
cd /d %D%
echo wait !DATE! !TIME! > %D%\law3_dmg_mu03_v71.queue.log
:waitloop
tasklist /FI "IMAGENAME eq guo2019_shot.exe" | find /I "guo2019_shot" >nul
if not errorlevel 1 (
  ping -n 121 127.0.0.1 >nul
  goto waitloop
)
echo start !DATE! !TIME! >> %D%\law3_dmg_mu03_v71.queue.log
%B% --config=1T --level=base --v0=71 --mu=0.3 --fs-contact=3 --ff-contact=3 --dt-scale=0.5 --fiber-inelastic=%D%\twaron_s25_dmg0.00415.fib --eps-break=0.05 --fiber-E=8e+10 --frames=off --save-t=2e-05 --threads=8 --n-mon=10 --n-diag=100 --repeat=1 --no-bracket --t-max=0.0015 --energy-log=100 --out=%D%\law3_dmg_mu03_v71 > %D%\law3_dmg_mu03_v71.log 2>&1
echo done rc=!ERRORLEVEL! !DATE! !TIME! >> %D%\law3_dmg_mu03_v71.queue.log
