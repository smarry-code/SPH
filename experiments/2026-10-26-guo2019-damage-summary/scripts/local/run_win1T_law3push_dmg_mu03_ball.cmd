@echo off
rem Owner 2026-10-25: continue the 1T window search (law 3 + push-out tau 0.5 us + damage 4 %%, mu 0.3, [65; 82]) with
rem the verdict by the ball only (as for 2T, no tolerance on broken yarns; --verdict-ball). Waits until the first driver
rem (pair 65 / 79, criterion <= 8 / 0) writes "done rc" to driver.out; finished shots are reused, verdicts recomputed.
set PYTHONIOENCODING=utf-8
cd /d C:\Users\Maria\AppData\Local\Temp\claude\c--dip-sph-SPH\9c93d159-f8d5-4f88-9269-5c9cdc59e2d2\scratchpad\local1T\win1T_law3push_dmg_mu03
:waitloop
find "done rc" C:\Users\Maria\AppData\Local\Temp\claude\c--dip-sph-SPH\9c93d159-f8d5-4f88-9269-5c9cdc59e2d2\scratchpad\local1T\win1T_law3push_dmg_mu03\driver.out >nul 2>&1
if errorlevel 1 (
  ping -n 61 127.0.0.1 >nul
  goto waitloop
)
echo restart with --verdict-ball %DATE% %TIME% >> C:\Users\Maria\AppData\Local\Temp\claude\c--dip-sph-SPH\9c93d159-f8d5-4f88-9269-5c9cdc59e2d2\scratchpad\local1T\win1T_law3push_dmg_mu03\driver2.out
"C:\Program Files\Python311\python.exe" C:\dip\sph\SPH\.claude\worktrees\contact-damped\experiments\2026-10-17-guo2019-1t-pinball-vbl\scripts\window_bisect.py --work=C:\Users\Maria\AppData\Local\Temp\claude\c--dip-sph-SPH\9c93d159-f8d5-4f88-9269-5c9cdc59e2d2\scratchpad\local1T\win1T_law3push_dmg_mu03 --bin=C:\dip\sph\SPH\.claude\worktrees\contact-damped\build-dmg\guo2019_shot.exe --fib=C:\Users\Maria\AppData\Local\Temp\claude\c--dip-sph-SPH\9c93d159-f8d5-4f88-9269-5c9cdc59e2d2\scratchpad\local1T\twaron_s25_dmg0.00415.fib --fiber-E=8e10 --eps-break=0.05 --config=1T --verdict-ball --contact=3 --mu=0.3 --vmin=65 --v-first=79 --vmax=82 --tol=1.5 --halves --max-rounds=6 "--extra=--dt-scale=0.5 --impulse-push-tau=5e-7" --threads=6 --threads-single=12 --t-max=1.5e-3 --t-max-long=3e-3 --save-t=2e-5 --energy-log=100 >> C:\Users\Maria\AppData\Local\Temp\claude\c--dip-sph-SPH\9c93d159-f8d5-4f88-9269-5c9cdc59e2d2\scratchpad\local1T\win1T_law3push_dmg_mu03\driver2.out 2>&1
echo done rc=%ERRORLEVEL% %DATE% %TIME% >> C:\Users\Maria\AppData\Local\Temp\claude\c--dip-sph-SPH\9c93d159-f8d5-4f88-9269-5c9cdc59e2d2\scratchpad\local1T\win1T_law3push_dmg_mu03\driver2.out
