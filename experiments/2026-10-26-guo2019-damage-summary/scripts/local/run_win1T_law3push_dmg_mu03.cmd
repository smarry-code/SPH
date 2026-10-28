@echo off
rem Owner 2026-10-25: window search Guo 1T base, law 3 WITH push-out (tau 0.5 us, owner's explicit request), dt x0.5,
rem mu = 0.3, fiber s25 + DAMAGE (full failure at 4 %%), E = 80 GPa; bisection in [65; 82] m/s (first pair 65 and 79;
rem at mu 0.25 with push 71 was already perforated, V_R 12.4). Criterion 1T: stopped, V_R = 0, broken <= 8 / 0.
rem Pair x 6 threads, single shot x 12 threads; frames 20 us, energy log.
set PYTHONIOENCODING=utf-8
cd /d C:\Users\Maria\AppData\Local\Temp\claude\c--dip-sph-SPH\9c93d159-f8d5-4f88-9269-5c9cdc59e2d2\scratchpad\local1T\win1T_law3push_dmg_mu03
"C:\Program Files\Python311\python.exe" C:\dip\sph\SPH\.claude\worktrees\contact-damped\experiments\2026-10-17-guo2019-1t-pinball-vbl\scripts\window_bisect.py --work=C:\Users\Maria\AppData\Local\Temp\claude\c--dip-sph-SPH\9c93d159-f8d5-4f88-9269-5c9cdc59e2d2\scratchpad\local1T\win1T_law3push_dmg_mu03 --bin=C:\dip\sph\SPH\.claude\worktrees\contact-damped\build-dmg\guo2019_shot.exe --fib=C:\Users\Maria\AppData\Local\Temp\claude\c--dip-sph-SPH\9c93d159-f8d5-4f88-9269-5c9cdc59e2d2\scratchpad\local1T\twaron_s25_dmg0.00415.fib --fiber-E=8e10 --eps-break=0.05 --config=1T --contact=3 --mu=0.3 --vmin=65 --v-first=79 --vmax=82 --tol=1.5 --halves --max-rounds=6 "--extra=--dt-scale=0.5 --impulse-push-tau=5e-7" --threads=6 --threads-single=12 --t-max=1.5e-3 --t-max-long=3e-3 --save-t=2e-5 --energy-log=100 > C:\Users\Maria\AppData\Local\Temp\claude\c--dip-sph-SPH\9c93d159-f8d5-4f88-9269-5c9cdc59e2d2\scratchpad\local1T\win1T_law3push_dmg_mu03\driver.out 2>&1
echo done rc=%ERRORLEVEL% %DATE% %TIME% >> C:\Users\Maria\AppData\Local\Temp\claude\c--dip-sph-SPH\9c93d159-f8d5-4f88-9269-5c9cdc59e2d2\scratchpad\local1T\win1T_law3push_dmg_mu03\driver.out
