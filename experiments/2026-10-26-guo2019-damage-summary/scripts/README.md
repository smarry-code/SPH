# Скрипты серий и визуализаций Guo 2019 (повреждённость, законы контакта 3–4)

Перенесены из временной папки сессии 28.10.2026. Пути внутри скриптов — абсолютные пути ноутбука
(`C:\Users\Maria\AppData\Local\Temp\...\scratchpad`, `C:\dip\sph\SPH\.claude\worktrees\contact-damped`) и сервера
(`~/SPH_dmg/build-elog/guo2019_shot`, `~/vbl*`); на другом устройстве их нужно поправить.

| Папка | Что внутри |
|---|---|
| `server/` | запуск серий на сервере: подбор окна (`run_win1T_*`, `run_win2T_*`, `run_law3*_window_2T*`), отдельные выстрелы и очередь после заполнения диска (`run_requeue_disk.sh`). Запуск: `setsid nohup bash <скрипт> &` |
| `local/` | запуск выстрелов и серий на ноутбуке (Windows, `.cmd` через Планировщик заданий) |
| `launch_old/` | запуски серий до повреждённости (Свигл, pinball, E = 78 ГПа, регрессия) |
| `viz/` | визуализации: `anim_2side.py` (GIF ткани с двух сторон, вся ткань), `anim_vz.py` (v_z(t)), `plot_energy_all.py` (энергия), `summary_dt05.py` (сводка серии), `ball_defl.py` (итоговый прогиб под шаром), `break_timing.py` / `break_where.py` (механизм разрывов); `make_viz_*.sh` / `redo_*.sh` — пакетные прогоны; `fetch_run.sh` — скачать выстрел с сервера (логи и каждый второй кадр) |

Драйверы подбора окна — `experiments/2026-10-17-guo2019-1t-pinball-vbl/scripts/window_bisect.py` и `window_pinball.py`.
Материал нити с повреждённостью — `../twaron_s25_dmg0.00415.fib` (полный разрыв при 4 %, E = 80 ГПа).
