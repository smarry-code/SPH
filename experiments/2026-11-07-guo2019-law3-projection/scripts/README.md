# Скрипты анализа серии «закон 3 как проекция скорости» (1T, V0=79, base)

Запускать только на `claude-compute` (python3.14 + numpy/scipy/matplotlib, установлены `--user`;
если нет: `curl -sS https://bootstrap.pypa.io/get-pip.py | python3 - --user --break-system-packages`,
затем `python3 -m pip install --user --break-system-packages numpy scipy matplotlib`).

| скрипт | вход | выход | команда |
|---|---|---|---|
| `common.py` | — | общие читатели (shot, monitor/energy с дедупликацией шагов, кадры, раскладка) | импортируется |
| `metrics.py` | `<root>/<RUN>/{shot.txt,monitor.csv,energy_log.csv}` | `metrics.json`, `tables.md` (таблицы + pass/FAIL по критериям ТЗ) | `python3 metrics.py ~/proj_analysis --out ~/proj_analysis_out` |
| `frames.py` | `<RUN>/{fiber_output_*,sph_output_*,guo2019_fabric_layout.txt,guo2019_fabric.fib,shot.txt}` | CSV: t, δ_ff/d, δ_fs/d, n_tunnel, проходы по кадрам; в stdout сверка с shot.txt | `python3 frames.py ~/proj_analysis/C2 --out ~/proj_analysis_out/frames/C2.csv` |
| `plots.py` | `<root>`, каталог CSV из frames.py | PNG (подписи по-русски) | `python3 plots.py ~/proj_analysis ~/proj_analysis_out/frames --out ~/proj_analysis_out/png` |
| `run_all.sh` | `~/proj_analysis` | всё выше в `~/proj_analysis_out` | `bash run_all.sh` |

## Данные

Источник: `claude-compute-8t:~/proj_runs/<RUN>/1T_base_v79.0_r0/`, RUN ∈ C1 C2 C3 C4 C5 C1_k1 C1_k2 C1_k3;
ключи запуска — `claude-compute-8t:~/proj_runs/run_<RUN>.sh`. Копируется плоско в `claude-compute:~/proj_analysis/<RUN>/`
без `checkpoint*` (и без `sph_gp.csv`); скрипты запуска — в `~/proj_analysis/run_scripts/`.

Канал mac -> серверы по скорости скачивания очень медленный (~10 кБ/с), поэтому копировать надо напрямую
сервер -> сервер с пробросом агента (на 8t только чтение, `nice`):

    ssh-add ~/.ssh/<ключ>          # локально
    ssh -A claude-compute 'for R in C1 C2 C3 C4 C5 C1_k1 C1_k2 C1_k3; do mkdir -p ~/proj_analysis/$R;
      ssh claude-compute-8t "cd ~/proj_runs/$R/1T_base_v79.0_r0 && nice -n 19 tar cf - --exclude=checkpoint* --exclude=sph_gp.csv ." \
        | tar xf - -C ~/proj_analysis/$R; done'

(агент должен жить всё время команды: не уводить её в `nohup`/фон на сервере.)

## Примечания к метрикам

- В monitor.csv / energy_log.csv при рестартах повторы шагов — берётся последнее вхождение.
- `W_proj_pos_max` в monitor.csv — НЕ накопленный, а максимум в окне между строками (немонотонен, в конце 0);
  итог прогона — `W_proj_pos_max` из shot.txt = максимум колонки. Ступени считаются по накопленному максимуму колонки.
- `frames.py` повторяет критерии `tests/guo2019_monitor.hpp::diagnostics()`; r_s = 0,5(m/max(rho, rho_0))^(1/3),
  rho_0 берётся из кадра шара первого шага. Кадры — подвыборка шагов, поэтому максимум по кадрам <= значения из shot.txt.
- Проходы (а) — по плану xy, см. докстринг `Crossings`; метрика ненадёжна около шара (крутые нити), даёт ложные
  срабатывания и в прогонах без проникания (C3) — использовать как справочную.
