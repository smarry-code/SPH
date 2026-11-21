# Воспроизведение: BC6, `bc6_shot`, типы контакта 0–3

## 1. Окружение

claude-compute-8t: Ubuntu 26.04.1, ядро 7.0.0-38, Intel Xeon Cascadelake (8 ядер), gcc 15.2.0,
CMake 4.2.3, Ninja 1.13.2. OpenMP: 4 потока на прогон; BC2 и BC6 шли параллельно (8 потоков всего).

## 2. Коммит

`61e24839f4de57904f4141f97f877a11db18d683` — «cases/bc2, cases/bc6: экспериментальные постановки и утилиты bc2_shot, bc6_shot».

## 3. Diff

Нет. Архив [`repro-diff.zip`](repro-diff.zip) содержит только `HEAD.txt` и `MANIFEST.md`.

## 4. Команды

```bash
git -c core.autocrlf=false archive --format=tar.gz -o src.tar.gz 61e2483
mkdir -p ~/ceramics/src_61e2483 && tar -xzf src.tar.gz -C ~/ceramics/src_61e2483
cd ~/ceramics/src_61e2483
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
ninja -C build -j4 bc6_shot test_e2e_bc6
for c in 3 0 1 2; do
  OMP_NUM_THREADS=4 ./build/bc6_shot --contact-type=$c --threads=4 --out=$HOME/ceramics/runs \
      --checkpoint-every=500 --overwrite
done
# сводный график
P=$HOME/ceramics/plot/bc6; mkdir -p $P
for c in 0 1 2 3; do cp $HOME/ceramics/runs/bc6_dx200_c$c/curve.txt $P/bc6_curve_$c.txt; done
./build/test_e2e_bc6 --full --plot-from=$P --emit-plot=$HOME/ceramics/plot/bc6_contact_types.svg
```

Тип 3 запускался задачей автозапуска `ceramics_bc6`, типы 0–2 — `ceramics_bc6_c012`
(`~/ceramics/run_bc6_c012.sh`: те же команды по очереди, при наличии чекпоинта — `--restart`).
Рестартов не было.

PDF-версия графика (на сервере — окружение с matplotlib 3.11.2, Python 3.14):

```bash
python3 tools/ceramic_contact_plot.py --case=bc6 --runs=$HOME/ceramics/runs --dx-um=200 \
    --reference=cases/bc6/experiment/reference.txt --out=bc6_contact_types.pdf
```

## 5. Стоимость

Типы 0–2: 90,1 + 91,3 + 94,1 с на 4 потоках (параллельно с BC2); кадры на сервере ~0,8 ГБ на прогон, в журнал не взяты.

## 6. Проверка

`data/c0/shot.txt`: `steps=1632`, `calc_plateau1=182.778`, `calc_ratio12=0.6861`; `data/c2/shot.txt`: `calc_dips=2`, `calc_dip2_t=3.49339e-06`.
