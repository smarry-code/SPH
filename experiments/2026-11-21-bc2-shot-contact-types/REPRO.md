# Воспроизведение: BC2, `bc2_shot`, типы контакта 0–3

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
ninja -C build -j4 bc2_shot test_e2e_bc2
for c in 3 0 1 2; do
  OMP_NUM_THREADS=4 ./build/bc2_shot --contact-type=$c --threads=4 --out=$HOME/ceramics/runs \
      --checkpoint-every=500 --overwrite
done
# сводный график
P=$HOME/ceramics/plot/bc2; mkdir -p $P
for c in 0 1 2 3; do cp $HOME/ceramics/runs/bc2_dx200_c$c/curve.txt $P/bc2_curve_$c.txt; done
./build/test_e2e_bc2 --full --plot-from=$P --emit-plot=$HOME/ceramics/plot/bc2_contact_types.svg
```

Тип 3 запускался задачей автозапуска `ceramics_bc2`, типы 0–2 — `ceramics_bc2_c012`
(`~/ceramics/run_bc2_c012.sh`: те же команды по очереди, при наличии чекпоинта — `--restart`).
Рестартов не было.

PDF-версия графика (на сервере — окружение с matplotlib 3.11.2, Python 3.14):

```bash
python3 tools/ceramic_contact_plot.py --case=bc2 --runs=$HOME/ceramics/runs --dx-um=200 \
    --reference=cases/bc2/experiment/reference.txt --out=bc2_contact_types.pdf
```

## 5. Стоимость

Типы 0–2: 85,0 + 85,1 + 110,5 с на 4 потоках (параллельно с BC6); кадры на сервере ~1 ГБ на прогон, в журнал не взяты.

## 6. Проверка

`data/c0/shot.txt`: `steps=1527`, `calc_peak=1221.11`, `amp_ratio=0.926874`; `data/c2/shot.txt`: `steps=2451`, `calc_peak=1083.32`.
