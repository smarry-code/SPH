# Воспроизведение: BC6, `bc6_shot`, тип 3

## 1. Окружение

claude-compute-8t: Ubuntu 26.04.1, ядро 7.0.0-38, Intel Xeon Cascadelake (8 ядер), gcc 15.2.0,
CMake 4.2.3, Ninja 1.13.2. OpenMP: 4 потока (`OMP_NUM_THREADS=4`, `--threads=4`); параллельно
шёл второй прогон на 4 потоках (BC2 и BC6 одновременно).

## 2. Коммит

`61e24839f4de57904f4141f97f877a11db18d683` — «cases/bc2, cases/bc6: экспериментальные постановки и утилиты bc2_shot, bc6_shot».

## 3. Diff

Нет: дерево на коммите без изменений. Архив [`repro-diff.zip`](repro-diff.zip) содержит только
`HEAD.txt` и `MANIFEST.md`.

## 4. Команды

```bash
git -c core.autocrlf=false archive --format=tar.gz -o src.tar.gz 61e2483
mkdir -p ~/ceramics/src_61e2483 && tar -xzf src.tar.gz -C ~/ceramics/src_61e2483
cd ~/ceramics/src_61e2483
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
ninja -C build -j4 bc6_shot
OMP_NUM_THREADS=4 ./build/bc6_shot --contact-type=3 --threads=4 --out=$HOME/ceramics/runs \
    --checkpoint-every=500 --overwrite
```

Запуск шёл через автозапуск сервера (`~/autostart/jobs/ceramics_bc6.job` → `~/ceramics/run_bc6.sh`:
та же команда, при наличии чекпоинта — с `--restart`). Рестарта не было.

## 5. Стоимость

95,2 с на 4 потоках; каталог прогона на сервере 790 МБ (кадры), в журнале — 0,1 МБ.

## 6. Проверка

`data/shot.txt`: `steps=1630`, `t_reached=4.30016e-06`, `calc_plateau1=182.333`, `calc_dips=1`, `calc_ratio12=0.693594`.
