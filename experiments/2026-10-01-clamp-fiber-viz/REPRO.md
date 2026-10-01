# Воспроизведение: закрепление нитей, EF1 / EF1-b / EF2 / EF4 / EF5

## 1. Окружение

| | |
|---|---|
| ОС | macOS (Darwin 24.6.0), Apple Silicon |
| Компилятор | Apple clang version 17.0.0 (clang-1700.0.13.5), `/usr/bin/c++` |
| Сборка | CMake 4.4.3, Unix Makefiles, `-DCMAKE_BUILD_TYPE=Release` |
| Потоки | `OMP_NUM_THREADS=1` |
| Python (анимации, графики) | 3.14.7, numpy 2.5.3, matplotlib 3.11.2, Pillow 12.3.0 (venv; в системном Python этих пакетов нет) |

## 2. Коммит

```
8ae109beaa7d34380efe163597c1f65f170da91b
clamp e2e: сцены Solid и Fiber (ES1-ES3, EF1-EF6, EK1), правки по ревью Codex
ветка clamp
```

## 3. Отличия рабочего дерева от этого коммита

Патчей нет: прогоны сняты на чистом дереве HEAD 8ae109b. [`repro-diff.zip`](repro-diff.zip) содержит `HEAD.txt` (хэш, автор, дата, заголовок) и `MANIFEST.md` с этим утверждением.

| Файл архива | Содержание |
|---|---|
| `HEAD.txt` | полный хэш и заголовок коммита |
| `MANIFEST.md` | «патчей нет», условие чистоты дерева |

Проверка чистоты: `git status --short -- include src tests tools CMakeLists.txt` пуст (проверено 2026-10-01 перед и после прогонов). Не закоммичены только `Clamp/`, `build_clamp/`, `fiber/*.md`, `experiments/` — на сборку и расчёт не влияют. Наложение не требуется; достаточно `git checkout 8ae109beaa7d34380efe163597c1f65f170da91b` (или `git worktree add --detach /tmp/wt <хэш>`).

## 4. Команды

```bash
# 1. сборка (в корне репозитория на HEAD 8ae109b; у автора каталог сборки — build_clamp)
cmake -S . -B build_clamp -DCMAKE_BUILD_TYPE=Release
cmake --build build_clamp --target test_e2e_clamp -j8

# 2. прогоны с кадрами; запускать НЕ из корня репозитория
mkdir -p /tmp/run && cd /tmp/run
export OMP_NUM_THREADS=1
B=<repo>/build_clamp/test_e2e_clamp
$B --case=EF1               --solver=pc --frames-dir=fr --frames=80 | tee EF1.log
$B --case=EF1 --variant=b   --solver=pc --frames-dir=fr --frames=80 | tee EF1-b.log
$B --case=EF2               --solver=pc --frames-dir=fr --frames=80 | tee EF2.log
$B --case=EF4               --solver=pc --frames-dir=fr --frames=80 | tee EF4.log   # сам делает два прогона: mu = 0 и mu > 0
$B --case=EF5               --solver=pc --frames-dir=fr --frames=80 | tee EF5.log

# 3. анимации и снимки (venv с numpy/matplotlib/pillow), из корня репозитория
python tools/clamp_animate.py --all /tmp/run/fr --out-dir /tmp/run/anim --max-frames 60 --dpi 80 --snapshot 15,40

# 4. графики и metrics.txt из data/*/reaction.csv
python experiments/2026-10-01-clamp-fiber-viz/scripts/plot_reaction.py
```

В `data/` лежат `reaction.csv` и `clamp_mask.txt` из `fr/<прогон>/` (по подкаталогу на прогон) и логи в `data/logs/`.

## 5. Стоимость

Каждый запуск теста — 0–1 с стенного времени на одном потоке (260–1462 шагов). Анимации всех шести прогонов (60 кадров, dpi 80) — около 1 мин суммарно. Исходные кадры 80 на прогон не хранятся. Каталог — 5,6 МБ (GIF 4,3 МБ, PNG 1,0 МБ, SVG 0,07 МБ, данные 0,24 МБ, архив 1 КБ); ориентир 1 МБ поднят владельцем до 10 МБ.

## 6. Проверка повтора

Проверено 2026-10-01 на HEAD 8ae109b: `git status --short -- include src tests tools CMakeLists.txt` пуст. Прогоны сняты и повторены (EF1-b, EF2, EF4 в отдельном каталоге): `reaction.csv` и `clamp_mask.txt` совпали побайтово для четырёх проверенных повторов (EF1-b, EF2, EF4 mu и mu0), а также побайтово совпали со снятыми в прошлой версии каталога (HEAD 2b7f1a3 + три патча, другая сборка) для всех шести прогонов — расчёт детерминирован при `OMP_NUM_THREADS=1`. Логи совпали с прежними с точностью до пути кадров, стенного времени и нового WARN в EF4 (пар «зажатая-зажатая» в сцене нет).

Должно совпасть: EF1 — плато −0,196348 Н (отклонение −7,609·10⁻⁴ %), падение ниже 10 % плато на 7,5958 мс (+1,277 %), импульс 0,462 %, энергия 0,652 %; EF1-b — разрыв связи (2,3) на шаге 2, R_x шага 2 = −0,27207 Н, затем 0; EF2 — T = 2,573214·10⁻² с (−1,87·10⁻³ %), среднее R_x ∓0,117810 Н, амплитуда R_y 3,6860·10⁻⁴ Н; EF4 — J_n = 1,2276·10⁻⁶ Н·с, J_f = 6,1360·10⁻⁷ Н·с, J_f/(mu·J_n) = 0,9997; EF5 — R_x ≡ 0 (0 ненулевых шагов из 260). Все пять запусков заканчиваются `PASS (провалов 0)`. Побитовое совпадение логов не обещается на другом компиляторе или числе потоков.
