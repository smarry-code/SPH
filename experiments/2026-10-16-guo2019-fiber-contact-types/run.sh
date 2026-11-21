#!/usr/bin/env bash
# Матрица шага 2.5 (notes/fiber-contact/plan-fiber-contact-types.md). Запуск из корня репозитория
# после сборки пресетом clang-windows: bash experiments/2026-10-16-guo2019-fiber-contact-types/run.sh <каталог_вывода>
R=${1:?каталог вывода}; X=build/guo2019_shot.exe; mkdir -p "$R"
for cfg in A 1T 1T_A A_1T; do
 for v in base fs1 fs2 fs3 ff1 ff2 ff3; do
  fl=""
  case $v in fs?) fl="--fs-contact=${v#fs}";; ff?) fl="--ff-contact=${v#ff}";; esac
  $X --config=$cfg --level=smoke --v0=150 --t-max=2e-5 --frames=off --threads=4 --out="$R/${cfg}_$v" $fl > "$R/${cfg}_$v.log" 2>&1
  echo "${cfg}_$v rc=$?"
 done
done
