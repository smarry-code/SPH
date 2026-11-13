#!/usr/bin/env python3
"""Вызывающие memmove/memset/memcpy (libc) в dwarf-профиле: по образцам главного потока, ключ -- первые кадры
проекта над libc (без inline-обёрток std::). Запуск: perf script -F tid,ip,sym,dso | python3 libc_callers.py <главный tid>"""
import collections
import re
import sys

main_tid = sys.argv[1]
SKIP = re.compile(r"^(__|std::|_M_|assign|copy|fill|operator\(\)|_Vector|vector|allocator|construct|resize|"
                  r"__copy|__fill|__uninit|uninitialized|_S_|_M_)")
cnt = collections.Counter(); tot = collections.Counter(); nmain = 0; nall = 0
tid, frames = None, []

def flush():
    global nmain, nall
    if tid is None:
        return
    nall += 1
    if tid != main_tid:
        return
    nmain += 1
    if not frames or "libc.so" not in frames[0][1]:
        return
    top = re.sub(r"_(evex|avx2|avx512|sse2)\w*", "", frames[0][0])
    chain = []
    par = any("libgomp" in d for _, d in frames)
    reg = None
    for s, d in frames[1:]:
        if reg is None and "parallel_for<" in s:
            m = re.search(r"parallel_for<([\w:]+(?:<[\w:]+)?)", s)
            reg = m.group(1) if m else None
    for s, d in frames[1:]:
        s = s.replace(" (inlined)", "")
        n = re.split(r"[<(]", s)[0].split("::")[-1]
        if not n or SKIP.match(n):
            continue
        if n not in chain:
            chain.append(n)
        if len(chain) == 3:
            break
    tot[top] += 1
    where = ("регион " + (reg or "?")) if par else "вне региона"
    cnt[(top, where, " < ".join(chain) or "?")] += 1

for line in sys.stdin:
    if not line.strip():
        flush(); tid, frames = None, []; continue
    if not line[0].isspace():  # заголовок образца (tid [время]...) -- без ведущих пробелов; кадры -- с пробелами
        flush(); tid, frames = line.split()[0], []; continue
    m = re.match(r"^\s+[0-9a-f]+\s+(.*?)\s+\((.*)\)\s*$", line)
    if m:
        frames.append((m.group(1), m.group(2)))
    else:
        m = re.match(r"^\s+[0-9a-f]+\s+(.*)$", line)
        if m:
            frames.append((m.group(1), ""))
flush()
if nmain == 0:
    print(f"образцов всего {nall}, образцов главного потока (tid {main_tid}) нет -- нечего разбирать")
    sys.exit(0)
print(f"образцов всего {nall}, главного потока {nmain}")
for k, v in tot.most_common():
    print(f"{k}: {v} образцов = {100 * v / nmain:.2f} % главного потока")
print("\nвызывающие (доля образцов главного потока):")
for (top, where, chain), v in cnt.most_common(30):
    print(f"  {100 * v / nmain:5.2f} %  {top:10s} {where:70s} {chain}")
