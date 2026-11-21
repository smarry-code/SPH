#!/usr/bin/env python3
"""Собирает самодостаточный отчёт report_session.html из report.html: картинки fig/*.png и анимации anim/*.gif
встраиваются как data: URI, файл открывается в любом браузере без внешних файлов (шрифты Google -- с запасными).
Запуск: python scripts/build_report_file.py (из каталога эксперимента).
"""
import base64
import os
import re

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
src = open(os.path.join(HERE, "report.html"), encoding="utf-8").read()


def inline(m):
    path = os.path.join(HERE, m.group(1))
    mime = "image/gif" if path.endswith(".gif") else "image/png"
    data = base64.b64encode(open(path, "rb").read()).decode()
    return 'src="data:%s;base64,%s"' % (mime, data)


i = src.index('<div class="wrap">')  # до него -- title, шрифты, стили: в <head>
head, body = src[:i], re.sub(r'src="((?:fig|anim)/[^"]+)"', inline, src[i:])
out = ('<!doctype html>\n<html lang="ru">\n<head>\n<meta charset="utf-8">\n'
       '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
       '<style>body{margin:0}img{max-width:100%}</style>\n' + head + '</head>\n<body>\n' + body +
       '\n</body>\n</html>\n')
dst = os.path.join(HERE, "report_session.html")
open(dst, "w", encoding="utf-8").write(out)
print(dst, "%.1f МБ" % (os.path.getsize(dst) / 1e6))
