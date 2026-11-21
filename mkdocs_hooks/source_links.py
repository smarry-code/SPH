# source_links.py
# MkDocs hook: публикация исходников проекта (src/, include/) в сайт и
# нормализация суффикса :N в ссылках на исходный код.
#
# Зачем это нужно:
#   Ссылки вида [`name()`](../../../src/linked_list.cpp:32) в разделах архитектуры
#   указывают на файлы ВНЕ docs_dir, поэтому в собранном сайте MkDocs они дают 404,
#   а суффикс :32 (номер строки) не является корректным веб-URL (браузер ищет файл
#   с двоеточием в имени). Hook решает обе проблемы:
#
#     1. on_files     — монтирует .cpp/.hpp/.h из src/ и include/ (со всеми подкаталогами) как виртуальные
#                        файлы сайта (site/src/..., site/include/...), чтобы
#                        относительные пути ../../../src/... резолвились;
#     2. on_post_page — заменяет в href суффикс :N на #LN, чтобы файл открывался
#                        в браузере, а не возвращал 404.
#
# Подход «Docs-as-Code»: при этом исходный текст ссылок (../../../src/x.cpp:N)
# остаётся неизменным и по-прежнему удобен для навигации в IDE/GitHub (VS Code
# открывает файл на нужной строке по суффиксу :N).

import os
import re
from glob import glob

from mkdocs.structure.files import File

# Расширения исходников, которые публикуются в сайт.
_SOURCE_EXTS = ("*.cpp", "*.hpp", "*.h", "*.cc")

# Каталоги исходников в корне проекта, монтируемые в сайт под тем же именем.
_SOURCE_DIRS = ("src", "include")


def on_files(files, config):
    """Монтирует src/ и include/ как виртуальные файлы сайта MkDocs."""
    docs_dir = config["docs_dir"]
    project_root = os.path.abspath(os.path.join(docs_dir, ".."))
    site_dir = config["site_dir"]
    use_directory_urls = config["use_directory_urls"]

    for sub in _SOURCE_DIRS:
        src_abs = os.path.join(project_root, sub)
        if not os.path.isdir(src_abs):
            continue
        for pattern in _SOURCE_EXTS:
            for abs_path in glob(os.path.join(src_abs, "**", pattern), recursive=True):
                # Относительный путь внутри сайта (например "src/linked_list.cpp").
                rel = os.path.relpath(abs_path, project_root).replace(os.sep, "/")
                files.append(
                    File(
                        path=rel,
                        src_dir=project_root,
                        dest_dir=site_dir,
                        use_directory_urls=use_directory_urls,
                    )
                )
    return files


# Шаблон: href="...src/linked_list.cpp:32"  →  href="...src/linked_list.cpp#L32"
_HREF_RE = re.compile(r'(href=")([^"]+?\.(?:cpp|hpp|h|cc)):(\d+)(")')


def on_post_page(output, **kwargs):
    """Переписывает суффикс :N → #LN в href на исходный код (для браузера)."""
    return _HREF_RE.sub(
        lambda m: m.group(1) + m.group(2) + "#L" + m.group(3) + m.group(4),
        output,
    )
