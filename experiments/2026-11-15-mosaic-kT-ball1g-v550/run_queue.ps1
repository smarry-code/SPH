# Запуск очереди queue.py отдельным процессом (не зависит от терминала и сессии Claude Code).
# Повторный запуск продолжает очередь с места остановки (готовые выстрелы не повторяются, незаконченный -- из чекпоинта).
# Использование: powershell -File experiments/2026-11-15-mosaic-kT-ball1g-v550/run_queue.ps1
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$p = Start-Process -FilePath "python" -ArgumentList "-u", "`"$here\queue.py`"" -WorkingDirectory $here `
    -WindowStyle Hidden -RedirectStandardOutput "$here\queue.out" -RedirectStandardError "$here\queue.err" -PassThru
Set-Content -Path "$here\queue.pid" -Value $p.Id
"queue.py PID $($p.Id)"
