# Ждёт завершения текущего процесса очереди (PID из queue.pid) и запускает очередь снова -- она выполнит EXTRA (17T).
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$qpid = [int](Get-Content "$here\queue.pid")
try { Wait-Process -Id $qpid -ErrorAction Stop } catch {}
& powershell -ExecutionPolicy Bypass -File "$here\run_queue.ps1"
