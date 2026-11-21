# Валидация кривых развёртки: монотонность по времени, длина, конечное время, пик.
import glob, os
D = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
for f in sorted(glob.glob(os.path.join(D, "curve_*.txt"))):
    T = []; V = []
    for L in open(f, encoding="utf-8", errors="replace"):
        L = L.strip()
        if not L or L[0] == "#": continue
        a = L.split(); T.append(float(a[0])); V.append(float(a[1]))
    mono = all(T[i] < T[i + 1] for i in range(len(T) - 1))
    print("%-22s tochek=%5d mono=%-5s t_end=%.4e peak=%7.1f" % (
        os.path.basename(f), len(T), mono, T[-1], max(V)))
