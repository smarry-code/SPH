"""Разлёт бусин: по кадрам fiber_output -- max|v|, число бусин с |v| > порогов, доля провисших связей
(l < l0), число бусин без связей; для pinball / Swegle / закона 4 (1T, mu = 0.1, V0 = 70)."""
import glob, math, os, sys


def frame_stats(path):
    t = None
    speeds = []
    deg = {}
    slack = taut = 0
    mode = None
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("TIME"):
                t = float(line.split()[1])
                continue
            if line.startswith("PARTICLES"):
                mode = "p"
                continue
            if line.startswith("BONDS"):
                mode = "b"
                continue
            if not line or line[0] == "#":
                continue
            s = line.split()
            if mode == "p":
                vx, vy, vz = float(s[6]), float(s[7]), float(s[8])
                speeds.append(math.sqrt(vx * vx + vy * vy + vz * vz))
            elif mode == "b":
                i, j = int(s[0]), int(s[1])
                deg[i] = deg.get(i, 0) + 1
                deg[j] = deg.get(j, 0) + 1
                if float(s[3]) < float(s[2]):
                    slack += 1
                else:
                    taut += 1
    n_free = sum(1 for k in range(len(speeds)) if deg.get(k, 0) == 0)
    return t, speeds, slack, taut, n_free


for label, root in zip(sys.argv[1::2], sys.argv[2::2]):
    frames = sorted(glob.glob(os.path.join(root, "fiber_output_*.txt")))
    print("==", label, len(frames), "кадров")
    print("  t, мкс | max|v| | |v|>100 | |v|>200 | |v|>300 | провисших связей, % | бусин без связей")
    for fp in frames[::4] + [frames[-1]]:
        t, sp, slack, taut, nfree = frame_stats(fp)
        print("  %6.0f | %6.0f | %6d | %6d | %6d | %5.1f | %d" % (
            t * 1e6, max(sp), sum(v > 100 for v in sp), sum(v > 200 for v in sp), sum(v > 300 for v in sp),
            100.0 * slack / max(slack + taut, 1), nfree))
