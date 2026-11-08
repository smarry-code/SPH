#!/usr/bin/env python3
"""Графики (SVG) и сводные числа каталога по data/ (см. thin_data.py).

  v_z.svg     -- v_z шара: вся история и крупно 1,3..2,2 мс (окно пробития)
  breaks.svg  -- число порванных нитей и реакция заделки ткани Ff_z, 350..560 мкс
  energy.svg  -- куда ушла энергия шара (energy_log.csv) для 69 и 70 м/с, доли E0
  data/metrics.txt -- числа для README

Usage: python plot_series.py <каталог эксперимента>
"""
import csv, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["svg.fonttype"] = "none"   # текст в SVG -- текстом (меньше файл)
plt.rcParams["font.size"] = 10
plt.rcParams["svg.hashsalt"] = "guo1t-mu027"  # SVG побайтово повторяем (без случайных id)
EXP = sys.argv[1]
D = os.path.join(EXP, "data")
VS = [69, 70, 71, 82]
COL = {69: "#2a78d6", 70: "#d6452a", 71: "#d69a2a", 82: "#555555"}
M_BALL = 0.007  # кг, шар 12 мм (shot.log)


def kv(path):
    with open(path, encoding="utf-8") as f:
        return dict(t.split("=", 1) for t in f.read().split() if "=" in t)


def table(path):
    with open(path, newline="") as f:
        r = list(csv.DictReader(f))
    return {k: [float(q[k]) for q in r] for k in r[0]}


shot = {v: kv(os.path.join(D, "shot_v%d_long.txt" % v)) for v in VS}
mon = {v: table(os.path.join(D, "monitor_v%d.csv" % v)) for v in VS}
en = {v: table(os.path.join(D, "energy_v%d.csv" % v)) for v in VS}


def label(v):
    s = shot[v]
    res = "не пробит" if s["outcome"] == "not_perforated" else "пробит, V_R = %.2f м/с" % float(s["V_R"])
    return "%d м/с: %s" % (v, res)


# --- v_z
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for v in VS:
    t = [x * 1e6 for x in mon[v]["t"]]
    for a in ax:
        a.plot(t, mon[v]["v_z"], color=COL[v], lw=1.3, label=label(v))
for a in ax:
    a.axhline(0, color="k", lw=0.6)
    a.axvline(1500, color="k", lw=0.6, ls=":")
    a.grid(alpha=0.3)
    a.set_xlabel("t, мкс")
ax[0].set(title="v_z шара", ylabel="м/с", xlim=(0, 2200))
ax[1].set(title="крупно: после 1,3 мс (пунктир -- t_max базового прогона)", xlim=(1300, 2200), ylim=(-10, 1))
ax[0].legend(fontsize=8, loc="lower right")
fig.tight_layout()
fig.savefig(os.path.join(EXP, "v_z.svg"), metadata={"Date": None})
plt.close(fig)

# --- разрывы и реакция заделки
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for v in VS:
    t = [x * 1e6 for x in mon[v]["t"]]
    nb = [a + b for a, b in zip(mon[v]["n_broken_v"], mon[v]["n_broken_h"])]
    ax[0].step(t, nb, where="post", color=COL[v], lw=1.3, label="%d м/с" % v)
    ax[1].plot(t, [abs(x) for x in mon[v]["Ff_z"]], color=COL[v], lw=1.1, label="%d м/с" % v)
ax[0].axhline(3, color="k", ls="--", lw=0.8, label="статья: 3 / 0 при 73 м/с")
ax[0].set(title="порвано нитей (верт. + гориз.)", xlabel="t, мкс", ylabel="нитей", xlim=(350, 560))
ax[1].set(title="|Ff_z| -- реакция заделки ткани", xlabel="t, мкс", ylabel="Н", xlim=(0, 2200))
for a in ax:
    a.grid(alpha=0.3)
    a.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(EXP, "breaks.svg"), metadata={"Date": None})
plt.close(fig)

# --- энергия: доли E0 для 69 и 70
CH = [("KE_solid", "KE шара"), ("D_friction", "трение"), ("D_contact_fiber", "контакт нитей (закон 3)"),
      ("KE_fiber", "KE ткани"), ("U_bond", "упругая энергия связей"), ("W_p_fiber", "пластичность нитей"),
      ("D_xsph", "XSPH"), ("residual", "невязка")]
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
for a, v in zip(ax, (69, 70)):
    e = en[v]
    e0 = 0.5 * M_BALL * v * v
    t = [x * 1e6 for x in e["time"]]
    for k, lab in CH:
        a.plot(t, [x / e0 for x in e[k]], lw=1.3, label=lab)
    a.set(title="%s (E0 = %.2f Дж)" % (label(v), e0), xlabel="t, мкс", xlim=(0, 2200))
    a.grid(alpha=0.3)
ax[0].set_ylabel("доля E0")
ax[1].legend(fontsize=8, loc="center right")
fig.tight_layout()
fig.savefig(os.path.join(EXP, "energy.svg"), metadata={"Date": None})
plt.close(fig)

# --- числа
out = []
for v in VS:
    s, m, e = shot[v], mon[v], en[v]
    e0 = 0.5 * M_BALL * v * v
    t, vz = m["t"], m["v_z"]
    i15 = min(range(len(t)), key=lambda i: abs(t[i] - 1.5e-3))
    cross = next((t[i] for i in range(1, len(t)) if vz[i - 1] < 0 <= vz[i]), None)
    imax = max(range(len(t)), key=lambda i: abs(m["Ff_z"][i]))
    out.append("V0=%d E0=%.4f Дж outcome=%s V_R=%s t_end=%.1f мкс steps_long=%s" % (
        v, e0, s["outcome"], s["V_R"], float(s["t_end"]) * 1e6, s["steps"]))
    out.append("  v_z(1,5 мс)=%.3f м/с (строка t=%.1f мкс); v_z(t_end)=%.3f; z_c(t_end)=%.2f мм; v_z=0 при %s" % (
        vz[i15], t[i15] * 1e6, vz[-1], m["z_c"][-1] * 1e3, "%.0f мкс" % (cross * 1e6) if cross else "нет"))
    out.append("  первый разрыв %.1f мкс, w_f при нём %.2f мм; порвано %s/%s; |Ff_z|max(выборка)=%.0f Н при %.0f мкс; "
               "F_clamp_f_z_max=%s Н" % (float(s["t_first_break"]) * 1e6, float(s["w_f_at_first_break"]) * 1e3,
                                         s["n_broken_v"], s["n_broken_h"], abs(m["Ff_z"][imax]), t[imax] * 1e6,
                                         s["F_clamp_f_z_max"]))
    out.append("  энергия на t_end, доли E0: " + ", ".join("%s %.4f" % (k, e[k][-1] / e0) for k, _ in CH))
    out.append("  e (поглощено шаром) = %s Дж = %.4f E0; E_tot_end/E_tot_0 = %s; max δff/d = %s, max δfs/d = %s, "
               "n_tunnel = %s" % (s["e"], float(s["e"]) / e0, s["E_tot_end/E_tot_0"], s["max_delta_ff/d"],
                                  s["max_delta_fs/d"], s["n_tunnel"]))
with open(os.path.join(D, "metrics.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(out) + "\n")
print("\n".join(out))
