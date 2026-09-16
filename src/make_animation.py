"""
Animation of impulse-wave propagation across Lake Hawea.

Water-surface anomaly is polarity data (crest above / trough below the still level),
so it uses a diverging blue-red ramp with a neutral grey midpoint - not a rainbow.
"""
import sys, os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, LightSource, TwoSlopeNorm
import matplotlib.patheffects as pe
import imageio.v2 as imageio

DATA, OUT = "/home/user/Lake-tsunami-scenarios/data", "/home/user/Lake-tsunami-scenarios/outputs"
RES, LAKE_LEVEL, XMIN, YMAX = 30.0, 345.0, 1292000.0, 5096000.0
SUB = (80, 1520, 100, 990)
INK, INK2 = "#0b0b0b", "#52514e"
DIV = LinearSegmentedColormap.from_list(
    "div", ["#0d366b", "#184f95", "#3987e5", "#9ec5f4", "#f0efec",
            "#f9bfa4", "#eb6834", "#c94f1f", "#782d0f"])
PLACES = {"Lake Hāwea township": (1302581, 5053435)}
LABELS = {"S1_moderate": "S1 · moderate · 5×10⁵ m³",
          "S2_large": "S2 · large · 5×10⁶ m³",
          "S3_extreme": "S3 · extreme · 2×10⁷ m³"}


def to_rc(x, y, step=1):
    return ((YMAX - y) / RES - SUB[0]) / step, ((x - XMIN) / RES - SUB[2]) / step


def build(name, fps=20):
    meta = json.load(open(f"{OUT}/{name}_meta.json"))
    snaps = np.load(f"{OUT}/{name}_snaps.npy")
    times = meta["times"]
    dem = np.load(f"{DATA}/dem.npy")[SUB[0]:SUB[1], SUB[2]:SUB[3]][::2, ::2]
    lake = np.load(f"{DATA}/lake_mask.npy")[SUB[0]:SUB[1], SUB[2]:SUB[3]][::2, ::2]
    ls = LightSource(azdeg=315, altdeg=45)
    hs = ls.hillshade(dem, vert_exag=1.5, dx=RES * 2, dy=RES * 2)

    vmax = float(np.nanpercentile(np.abs(snaps), 99.7))
    vmax = max(vmax, 1.0)
    print(f"{name}: {len(snaps)} frames, colour scale ±{vmax:.1f} m")

    ny, nx = lake.shape
    fig, ax = plt.subplots(figsize=(6.0, 9.0), dpi=110)
    fig.subplots_adjust(left=.02, right=.86, top=.93, bottom=.02)
    ax.imshow(hs, cmap="gray", vmin=0, vmax=1.4, origin="upper", interpolation="bilinear")
    im = ax.imshow(snaps[0], cmap=DIV, norm=Normalize(-vmax, vmax), origin="upper",
                   interpolation="bilinear")
    r, c = to_rc(*meta["impact_xy"], step=2)
    ax.plot(c, r, marker="*", ms=15, mfc="#e34948", mec="white", mew=1.1, zorder=6)
    for nm, (x, y) in PLACES.items():
        rr, cc = to_rc(x, y, step=2)
        ax.plot(cc, rr, "o", ms=5, mfc="white", mec=INK, mew=1.3, zorder=6)
        ax.annotate(nm, (cc, rr), fontsize=8, color=INK, xytext=(7, 4),
                    textcoords="offset points", zorder=6,
                    path_effects=[pe.withStroke(linewidth=2.5, foreground="white")])
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_color("#c9c8c3")
    ttl = ax.set_title(f"Lake Hāwea · {LABELS[name]}", fontsize=11.5, color=INK, pad=10)
    clock = ax.text(.03, .022, "", transform=ax.transAxes, fontsize=13, color=INK,
                    va="bottom", family="monospace",
                    path_effects=[pe.withStroke(linewidth=3, foreground="white")])
    note = ax.text(.03, .075, "", transform=ax.transAxes, fontsize=8.5, color=INK2,
                   va="bottom",
                   path_effects=[pe.withStroke(linewidth=2.5, foreground="white")])
    cb = fig.colorbar(im, ax=ax, fraction=.038, pad=.02)
    cb.set_label("water surface relative to normal lake level (m)", fontsize=8.5)
    cb.ax.tick_params(labelsize=8, colors=INK2)

    frames = []
    fig.canvas.draw()
    for k, (fr, t) in enumerate(zip(snaps, times)):
        im.set_data(fr)
        clock.set_text(f"t = {int(t)//60:02d}:{int(t)%60:02d}")
        peak = np.nanmax(fr)
        note.set_text(f"peak crest now: {peak:5.1f} m")
        fig.canvas.draw()
        buf = np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()
        frames.append(buf)
    plt.close(fig)

    mp4 = f"{OUT}/{name}_animation.mp4"
    imageio.mimwrite(mp4, frames, fps=fps, codec="libx264", quality=8,
                     macro_block_size=1, output_params=["-pix_fmt", "yuv420p"])
    print(f"wrote {mp4}  ({os.path.getsize(mp4)/1e6:.1f} MB)")
    gif = f"{OUT}/{name}_animation.gif"
    imageio.mimwrite(gif, frames[::2], fps=max(fps // 2, 8), loop=0)
    print(f"wrote {gif}  ({os.path.getsize(gif)/1e6:.1f} MB)")
    return mp4


if __name__ == "__main__":
    for n in (sys.argv[1:] or ["S2_large"]):
        build(n)
