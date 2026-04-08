#!/usr/bin/env python3

import argparse
import mpmath as mp
import matplotlib.pyplot as plt

mp.mp.dps = 80


# --------------------------------------------------
# User inputs
# --------------------------------------------------
START_ZERO_INDEX = 1      # first nontrivial zeta zero
NUM_CHAMBERS = 5
NUM_LEVELS = 5
CRAWL_STEPS = 2000
ROOT_TOL = mp.mpf("1e-30")
T_SPACE_SAMPLES = 1200
ENERGY_SPACE_SAMPLES = 600
SABOTAGE_WINDOW_RADIUS = mp.mpf("0.35")
SABOTAGE_WINDOW_SAMPLES = 700
DEFAULT_SABOTAGE_DELTA = mp.mpf("0.10")


# --------------------------------------------------
# Core function
# --------------------------------------------------
def energy(t):
    return abs(mp.zeta(mp.mpf("0.5") + 1j * t))


# --------------------------------------------------
# Helpers
# --------------------------------------------------
def get_zero_t(n):
    z = mp.zetazero(n)
    return mp.im(z)


def get_consecutive_zeros(start_index, num_chambers):
    # need num_chambers + 1 zeros for num_chambers intervals
    return [get_zero_t(start_index + i) for i in range(num_chambers + 1)]


# --------------------------------------------------
# Crawl to get a good antinode guess
# --------------------------------------------------
def crawl_antinode(a, b, steps=2000):
    width = b - a
    best_t = a
    best_e = energy(a)

    for i in range(steps + 1):
        t = a + width * mp.mpf(i) / mp.mpf(steps)
        e = energy(t)
        if e > best_e:
            best_e = e
            best_t = t

    return best_t, best_e


# --------------------------------------------------
# Refine antinode by solving dE/dt = 0
# --------------------------------------------------
def refine_antinode(t_guess):
    fprime = lambda t: mp.diff(energy, t)
    try:
        t_star = mp.findroot(fprime, t_guess)
    except Exception:
        t_star = t_guess
    return t_star, energy(t_star)


# --------------------------------------------------
# Find sign-change bracket for E(t) - target on [a,b]
# --------------------------------------------------
def bracket_level(target, a, b, steps=500):
    f = lambda t: energy(t) - target

    prev_t = a
    prev_v = f(a)

    for i in range(1, steps + 1):
        t = a + (b - a) * mp.mpf(i) / mp.mpf(steps)
        v = f(t)

        if prev_v == 0:
            return prev_t, prev_t
        if v == 0:
            return t, t
        if prev_v * v < 0:
            return prev_t, t

        prev_t = t
        prev_v = v

    raise RuntimeError("No bracket found for target energy {}".format(target))


# --------------------------------------------------
# Simple bisection
# --------------------------------------------------
def bisect_root(f, a, b, tol=mp.mpf("1e-30"), max_iter=300):
    fa = f(a)
    fb = f(b)

    if fa == 0:
        return a
    if fb == 0:
        return b
    if fa * fb > 0:
        raise ValueError("Bisection requires opposite signs")

    for _ in range(max_iter):
        m = (a + b) / 2
        fm = f(m)

        if fm == 0 or abs(b - a) < tol:
            return m

        if fa * fm < 0:
            b = m
            fb = fm
        else:
            a = m
            fa = fm

    return (a + b) / 2


# --------------------------------------------------
# Chamber analysis
# --------------------------------------------------
def analyze_chamber(left_zero, right_zero, num_levels=5):
    t_guess, _ = crawl_antinode(left_zero, right_zero, steps=CRAWL_STEPS)
    antinode_t, antinode_e = refine_antinode(t_guess)

    levels = [
        antinode_e * mp.mpf(k) / mp.mpf(num_levels + 1)
        for k in range(1, num_levels + 1)
    ]

    rows = []

    for idx, target in enumerate(levels, start=1):
        f = lambda t: energy(t) - target

        la, lb = bracket_level(target, left_zero, antinode_t)
        t_left = bisect_root(f, la, lb, tol=ROOT_TOL)

        ra, rb = bracket_level(target, antinode_t, right_zero)
        t_right = bisect_root(f, ra, rb, tol=ROOT_TOL)

        rows.append({
            "index": idx,
            "target_energy": target,
            "normalized_energy": target / antinode_e,
            "t_left": t_left,
            "t_right": t_right,
        })

    return {
        "left_zero": left_zero,
        "right_zero": right_zero,
        "antinode_t": antinode_t,
        "antinode_e": antinode_e,
        "rows": rows,
    }


# --------------------------------------------------
# Build raw T-space samples
# --------------------------------------------------
def build_t_space_samples(left_zero, right_zero, antinode_e, sample_count=T_SPACE_SAMPLES):
    xs = []
    ys = []
    for i in range(sample_count + 1):
        t = left_zero + (right_zero - left_zero) * mp.mpf(i) / mp.mpf(sample_count)
        e = energy(t)
        xs.append(float(t))
        ys.append(float(e / antinode_e))
    return xs, ys


# --------------------------------------------------
# Build energy-space samples for one chamber
# --------------------------------------------------
def build_energy_space_samples(left_zero, antinode_t, right_zero, antinode_e, sample_count=ENERGY_SPACE_SAMPLES):
    left_x = []
    left_y = []
    left_t = []

    for i in range(sample_count + 1):
        t = left_zero + (antinode_t - left_zero) * mp.mpf(i) / mp.mpf(sample_count)
        u = energy(t) / antinode_e
        x = u - 1
        left_x.append(float(x))
        left_y.append(float(u))
        left_t.append(t)

    right_x = []
    right_y = []
    right_t = []

    for i in range(sample_count + 1):
        t = antinode_t + (right_zero - antinode_t) * mp.mpf(i) / mp.mpf(sample_count)
        u = energy(t) / antinode_e
        x = 1 - u
        right_x.append(float(x))
        right_y.append(float(u))
        right_t.append(t)

    return left_x, left_y, left_t, right_x, right_y, right_t


def build_forced_energy_space_samples(left_zero, antinode_t, right_zero, antinode_e, correction, sample_count=ENERGY_SPACE_SAMPLES):
    left_x = []
    left_y = []
    right_x = []
    right_y = []

    correction_norm = correction / antinode_e

    for i in range(sample_count + 1):
        t = left_zero + (antinode_t - left_zero) * mp.mpf(i) / mp.mpf(sample_count)
        u = (energy(t) - correction) / antinode_e
        x = u - 1
        left_x.append(float(x))
        left_y.append(float(u))

    for i in range(sample_count + 1):
        t = antinode_t + (right_zero - antinode_t) * mp.mpf(i) / mp.mpf(sample_count)
        u = (energy(t) - correction) / antinode_e
        x = 1 - u
        right_x.append(float(x))
        right_y.append(float(u))

    return left_x, left_y, right_x, right_y, correction_norm


# --------------------------------------------------
# Sabotage model
# --------------------------------------------------
def build_sabotage_event(true_zero_index, delta_t, window_radius=SABOTAGE_WINDOW_RADIUS, sample_count=SABOTAGE_WINDOW_SAMPLES):
    true_zero_t = get_zero_t(true_zero_index)
    bad_zero_t = true_zero_t + delta_t

    real_zero_t = get_zero_t(true_zero_index)
    real_energy_at_true_zero = energy(real_zero_t)
    energy_at_bad_zero = energy(bad_zero_t)
    correction = energy_at_bad_zero - real_energy_at_true_zero

    left = bad_zero_t - window_radius
    right = bad_zero_t + window_radius

    xs = []
    real_ys = []
    fixed_ys = []

    for i in range(sample_count + 1):
        t = left + (right - left) * mp.mpf(i) / mp.mpf(sample_count)
        e = energy(t)
        xs.append(float(t))
        real_ys.append(float(e))
        fixed_ys.append(float(e - correction))

    return {
        "true_zero_index": true_zero_index,
        "true_zero_t": true_zero_t,
        "bad_zero_t": bad_zero_t,
        "real_zero_t": real_zero_t,
        "real_energy_at_true_zero": real_energy_at_true_zero,
        "energy_at_bad_zero": energy_at_bad_zero,
        "correction": correction,
        "window_left": left,
        "window_right": right,
        "xs": xs,
        "real_ys": real_ys,
        "fixed_ys": fixed_ys,
    }


# --------------------------------------------------
# Plot all chambers, plus sabotage if requested
# --------------------------------------------------
def make_plot(chambers, sabotage_event=None):
    if sabotage_event is None:
        fig = plt.figure(figsize=(16, 12))
        ax1 = fig.add_subplot(2, 1, 1)
        ax2 = fig.add_subplot(2, 1, 2)
    else:
        fig = plt.figure(figsize=(16, 16))
        ax1 = fig.add_subplot(3, 1, 1)
        ax2 = fig.add_subplot(3, 1, 2)
        ax3 = fig.add_subplot(3, 1, 3)

    # ----------------------------------------------
    # Plot 1: raw chambers in T space
    # ----------------------------------------------
    for i, chamber in enumerate(chambers, start=1):
        left_zero = chamber["left_zero"]
        right_zero = chamber["right_zero"]
        antinode_t = chamber["antinode_t"]
        antinode_e = chamber["antinode_e"]
        rows = chamber["rows"]

        xs, ys = build_t_space_samples(left_zero, right_zero, antinode_e)
        ax1.plot(xs, ys, linewidth=2)

        ax1.axvline(float(left_zero), linestyle="--", linewidth=0.8)
        ax1.axvline(float(antinode_t), linestyle="--", linewidth=0.8)
        ax1.axvline(float(right_zero), linestyle="--", linewidth=0.8)

        ax1.text(float(antinode_t), 1.01, f"C{i}", ha="center", va="bottom", fontsize=9)

        for row in rows:
            e = float(row["normalized_energy"])
            tl = float(row["t_left"])
            tr = float(row["t_right"])
            ax1.plot([tl, tr], [e, e], linewidth=0.8)
            ax1.plot([tl], [e], marker="o", markersize=3)
            ax1.plot([tr], [e], marker="o", markersize=3)

    if sabotage_event is not None:
        bad_zero_t = float(sabotage_event["bad_zero_t"])
        true_zero_t = float(sabotage_event["true_zero_t"])
        chamber_antinode_e = None
        for chamber in chambers:
            if abs(chamber["left_zero"] - sabotage_event["true_zero_t"]) < mp.mpf("1e-20"):
                chamber_antinode_e = chamber["antinode_e"]
                break
            if abs(chamber["right_zero"] - sabotage_event["true_zero_t"]) < mp.mpf("1e-20"):
                chamber_antinode_e = chamber["antinode_e"]
                break
        if chamber_antinode_e is None:
            chamber_antinode_e = chambers[0]["antinode_e"]

        ax1.axvline(true_zero_t, linestyle=":", linewidth=1.8)
        ax1.axvline(bad_zero_t, linestyle="-.", linewidth=1.8)
        ax1.plot([bad_zero_t], [float(sabotage_event["energy_at_bad_zero"] / chamber_antinode_e)], marker="x", markersize=8)
        ax1.text(true_zero_t, 0.08, "true zero", rotation=90, va="bottom", ha="right", fontsize=8)
        ax1.text(bad_zero_t, 0.08, "sabotaged zero", rotation=90, va="bottom", ha="left", fontsize=8)

    ax1.set_title("Five chambers in T space")
    ax1.set_xlabel("T")
    ax1.set_ylabel("Normalized energy (per chamber)")
    ax1.set_ylim(-0.02, 1.05)

    # ----------------------------------------------
    # Plot 2: all chambers in energy space
    # ----------------------------------------------
    chamber_spacing = 2.8

    for i, chamber in enumerate(chambers):
        center = i * chamber_spacing

        left_zero = chamber["left_zero"]
        right_zero = chamber["right_zero"]
        antinode_t = chamber["antinode_t"]
        antinode_e = chamber["antinode_e"]
        rows = chamber["rows"]

        left_x, left_y, left_t, right_x, right_y, right_t = build_energy_space_samples(
            left_zero, antinode_t, right_zero, antinode_e
        )

        left_x = [x + center for x in left_x]
        right_x = [x + center for x in right_x]

        ax2.plot(left_x, left_y, linewidth=2)
        ax2.plot(right_x, right_y, linewidth=2)

        ax2.axvline(center - 1.0, linestyle="--", linewidth=0.8)
        ax2.axvline(center + 0.0, linestyle="--", linewidth=0.8)
        ax2.axvline(center + 1.0, linestyle="--", linewidth=0.8)

        ax2.text(center - 1.0, -0.04, "node", ha="center", va="top", fontsize=8)
        ax2.text(center + 0.0, -0.04, "antinode", ha="center", va="top", fontsize=8)
        ax2.text(center + 1.0, -0.04, "node", ha="center", va="top", fontsize=8)
        ax2.text(center, 1.03, f"C{i+1}", ha="center", va="bottom", fontsize=9)

        for row in rows:
            u = float(row["normalized_energy"])
            x_left = center + (u - 1.0)
            x_right = center + (1.0 - u)

            ax2.plot([x_left, x_right], [u, u], linewidth=0.8)
            ax2.plot([x_left], [u], marker="o", markersize=3)
            ax2.plot([x_right], [u], marker="o", markersize=3)

            ax2.text(
                x_left - 0.03,
                u,
                mp.nstr(row["t_left"], 10),
                ha="right",
                va="center",
                fontsize=6
            )
            ax2.text(
                x_right + 0.03,
                u,
                mp.nstr(row["t_right"], 10),
                ha="left",
                va="center",
                fontsize=6
            )

    if sabotage_event is not None:
        correction = sabotage_event["correction"]
        deficit_depths = []

        for i, chamber in enumerate(chambers):
            center = i * chamber_spacing
            antinode_e = chamber["antinode_e"]
            left_zero = chamber["left_zero"]
            right_zero = chamber["right_zero"]
            antinode_t = chamber["antinode_t"]

            forced_left_x, forced_left_y, forced_right_x, forced_right_y, correction_norm = build_forced_energy_space_samples(
                left_zero, antinode_t, right_zero, antinode_e, correction
            )

            forced_left_x = [x + center for x in forced_left_x]
            forced_right_x = [x + center for x in forced_right_x]
            deficit_depths.append(float(correction_norm))

            ax2.plot(
                forced_left_x,
                forced_left_y,
                linestyle='--',
                linewidth=1.6,
                color='tab:red',
                alpha=0.95
            )
            ax2.plot(
                forced_right_x,
                forced_right_y,
                linestyle='--',
                linewidth=1.6,
                color='tab:red',
                alpha=0.95
            )

            y_def = float(-correction_norm)
            x_left_def = center - 1.0 - float(correction_norm)
            x_right_def = center + 1.0 + float(correction_norm)

            ax2.plot([center - 1.0, x_left_def], [0.0, y_def], linestyle=':', linewidth=1.0, color='tab:red')
            ax2.plot([center + 1.0, x_right_def], [0.0, y_def], linestyle=':', linewidth=1.0, color='tab:red')
            ax2.plot([x_left_def], [y_def], marker='x', markersize=6, color='tab:red')
            ax2.plot([x_right_def], [y_def], marker='x', markersize=6, color='tab:red')

        affected_index = sabotage_event["true_zero_index"] - START_ZERO_INDEX
        if 0 <= affected_index < len(chambers) - 1:
            shared_center_left = affected_index * chamber_spacing
            shared_center_right = (affected_index + 1) * chamber_spacing
            left_depth = float(correction / chambers[affected_index]["antinode_e"])
            right_depth = float(correction / chambers[affected_index + 1]["antinode_e"])
            ax2.text(
                shared_center_left + 1.02 + left_depth,
                -left_depth,
                'shared-node deficit',
                fontsize=7,
                ha='left',
                va='bottom',
                color='tab:red'
            )
            ax2.text(
                shared_center_right - 1.02 - right_depth,
                -right_depth,
                'shared-node deficit',
                fontsize=7,
                ha='right',
                va='bottom',
                color='tab:red'
            )

        deficit_floor = -max(deficit_depths) if deficit_depths else -0.02
        ax2.set_ylim(min(deficit_floor - 0.08, -0.02), 1.08)
    else:
        ax2.set_ylim(-0.02, 1.08)

    ax2.set_title("Five chambers in energy space")
    ax2.set_xlabel("Energy-space branch position (offset by chamber)")
    ax2.set_ylabel("Normalized energy")
    ax2.set_xlim(-1.4, (len(chambers) - 1) * chamber_spacing + 1.4)

    if sabotage_event is not None:
        correction = sabotage_event["correction"]
        energy_at_bad_zero = sabotage_event["energy_at_bad_zero"]

        ax3.plot(sabotage_event["xs"], sabotage_event["real_ys"], linewidth=2, label="real energy")
        ax3.plot(sabotage_event["xs"], sabotage_event["fixed_ys"], linewidth=2, label="corrected / forced profile")
        ax3.axvline(float(sabotage_event["true_zero_t"]), linestyle=":", linewidth=1.5, label="true zero")
        ax3.axvline(float(sabotage_event["bad_zero_t"]), linestyle="-.", linewidth=1.5, label="sabotaged zero")
        ax3.axhline(float(correction), linestyle="--", linewidth=1.2, label="correction energy")
        ax3.plot([float(sabotage_event["bad_zero_t"])] , [float(energy_at_bad_zero)], marker="o", markersize=6)
        ax3.plot([float(sabotage_event["bad_zero_t"])] , [0.0], marker="x", markersize=8)
        ax3.set_title("Sabotage event: forced node correction")
        ax3.set_xlabel("T")
        ax3.set_ylabel("Absolute energy")
        ax3.legend(loc="best")

    plt.tight_layout()
    plt.show()


# --------------------------------------------------
# CLI
# --------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot zeta chambers and optionally inject a sabotaged zero event."
    )
    parser.add_argument(
        "-sabotage",
        action="store_true",
        help="Inject a bad zero near a true zero and show the correction event."
    )
    parser.add_argument(
        "--sabotage-index",
        type=int,
        default=START_ZERO_INDEX,
        help="1-based zeta zero index to sabotage. Default: 1"
    )
    parser.add_argument(
        "--sabotage-delta",
        type=str,
        default=str(DEFAULT_SABOTAGE_DELTA),
        help="Additive T shift for the injected bad zero. Example: 0.10 or -0.05"
    )
    return parser.parse_args()


# --------------------------------------------------
# Main
# --------------------------------------------------
def main():
    args = parse_args()
    sabotage_delta = mp.mpf(args.sabotage_delta)

    zeros = get_consecutive_zeros(START_ZERO_INDEX, NUM_CHAMBERS)

    chambers = []
    for i in range(NUM_CHAMBERS):
        chamber = analyze_chamber(
            zeros[i],
            zeros[i + 1],
            num_levels=NUM_LEVELS
        )
        chambers.append(chamber)

    for i, chamber in enumerate(chambers, start=1):
        print("Chamber", i)
        print("Left zero     :", mp.nstr(chamber["left_zero"], 30))
        print("Antinode T    :", mp.nstr(chamber["antinode_t"], 30))
        print("Right zero    :", mp.nstr(chamber["right_zero"], 30))
        print("Antinode Emax :", mp.nstr(chamber["antinode_e"], 30))
        print()

        for row in chamber["rows"]:
            print({
                "energy_norm": mp.nstr(row["normalized_energy"], 12),
                "t_left": mp.nstr(row["t_left"], 18),
                "t_right": mp.nstr(row["t_right"], 18),
            })
        print()

    sabotage_event = None
    if args.sabotage:
        sabotage_event = build_sabotage_event(args.sabotage_index, sabotage_delta)
        print("SABOTAGE EVENT")
        print("True zero index           :", sabotage_event["true_zero_index"])
        print("True zero T               :", mp.nstr(sabotage_event["true_zero_t"], 30))
        print("Injected bad zero T       :", mp.nstr(sabotage_event["bad_zero_t"], 30))
        print("Energy at bad zero        :", mp.nstr(sabotage_event["energy_at_bad_zero"], 30))
        print("Energy at true zero       :", mp.nstr(sabotage_event["real_energy_at_true_zero"], 30))
        print("Correction debt           :", mp.nstr(sabotage_event["correction"], 30))
        print("Forced profile at bad zero:", mp.nstr(energy(sabotage_event["bad_zero_t"]) - sabotage_event["correction"], 30))
        print("Energy-space deficit notes:")
        for i, chamber in enumerate(chambers, start=1):
            deficit_norm = sabotage_event["correction"] / chamber["antinode_e"]
            print("  Chamber {} deficit depth : {}".format(i, mp.nstr(deficit_norm, 16)))
        print()

    make_plot(chambers, sabotage_event=sabotage_event)


if __name__ == "__main__":
    main()
