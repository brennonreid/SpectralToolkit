#!/usr/bin/env python3
from mpmath import mp
import argparse

# Default precision (override via CLI)
mp.dps = 600


# -----------------------------
# Core Functions
# -----------------------------
def get_abs_zeta(t_val):
    return mp.fabs(mp.zeta(mp.mpc("0.5", t_val)))


def strict_truncate(val, decimals):
    """
    True strict truncation to `decimals` digits after the decimal point.
    No rounding, ever. Truncates toward 0.
    """
    if decimals < 0:
        raise ValueError("decimals must be >= 0")

    if decimals == 0:
        return mp.sign(val) * mp.floor(mp.fabs(val))

    scale = mp.power(10, decimals)
    return mp.sign(val) * (mp.floor(mp.fabs(val) * scale) / scale)


def fixed_trunc_str(val, decimals):
    """
    Render `val` as a fixed-point decimal string with exactly `decimals` digits
    after the decimal point, using the same truncation rule (no rounding).
    """
    if decimals < 0:
        raise ValueError("decimals must be >= 0")

    neg = (val < 0)
    aval = mp.fabs(val)

    if decimals == 0:
        iv = int(mp.floor(aval))
        return ("-" if neg and iv != 0 else "") + str(iv)

    scale = mp.power(10, decimals)
    scaled = mp.floor(aval * scale)
    scaled_int = int(scaled)

    base = 10 ** decimals
    int_part = scaled_int // base
    frac_part = scaled_int % base

    frac_str = str(frac_part).rjust(decimals, "0")
    sign = "-" if neg else ""
    return f"{sign}{int_part}.{frac_str}"


# -----------------------------
# CLI
# -----------------------------
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--t", default="14.134", help="Starting t address (string). Can be low precision like 14.1")
    p.add_argument("--max_steps", type=int, default=200, help="Global iteration limit.")
    p.add_argument("--mp_dps", type=int, default=600, help="mpmath precision (mp.dps).")
    return p.parse_args()


def main():
    args = parse_args()
    mp.dps = int(args.mp_dps)

    # INPUT
    t_string = str(args.t).strip()
    t_current = mp.mpf(t_string)

    # Auto-detect starting decimals
    if "." in t_string:
        current_decimals = len(t_string.split(".", 1)[1])
    else:
        current_decimals = 0

    print(f"{'Lvl':<4} | {'Step':<6} | {'Grow':<6} | {'Status':<8} | {'t_address'} | {'Energy'}")
    print("-" * 200)

    total_steps = 0
    max_steps = int(args.max_steps)

    while total_steps < max_steps:
        total_steps += 1

        # ------------------------------------------------------------
        # NEW LOOP RULES (per your spec)
        #
        # Phase A: while decimals < 5, grow by +1 each iteration until 5
        # Phase B: once decimals >= 5, grow by +current_decimals each time
        # Cap the growth step at 25 (max step_size)
        # ------------------------------------------------------------
        if current_decimals < 5:
            step_size = 1
            grow_mode = "ramp"
        else:
            step_size = current_decimals
            grow_mode = "pow"

        if step_size > 25:
            step_size = 25

        target_decimals = current_decimals + step_size

        # 2. PROBE & SLOPE (same math as before)
        probe_depth = target_decimals + 5
        h = mp.power(10, -probe_depth)

        e_now = get_abs_zeta(t_current)
        if e_now == 0:
            print(">>> PERFECT ZERO <<<")
            break

        e_plus = get_abs_zeta(t_current + h)

        diff = e_plus - e_now
        if diff == 0:
            slope = mp.mpf("1e-300")
        else:
            slope = diff / h

        raw_jump = (e_now / slope)

        # 3. THE CLAMP (same math as before)
        movement_cap = mp.power(10, -current_decimals)

        if mp.fabs(raw_jump) > movement_cap:
            # HIT THE WALL
            safe_jump = movement_cap * mp.sign(raw_jump)
            jump_note = "[CAP]"
            next_decimals_state = current_decimals
            display_decimals = current_decimals
            truncate_decimals = current_decimals  # never expand on CAP
        else:
            # FITS IN TRENCH
            safe_jump = raw_jump
            jump_note = "[FIT]"
            next_decimals_state = target_decimals
            display_decimals = target_decimals
            truncate_decimals = target_decimals

        t_next = t_current - safe_jump

        # 4. TRUNCATION
        t_current = strict_truncate(t_next, truncate_decimals)

        # OUTPUT
        e_final = get_abs_zeta(t_current)
        t_print = fixed_trunc_str(t_current, display_decimals)
        e_print = mp.nstr(e_final, 45)
        print(f"{current_decimals:<4} | {step_size:<6} | {grow_mode:<6} | {jump_note:<8} | {t_print} | {e_print}")

        # Update State
        current_decimals = next_decimals_state

        # Stop immediately on CAP (kept as your current policy)
        if jump_note == "[CAP]":
            print(">>> HIT CAP: stopping (resolution wall reached at current_decimals) <<<")
            break

        # Exit conditions (unchanged)
        if current_decimals >= mp.dps - 20:
            print(">>> REACHED MAX PRECISION LIMIT <<<")
            break

        if e_final < mp.power(10, -(mp.dps - 20)):
            print(">>> CONVERGED <<<")
            break

    print("-" * 200)
    final_e = get_abs_zeta(t_current)
    print(f"FINAL SINK ENERGY: {mp.nstr(final_e, 150)}")
    print(f"FINAL T VALUE:     {fixed_trunc_str(t_current, current_decimals)}")


if __name__ == "__main__":
    main()
