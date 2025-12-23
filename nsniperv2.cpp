// nsniper.cpp
//
// Behavior (per your spec):
// - Always apply the jump (CAP or FIT), regardless of size.
// - Track fractional digits only (Lvl column).
// - Internal working fractional precision can overshoot depth.
// - After each jump, if work_frac_digits >= depth:
//     * If == depth: compute energy at t_work and finish.
//     * If  > depth: strict-truncate (no rounding) to exactly depth frac digits,
//                    compute energy at truncated value, and finish.
//
// Additional rule (your last request):
// - If launched with --dir down (or up), force that direction until the first "lock" occurs.
// - Definition of "lock": first time CAP engages (|raw_jump| > cap).
// - Once lock is found, permanently switch direction to auto.
//
// Important directional semantics:
// - We apply direction forcing to BOTH CAP and FIT before lock.
//   This guarantees "it should go down" until lock.
// - Update step is: t_next = t_work - safe_jump
//     * For DOWN: safe_jump must be +|jump|  (t decreases)
//     * For UP:   safe_jump must be -|jump|  (t increases)
//     * For AUTO: safe_jump is chosen by CAP rule or Newton FIT sign
//
// Notes:
// - No rounding anywhere in the final truncation.
// - No "stall" logic at depth; overshoot mode makes "jump < 1e-depth" irrelevant.
// - Minimal changes; no speculative refactors.
//
// CHANGE (per your request):
// 1) Replace hard-coded step cap (25) with dynamic --max_step (default 25).
// 2) Reinstate truncation rule: in pow-mode growth is "double minus 1 digit".
// 3) Add FIT "shortcut": if a FIT jump magnitude implies far more precision than
//    the planned step schedule, accept it and mark output status as [FIT*].
// 4) FINAL WRAP-UP: when depth is met, print a qualified FINAL_SLOPE(mid) using
//    a depth-tied probe step, without changing the per-row slope display.
//
// DPS rigor preserved:
// - compute_dps still based on integer_digits + internal_max_frac + guard(64).
// - any output conversion sig-hint guard tied to args.max_step.

#include <iostream>
#include <string>
#include <iomanip>
#include <cstdlib>

#include "decimal.hpp"
#include "arb_interface.hpp"

// -----------------------------
// Helpers
// -----------------------------

// strict_truncate(val, decimals): truncates toward 0, never rounds.
static Decimal strict_truncate_toward0(const Decimal &x, int decimals) {
    if (decimals < 0) {
        std::cerr << "ERROR: decimals must be >= 0\n";
        std::exit(1);
    }

    if (decimals == 0) {
        mpz_class n = x.get_num();
        mpz_class d = x.get_den();
        if (n >= 0) return Decimal(n / d);
        mpz_class nn = -n;
        mpz_class fl = nn / d;
        return Decimal(-fl);
    }

    mpz_class scale = pow10_z(decimals);
    Decimal y = x * Decimal(scale);

    mpz_class n = y.get_num();
    mpz_class d = y.get_den();

    mpz_class q;
    if (n >= 0) {
        q = n / d;
    } else {
        mpz_class nn = -n;
        mpz_class fl = nn / d;
        q = -fl;
    }

    return Decimal(q, scale);
}

// fixed-point string with exactly `decimals` digits (no rounding).
static std::string fixed_trunc_str_decimal_exact(const Decimal &val, int decimals) {
    if (decimals < 0) {
        std::cerr << "ERROR: decimals must be >= 0\n";
        std::exit(1);
    }

    Decimal x = val;
    bool neg = (x < Decimal(0));
    if (neg) x = -x;

    mpz_class int_part;
    mpz_class frac_part;
    mpz_class scale;

    if (decimals == 0) {
        Decimal t = strict_truncate_toward0(x, 0);
        std::string s = decimal_to_string(t, 0);
        if (neg && s != "0") return "-" + s;
        return s;
    }

    scale = pow10_z(decimals);

    // scaled_int = trunc(x * scale)
    Decimal scaled = strict_truncate_toward0(x * Decimal(scale), 0);
    mpz_class scaled_int = scaled.get_num() / scaled.get_den();

    int_part = scaled_int / scale;
    frac_part = scaled_int % scale;

    std::string int_str = int_part.get_str();
    std::string frac_str = frac_part.get_str();
    if ((int)frac_str.size() < decimals) {
        frac_str.insert(frac_str.begin(), (size_t)(decimals - (int)frac_str.size()), '0');
    }

    std::string out = int_str + "." + frac_str;
    if (neg) out = "-" + out;
    return out;
}

// Convert decimal digits to Arb bits precision.
static slong dps_to_bits(int dps) {
    long double bits = (long double)dps * 3.32192809488736234787L;
    bits += 64.0L;
    if (bits < 128.0L) bits = 128.0L;
    return (slong)bits;
}

// Scientific-notation aware Decimal parser wrapper.
// Accepts: [-]digits[.digits][e[+-]digits]
static Decimal decimal_from_string_allow_exp(const std::string &s) {
    return decimal_from_string(s);
}

// Convert arb midpoint to a string using arf_get_str(mid, sig_digits).
static std::string arb_mid_to_string_sig(const arb_t x, int sig_digits) {
    const arf_struct *mid = arb_midref(x);
    char *cstr = arf_get_str(mid, sig_digits);
    std::string s(cstr);
    flint_free(cstr);
    return s;
}

// Convert arb midpoint to Decimal and strict-truncate toward 0 to `decimals`.
static Decimal arb_mid_to_decimal_strict_trunc(const arb_t x, int decimals, int sig_digits_hint) {
    const arf_struct *mid = arb_midref(x);

    char *cstr = arf_get_str(mid, sig_digits_hint);
    std::string s(cstr);
    flint_free(cstr);

    Decimal v = decimal_from_string_allow_exp(s);
    return strict_truncate_toward0(v, decimals);
}

// Count integer digits (no sign) of current Decimal.
static int count_int_digits(const Decimal &x) {
    Decimal t_int = strict_truncate_toward0(x, 0);
    std::string is = decimal_to_string(t_int, 0);
    if (!is.empty() && is[0] == '-') is.erase(0, 1);
    int d = (int)is.size();
    if (d < 1) d = 1;
    return d;
}

// Parse base-10 exponent from an "m e exp" style string (e.g. "-7.27e-841").
// Returns true if exponent parsed; false if not found.
static bool parse_exp10_from_sci(const std::string &s, long long &out_exp10) {
    size_t epos = s.find('e');
    if (epos == std::string::npos) epos = s.find('E');
    if (epos == std::string::npos) return false;

    std::string exp_part = s.substr(epos + 1);
    if (exp_part.empty()) return false;

    try {
        out_exp10 = std::stoll(exp_part);
        return true;
    } catch (...) {
        return false;
    }
}

// -----------------------------
// CLI
// -----------------------------
struct Args {
    std::string t = "14.134";
    int max_steps = 200;

    // user-facing: exact fractional digits requested
    int depth = 400;

    // Direction control until lock; after lock becomes auto
    std::string dir = "auto";   // auto|up|down

    // Overshoot fractional digits beyond depth for internal work state
    int overshoot_frac = 32;

    // Dynamic cap for step_size (replaces hard-coded 25)
    int max_step = 25;
};

static Args parse_args(int argc, char **argv) {
    Args a;
    for (int i = 1; i < argc; i++) {
        std::string k = argv[i];

        auto require_value = [&](const std::string &flag) {
            if (i + 1 >= argc) {
                std::cerr << "ERROR: Missing value after " << flag << "\n";
                std::exit(1);
            }
        };

        if (k == "--t") {
            require_value(k);
            a.t = argv[++i];
        } else if (k == "--max_steps") {
            require_value(k);
            a.max_steps = std::stoi(argv[++i]);
        } else if (k == "--depth") {
            require_value(k);
            a.depth = std::stoi(argv[++i]);
            if (a.depth < 0) {
                std::cerr << "ERROR: --depth must be >= 0\n";
                std::exit(1);
            }
        } else if (k == "--dir") {
            require_value(k);
            a.dir = argv[++i];
            if (a.dir != "auto" && a.dir != "up" && a.dir != "down") {
                std::cerr << "ERROR: --dir must be one of: auto | up | down\n";
                std::exit(1);
            }
        } else if (k == "--overshoot_frac") {
            require_value(k);
            a.overshoot_frac = std::stoi(argv[++i]);
            if (a.overshoot_frac < 0) {
                std::cerr << "ERROR: --overshoot_frac must be >= 0\n";
                std::exit(1);
            }
        } else if (k == "--max_step") {
            require_value(k);
            a.max_step = std::stoi(argv[++i]);
            if (a.max_step < 1) {
                std::cerr << "ERROR: --max_step must be >= 1\n";
                std::exit(1);
            }
        } else {
            std::cerr << "ERROR: Unknown argument: " << k << "\n";
            std::exit(1);
        }
    }
    return a;
}

// -----------------------------
// Main
// -----------------------------
int main(int argc, char **argv) {
    Args args = parse_args(argc, argv);

    std::string t_string = trim(args.t);
    Decimal t_work = decimal_from_string_allow_exp(t_string);

    // starting fractional digits tracked (from input string)
    int work_frac_digits = count_decimals_str(t_string);

    // Internal max fractional digits we may keep in t_work
    const int internal_max_frac = args.depth + args.overshoot_frac;

    // Choose compute dps based on integer digits + internal_max_frac + guard(64)
    int start_int_digits = count_int_digits(t_work);
    int compute_dps = (start_int_digits + internal_max_frac) + 64;

    slong prec_bits = dps_to_bits(compute_dps);
    ArbZetaContext ctx(prec_bits);

    const Decimal sigma = Decimal(1, 2);

    // Direction starts as user-specified, but switches to auto after first CAP ("lock").
    std::string dir_effective = args.dir;
    bool lock_found = false;

    // Startup banner
    std::cout << "DEPTH (exact frac digits): " << args.depth << "\n";
    std::cout << "OVERSHOOT_FRAC:            " << args.overshoot_frac << "\n";
    std::cout << "INTERNAL_MAX_FRAC:         " << internal_max_frac << " (depth + overshoot)\n";
    std::cout << "MAX_STEP:                  " << args.max_step << "\n";
    std::cout << "START_INT_DIGITS:          " << start_int_digits << "\n";
    std::cout << "COMPUTE_DPS:               " << compute_dps << "\n";
    std::cout << "ARB PREC (bits):           " << prec_bits << "\n\n";

    std::cout << std::left
              << std::setw(4)  << "Lvl" << " | "
              << std::setw(6)  << "Step" << " | "
              << std::setw(6)  << "Grow" << " | "
              << std::setw(8)  << "Status" << " | "
              << std::setw(28) << "Slope(mid)" << " | "
              << std::setw(28) << "Jump(mid)" << " | "
              << "t_address"
              << " | "
              << "Energy"
              << "\n";
    std::cout << std::string(220, '-') << "\n";

    int total_steps = 0;

    arb_t e_now, e_plus, diff, slope, raw_jump, abs_raw_jump;
    arb_t h_arb, cap_arb, safe_jump, t_work_arb, t_next_arb;
    arb_init(e_now);
    arb_init(e_plus);
    arb_init(diff);
    arb_init(slope);
    arb_init(raw_jump);
    arb_init(abs_raw_jump);
    arb_init(h_arb);
    arb_init(cap_arb);
    arb_init(safe_jump);
    arb_init(t_work_arb);
    arb_init(t_next_arb);

    while (total_steps < args.max_steps) {
        total_steps += 1;

        // Growth schedule based on current fractional digits tracked
        // IMPORTANT: pow-mode uses "double minus 1 digit" rule.
        int step_size_raw;
        std::string grow_mode;
        if (work_frac_digits < 5) {
            step_size_raw = 1;
            grow_mode = "ramp";
        } else {
            step_size_raw = work_frac_digits;
            grow_mode = "pow";
        }

        // Cap intended step by --max_step
        int step_size_capped = step_size_raw;
        if (step_size_capped > args.max_step) step_size_capped = args.max_step;

        // Apply truncation rule:
        // - ramp: +1
        // - pow:  +(capped - 1)  (double minus 1 digit)
        int step_size_eff;
        if (work_frac_digits < 5) {
            step_size_eff = 1;
        } else {
            step_size_eff = step_size_capped - 1;
            if (step_size_eff < 1) step_size_eff = 1;
        }

        int target_frac = work_frac_digits + step_size_eff;
        if (target_frac > internal_max_frac) target_frac = internal_max_frac;

        // Probe & slope depth (use target_frac + small guard)
        int probe_depth = target_frac + 5;
        Decimal h_dec(1, pow10_z(probe_depth));
        ctx.arb_set_decimal_safe(h_arb, h_dec);

        ctx.compute(t_work, sigma);
        arb_set(e_now, ctx.abs_z);

        if (arb_is_zero(e_now)) {
            // Perfect zero at current working state; still honor depth rule for output.
            break;
        }

        Decimal t_plus = t_work + h_dec;
        ctx.compute(t_plus, sigma);
        arb_set(e_plus, ctx.abs_z);

        arb_sub(diff, e_plus, e_now, prec_bits);
        if (arb_is_zero(diff)) {
            arb_set_str(slope, "1e-300", prec_bits);
        } else {
            arb_div(slope, diff, h_arb, prec_bits);
        }

        arb_div(raw_jump, e_now, slope, prec_bits);

        // Clamp (movement cap = 10^-work_frac_digits; integer start uses 0.1)
        Decimal cap_dec(1, pow10_z(work_frac_digits));
        if (work_frac_digits == 0) {
            cap_dec = Decimal(1, 10); // 0.1
        }
        ctx.arb_set_decimal_safe(cap_arb, cap_dec);

        arb_abs(abs_raw_jump, raw_jump);
        bool is_cap = arb_gt(abs_raw_jump, cap_arb);

        // "Lock" rule: first time CAP engages, switch to auto permanently.
        if (is_cap && !lock_found) {
            lock_found = true;
            dir_effective = "auto";
        }

        std::string status;
        int next_work_frac = work_frac_digits;

        if (is_cap) {
            status = "[CAP]";

            // Direction control for CAP
            if (dir_effective == "down") {
                // want t_next = t - cap  => safe_jump = +cap
                arb_set(safe_jump, cap_arb);
            } else if (dir_effective == "up") {
                // want t_next = t + cap  => safe_jump = -cap
                arb_neg(safe_jump, cap_arb);
            } else {
                // auto: compare energy at +/- cap
                arb_t e_cap_plus, e_cap_minus;
                arb_init(e_cap_plus);
                arb_init(e_cap_minus);

                Decimal t_cap_plus  = t_work + cap_dec;
                Decimal t_cap_minus = t_work - cap_dec;

                ctx.compute(t_cap_plus, sigma);
                arb_set(e_cap_plus, ctx.abs_z);

                ctx.compute(t_cap_minus, sigma);
                arb_set(e_cap_minus, ctx.abs_z);

                if (arb_lt(e_cap_plus, e_cap_minus)) {
                    // choose t_next = t + cap => safe_jump = -cap
                    arb_neg(safe_jump, cap_arb);
                } else {
                    // choose t_next = t - cap => safe_jump = +cap
                    arb_set(safe_jump, cap_arb);
                }

                arb_clear(e_cap_plus);
                arb_clear(e_cap_minus);
            }

            // On integer CAP, ensure at least 1 fractional digit is tracked going forward.
            if (work_frac_digits == 0) {
                next_work_frac = 1;
            } else {
                // CAP does not automatically increase precision; keep as-is.
                next_work_frac = work_frac_digits;
            }
        } else {
            status = "[FIT]";

            // Default Newton FIT jump
            arb_set(safe_jump, raw_jump);
            next_work_frac = target_frac;

            // Before lock, force user direction for FIT as well (guarantees monotone down/up).
            if (!lock_found) {
                if (dir_effective == "down") {
                    // safe_jump = +|raw_jump|
                    arb_abs(safe_jump, raw_jump);
                } else if (dir_effective == "up") {
                    // safe_jump = -|raw_jump|
                    arb_abs(safe_jump, raw_jump);
                    arb_neg(safe_jump, safe_jump);
                } else {
                    // auto: keep Newton sign
                }
            }

            // --- FIT SHORTCUT LOGIC ---
            // If the FIT jump magnitude implies much deeper decimal addressability than the
            // planned next_work_frac, allow it (still a FIT, still within cap criteria).
            //
            // We infer implied fractional digits from scientific exponent of jump midpoint:
            //   |jump| ~= m * 10^(exp10)
            // If exp10 is negative, it suggests the jump "touches" around -exp10 decimals.
            {
                arb_t abs_jump;
                arb_init(abs_jump);
                arb_abs(abs_jump, safe_jump);

                // Use a small sig digit count; we only need the exponent.
                std::string aj = arb_mid_to_string_sig(abs_jump, 20);
                arb_clear(abs_jump);

                long long exp10 = 0;
                if (parse_exp10_from_sci(aj, exp10)) {
                    if (exp10 < 0) {
                        long long implied = -exp10;
                        if (implied > (long long)next_work_frac) {
                            if (implied > (long long)internal_max_frac) implied = (long long)internal_max_frac;
                            next_work_frac = (int)implied;
                            status = "[FIT*]";
                        }
                    }
                }
            }
        }

        // Apply step: t_next = t_work - safe_jump
        ctx.arb_set_decimal_safe(t_work_arb, t_work);
        arb_sub(t_next_arb, t_work_arb, safe_jump, prec_bits);

        // Convert next arb midpoint to Decimal at INTERNAL working precision (next_work_frac),
        // NOT at args.depth. This is the critical separation.
        {
            int int_digits_now = count_int_digits(t_work);

            // guard digits: tie to max_step (matches your earlier dps/guard reasoning).
            // We keep this policy even with shortcut; shortcut precision is bounded by internal_max_frac.
            int sig_hint = int_digits_now + next_work_frac + args.max_step;

            t_work = arb_mid_to_decimal_strict_trunc(t_next_arb, next_work_frac, sig_hint);
        }

        // Output: address is shown at min(next_work_frac, args.depth) for readability.
        int display_decimals = next_work_frac;
        if (display_decimals > args.depth) display_decimals = args.depth;

        std::string t_addr_print = fixed_trunc_str_decimal_exact(t_work, display_decimals);

        // Keep slope/jump prints stable; slope shown with 50 sig digits for run display.
        std::string slope_mid = arb_mid_to_string_sig(slope, 50);
        std::string jump_mid  = arb_mid_to_string_sig(safe_jump, 20);

        ctx.compute(t_work, sigma);
        std::string e_print = ctx.get_energy_str(45);

        // Step column prints the planned "velocity" step (effective step schedule),
        // not the shortcut-implied digits. Shortcut is indicated via Status [FIT*].
        std::cout << std::left
                  << std::setw(4)  << work_frac_digits << " | "
                  << std::setw(6)  << step_size_eff << " | "
                  << std::setw(6)  << grow_mode << " | "
                  << std::setw(8)  << status << " | "
                  << std::setw(28) << slope_mid << " | "
                  << std::setw(28) << jump_mid << " | "
                  << t_addr_print << " | "
                  << e_print
                  << "\n";

        // Update fractional digit tracking state
        work_frac_digits = next_work_frac;

        // --- YOUR DEPTH TERMINATION RULE (exactly as requested) ---
        // After jump, if fractional digits tracked >= depth:
        //   - if == depth: energy at t_work (already computed) and finish
        //   - if  > depth: truncate to depth (no rounding), compute energy at truncated, finish
        if (work_frac_digits >= args.depth) {
            Decimal t_final;
            if (work_frac_digits == args.depth) {
                t_final = t_work;
            } else {
                t_final = strict_truncate_toward0(t_work, args.depth);
            }

            // Keep your existing final energy evaluation.
            ctx.compute(t_final, sigma);

            // --- FINAL SLOPE (qualified; NOT the 50-sig per-row display) ---
            // Compute slope at t_final using a depth-tied probe: h = 10^-(depth+5),
            // clamped to internal_max_frac.
            int final_probe_depth = args.depth + 5;
            if (final_probe_depth > internal_max_frac) final_probe_depth = internal_max_frac;
            if (final_probe_depth < 1) final_probe_depth = 1;

            Decimal h_final_dec(1, pow10_z(final_probe_depth));

            arb_t e0f, e1f, difff, hf, slope_final;
            arb_init(e0f);
            arb_init(e1f);
            arb_init(difff);
            arb_init(hf);
            arb_init(slope_final);

            ctx.arb_set_decimal_safe(hf, h_final_dec);

            // e0 = E(t_final)
            ctx.compute(t_final, sigma);
            arb_set(e0f, ctx.abs_z);

            // e1 = E(t_final + h)
            Decimal t_final_plus = t_final + h_final_dec;
            ctx.compute(t_final_plus, sigma);
            arb_set(e1f, ctx.abs_z);

            // slope_final = (e1 - e0) / h
            arb_sub(difff, e1f, e0f, prec_bits);
            if (arb_is_zero(difff)) {
                arb_set_str(slope_final, "1e-300", prec_bits);
            } else {
                arb_div(slope_final, difff, hf, prec_bits);
            }

            // Print size for final slope (qualified)
            // Policy: use max(80, depth/2) sig digits, capped by compute_dps.
            int final_slope_sig = args.depth / 2;
            if (final_slope_sig < 80) final_slope_sig = 80;
            if (final_slope_sig > compute_dps) final_slope_sig = compute_dps;

            std::cout << ">>> REACHED TARGET DEPTH (SETTLED) <<<\n";
            std::cout << std::string(220, '-') << "\n";
            std::cout << "FINAL SINK ENERGY: " << ctx.get_energy_str(150) << "\n";
            std::cout << "FINAL T VALUE:     " << fixed_trunc_str_decimal_exact(t_final, args.depth) << "\n";
            std::cout << "FINAL_SLOPE(mid):  " << arb_mid_to_string_sig(slope_final, final_slope_sig) << "\n";
            std::cout << "FINAL_SLOPE_SIG:   " << final_slope_sig << "\n";
            std::cout << "FINAL_PROBE_H:     1e-" << final_probe_depth << "\n";

            arb_clear(e0f);
            arb_clear(e1f);
            arb_clear(difff);
            arb_clear(hf);
            arb_clear(slope_final);

            arb_clear(e_now);
            arb_clear(e_plus);
            arb_clear(diff);
            arb_clear(slope);
            arb_clear(raw_jump);
            arb_clear(abs_raw_jump);
            arb_clear(h_arb);
            arb_clear(cap_arb);
            arb_clear(safe_jump);
            arb_clear(t_work_arb);
            arb_clear(t_next_arb);

            return 0;
        }
    }

    // If we fall out (max_steps) without reaching depth, report best we have.
    std::cout << std::string(220, '-') << "\n";
    ctx.compute(t_work, sigma);
    std::cout << "FINAL SINK ENERGY: " << ctx.get_energy_str(150) << "\n";
    std::cout << "FINAL T VALUE:     " << fixed_trunc_str_decimal_exact(t_work, work_frac_digits) << "\n";

    arb_clear(e_now);
    arb_clear(e_plus);
    arb_clear(diff);
    arb_clear(slope);
    arb_clear(raw_jump);
    arb_clear(abs_raw_jump);
    arb_clear(h_arb);
    arb_clear(cap_arb);
    arb_clear(safe_jump);
    arb_clear(t_work_arb);
    arb_clear(t_next_arb);

    return 0;
}
