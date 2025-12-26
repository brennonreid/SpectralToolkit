// nsniperv2.cpp
//
// Closure-stable Johansson-digit matching version.
// Fix: DO NOT convert Arb values through arf_get_str()/scientific notation for state updates.
// Use arf_get_fmpz_2exp() to convert jump midpoints exactly to Decimal,
// then apply update in Decimal and strict-truncate (never round).
//
// Policy:
// - State update is ALWAYS Decimal arithmetic.
// - CAP is exact Decimal +/- 10^-k (integer start uses 0.1).
// - FIT jump computed in Arb, converted exactly via mantissa+2exp, then truncated.
// - After each update, strict truncate t_work to intended fractional digits.
// - When reaching depth: jump, then truncate to exactly depth, stop. No extra adjustments.
//
// NEW (your request):
// - Print FINAL_SLOPE(mid) at FULL FIXED DECIMAL precision (exactly args.depth decimals),
//   with NO rounding and NO scientific notation.
//
// NEW (your request #2):
// - Track ratio of energy loss each step:
//     E_ratio = E_next / E_now
//     Drop%   = (1 - E_ratio) * 100
// - Display per step
// - Cache last ratio + drop% and print in wrap-up (fixed, args.depth decimals, truncate only)
//
// NEW (print fixes, your request now):
// 1) FINAL SINK ENERGY should not print as 0 due to args.depth truncation.
//    Use DISPLAY-ONLY adaptive fixed formatting (no scientific notation) so small values show.
// 2) Final output should avoid trailing zeros unless necessary (DISPLAY-ONLY).
//    Keep FINAL_SLOPE(mid) and FINAL T VALUE strict fixed-width; trim ratio/drop and final energy.
//
// Build:
//   g++ -O2 -std=c++17 nsniperv2.cpp -lflint -lgmp -lmpfr -o nsniperv2

#include <iostream>
#include <string>
#include <iomanip>
#include <cstdlib>
#include <cctype>

#include "decimal.hpp"
#include "arb_interface.hpp"

#include <flint/arf.h>
#include <flint/fmpz.h>

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

    if (decimals == 0) {
        Decimal t = strict_truncate_toward0(x, 0);
        std::string s = decimal_to_string(t, 0);
        if (neg && s != "0") return "-" + s;
        return s;
    }

    mpz_class scale = pow10_z(decimals);

    // scaled_int = trunc(x * scale)
    Decimal scaled = strict_truncate_toward0(x * Decimal(scale), 0);
    mpz_class scaled_int = scaled.get_num() / scaled.get_den();

    mpz_class int_part = scaled_int / scale;
    mpz_class frac_part = scaled_int % scale;

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

// Count integer digits (no sign) of current Decimal.
static int count_int_digits(const Decimal &x) {
    Decimal t_int = strict_truncate_toward0(x, 0);
    std::string is = decimal_to_string(t_int, 0);
    if (!is.empty() && is[0] == '-') is.erase(0, 1);
    int d = (int)is.size();
    if (d < 1) d = 1;
    return d;
}

// Arb midpoint to string (DISPLAY ONLY).
static std::string arb_mid_to_string_sig(const arb_t x, int sig_digits) {
    const arf_struct *mid = arb_midref(x);
    char *cstr = arf_get_str(mid, sig_digits);
    if (!cstr) {
        std::cerr << "ERROR: arf_get_str failed\n";
        std::exit(1);
    }
    std::string s(cstr);
    flint_free(cstr);
    return s;
}

// Convert arf exactly via mantissa * 2^exp2, then strict-truncate toward 0 to `decimals`.
static Decimal arf_mid_to_decimal_exact_trunc_pow10(const arf_struct *mid, int decimals) {
    if (decimals < 0) {
        std::cerr << "ERROR: decimals must be >= 0\n";
        std::exit(1);
    }

    if (arf_is_zero(mid)) return Decimal(0);

    fmpz_t m_f;
    fmpz_init(m_f);

    slong e2 = 0;
    arf_get_fmpz_2exp(m_f, &e2, mid);

    char *m_str = fmpz_get_str(NULL, 10, m_f);
    if (!m_str) {
        std::cerr << "ERROR: fmpz_get_str failed\n";
        std::exit(1);
    }
    mpz_class m_mpz(m_str);
    flint_free(m_str);
    fmpz_clear(m_f);

    Decimal r;
    if (e2 >= 0) {
        mpz_class two_pow = mpz_class(1) << (unsigned long)e2;
        r = Decimal(m_mpz * two_pow);
    } else {
        mpz_class two_pow = mpz_class(1) << (unsigned long)(-e2);
        r = Decimal(m_mpz, two_pow);
    }

    return strict_truncate_toward0(r, decimals);
}

static Decimal arb_mid_to_decimal_exact_trunc(const arb_t x, int decimals) {
    const arf_struct *mid = arb_midref(x);
    return arf_mid_to_decimal_exact_trunc_pow10(mid, decimals);
}

// -----------------------------
// DISPLAY helpers (no effect on state updates)
// -----------------------------

// Trim trailing zeros after decimal point; remove '.' if it becomes the last char.
// Keeps sign intact.
static std::string trim_trailing_zeros_fixed(std::string s) {
    if (s.find('e') != std::string::npos || s.find('E') != std::string::npos) return s;

    size_t dot = s.find('.');
    if (dot == std::string::npos) return s;

    while (!s.empty() && s.back() == '0') s.pop_back();
    if (!s.empty() && s.back() == '.') s.pop_back();

    if (s == "-0") s = "0";
    return s;
}

// Convert scientific string like "-3.8202792979e-40" to fixed decimal (no exponent).
// DISPLAY ONLY (arf_get_str rounds). If max_decimals >= 0, truncate to at most that many.
static std::string sci_to_fixed_noexp_trunc(const std::string &sci, int max_decimals) {
    std::string s = sci;
    if (s == "0" || s == "+0" || s == "-0") return "0";

    bool neg = false;
    if (!s.empty() && (s[0] == '+' || s[0] == '-')) {
        neg = (s[0] == '-');
        s.erase(0, 1);
    }

    std::string mant = s;
    long long exp10 = 0;

    size_t epos = s.find_first_of("eE");
    if (epos != std::string::npos) {
        mant = s.substr(0, epos);
        std::string ex = s.substr(epos + 1);
        exp10 = std::stoll(ex);
    }

    long long frac_digits = 0;
    size_t dpos = mant.find('.');
    if (dpos != std::string::npos) {
        frac_digits = (long long)(mant.size() - dpos - 1);
        mant.erase(dpos, 1);
    }

    while (mant.size() > 1 && mant[0] == '0') mant.erase(0, 1);

    long long shift = exp10 - frac_digits;

    std::string out;
    if (shift >= 0) {
        out = mant;
        out.append((size_t)shift, '0');
    } else {
        long long k = -shift;
        if ((long long)mant.size() > k) {
            size_t cut = mant.size() - (size_t)k;
            out = mant.substr(0, cut) + "." + mant.substr(cut);
        } else {
            out = "0.";
            out.append((size_t)(k - (long long)mant.size()), '0');
            out += mant;
        }
    }

    if (neg && out != "0") out = "-" + out;

    if (max_decimals >= 0) {
        size_t dot = out.find('.');
        if (dot != std::string::npos) {
            size_t have = out.size() - dot - 1;
            if ((int)have > max_decimals) {
                out.erase(dot + 1 + (size_t)max_decimals);
                if (max_decimals == 0) out.erase(dot, 1);
            }
        }
    }

    return out;
}

// Prints: LABEL: [mid +/- rad]   (with `digits` significant digits)
// NOTE: display-only; does not affect state updates.
static void print_arb_bracket(const char* label, const arb_t x, slong digits)
{
    std::cout << label << ": ";
    arb_printn(x, digits, 0); // [mid +/- rad]
    std::cout << "\n";
}

// Arb midpoint -> fixed decimal, adaptive enough to avoid printing as 0
// when truncating to args.depth would underflow.
// DISPLAY ONLY.
static std::string arb_mid_to_fixed_decimal_adaptive_trim(const arb_t x, int sig_digits, int max_decimals) {
    const arf_struct *mid = arb_midref(x);

    if (arf_is_zero(mid)) return "0";

    // 1) First attempt: existing behavior
    {
        char *cstr = arf_get_str(mid, sig_digits);
        if (!cstr) {
            std::cerr << "ERROR: arf_get_str failed\n";
            std::exit(1);
        }
        std::string sci(cstr);
        flint_free(cstr);

        std::string fixed = sci_to_fixed_noexp_trunc(sci, max_decimals);
        fixed = trim_trailing_zeros_fixed(fixed);

        if (fixed != "0") return fixed;
    }

    // 2) If it would show as 0 but is nonzero, estimate decimal exponent
    //    using mantissa * 2^e2 from arf_get_fmpz_2exp().
    fmpz_t m_f;
    fmpz_init(m_f);
    slong e2 = 0;
    arf_get_fmpz_2exp(m_f, &e2, mid);

    // magnitude(m) in bits
    slong m_bits = fmpz_bits(m_f);
    if (m_bits < 1) m_bits = 1;

    // Rough log10(|m|) from bit-length: log10(m) ~= (bits-1)*log10(2)
    // Then log10(|mid|) ~= log10(|m|) + e2*log10(2)
    const long double LOG10_2 = 0.30102999566398119521L;

    long double log10_m_lo = (long double)(m_bits - 1) * LOG10_2;
    long double log10_mid  = log10_m_lo + (long double)e2 * LOG10_2;

    // decimal exponent floor
    long long exp10 = (long long) (log10_mid >= 0 ? (long long)log10_mid : (long long)(log10_mid - 1.0L));

    // If exp10 is very negative, we need about (-exp10 + guard) decimals to show first nonzero.
    long long needed_decimals = 0;
    if (exp10 < 0) {
        needed_decimals = (-exp10) + 8; // guard digits for readability
    } else {
        needed_decimals = 0;
    }

    fmpz_clear(m_f);

    // Respect a minimum of max_decimals (your old behavior), but allow more when needed.
    // Put a hard cap so logs don't explode accidentally.
    long long use_decimals_ll = (long long)max_decimals;
    if (needed_decimals > use_decimals_ll) use_decimals_ll = needed_decimals;

    const long long HARD_CAP = (long long)max_decimals + 1000; // display-only, bounded growth
    if (use_decimals_ll > HARD_CAP) use_decimals_ll = HARD_CAP;

    int use_decimals = (use_decimals_ll < 0) ? 0 : (int)use_decimals_ll;

    // 3) Re-render with more decimals
    char *c2 = arf_get_str(mid, sig_digits * 2);
    if (!c2) {
        std::cerr << "ERROR: arf_get_str failed\n";
        std::exit(1);
    }
    std::string sci2(c2);
    flint_free(c2);

    std::string fixed2 = sci_to_fixed_noexp_trunc(sci2, use_decimals);
    fixed2 = trim_trailing_zeros_fixed(fixed2);

    // If still "0" (possible if HARD_CAP too low), fall back to scientific string.
    if (fixed2 == "0") {
        return sci2; // last resort display-only
    }

    return fixed2;
}

// -----------------------------
// CLI
// -----------------------------
struct Args {
    std::string t = "14.134";
    int max_steps = 200;
    int depth = 400;
    int overshoot_frac = 32;
    int max_step = 25;
    std::string dir = "auto";
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
        } else if (k == "--dir") {
            require_value(k);
            a.dir = argv[++i];
            if (a.dir != "auto" && a.dir != "up" && a.dir != "down") {
                std::cerr << "ERROR: --dir must be one of: auto | up | down\n";
                std::exit(1);
            }
        } else {
            std::cerr << "ERROR: Unknown argument: " << k << "\n";
            std::exit(1);
        }
    }
    return a;
}

// Print arb midpoint as fixed decimal with exactly `decimals` digits after the point,
// truncation only, no rounding, no scientific notation. (May underflow to "0....0".)
static std::string arb_mid_to_fixed_decimal_trunc(const arb_t x, int decimals) {
    Decimal mid_dec = arb_mid_to_decimal_exact_trunc(x, decimals);
    return fixed_trunc_str_decimal_exact(mid_dec, decimals);
}

// -----------------------------
// Main
// -----------------------------
int main(int argc, char **argv) {
    Args args = parse_args(argc, argv);

    std::string t_string = trim(args.t);                 // from decimal.hpp
    Decimal t_work = decimal_from_string(t_string);

    int work_frac_digits = count_decimals_str(t_string); // from decimal.hpp

    const int internal_max_frac = args.depth + args.overshoot_frac;

    int start_int_digits = count_int_digits(t_work);
    int compute_dps = (start_int_digits + internal_max_frac) + args.max_step + 256;
    slong prec_bits = dps_to_bits(compute_dps);

    ArbZetaContext ctx(prec_bits);
    const Decimal sigma = Decimal(1, 2);

    std::string dir_effective = args.dir;
    bool lock_found = false;

    // cache slope(mid) for final report (Decimal, fixed output, no rounding).
    Decimal final_slope_mid_dec(0);
    bool have_final_slope = false;

    // cache final energy ratio/drop (Decimal, fixed output, no rounding).
    Decimal final_energy_ratio_dec(0);    // E_next / E_now
    Decimal final_energy_drop_pct_dec(0); // (1 - ratio) * 100
    bool have_final_energy_ratio = false;

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
              << std::setw(28) << "Jump(mid)"  << " | "
              << std::setw(18) << "E_ratio"    << " | "
              << std::setw(18) << "Drop%"      << " | "
              << "t_address"
              << " | "
              << "Energy"
              << "\n";
    std::cout << std::string(220, '-') << "\n";

    int total_steps = 0;

    arb_t e_now, e_plus, diff, slope, raw_jump, abs_raw_jump;
    arb_t h_arb, cap_arb, safe_jump;

    // energy-after + ratio/drop metrics
    arb_t e_after, e_ratio, e_drop_pct, one_arb, hundred_arb, tmp_arb;

    arb_init(e_now);
    arb_init(e_plus);
    arb_init(diff);
    arb_init(slope);
    arb_init(raw_jump);
    arb_init(abs_raw_jump);
    arb_init(h_arb);
    arb_init(cap_arb);
    arb_init(safe_jump);

    arb_init(e_after);
    arb_init(e_ratio);
    arb_init(e_drop_pct);
    arb_init(one_arb);
    arb_init(hundred_arb);
    arb_init(tmp_arb);

    arb_one(one_arb);
    arb_set_si(hundred_arb, 100);

    while (total_steps < args.max_steps) {
        total_steps += 1;

        int step_size_raw;
        std::string grow_mode;
        if (work_frac_digits < 5) {
            step_size_raw = 1;
            grow_mode = "ramp";
        } else {
            step_size_raw = work_frac_digits;
            grow_mode = "pow";
        }

        int step_size_capped = step_size_raw;
        if (step_size_capped > args.max_step) step_size_capped = args.max_step;

        int step_size_eff;
        if (work_frac_digits < 5) {
            step_size_eff = 1;
        } else {
            step_size_eff = step_size_capped - 1;
            if (step_size_eff < 1) step_size_eff = 1;
        }

        int target_frac = work_frac_digits + step_size_eff;
        if (target_frac > internal_max_frac) target_frac = internal_max_frac;

        int probe_depth;
        if (step_size_capped == args.max_step) {
            probe_depth = args.depth + 64;
        } else {
            probe_depth = target_frac + 5;
        }
        if (probe_depth > internal_max_frac) probe_depth = internal_max_frac;
        if (probe_depth < 1) probe_depth = 1;

        Decimal h_dec(1, pow10_z(probe_depth));
        ctx.arb_set_decimal_safe(h_arb, h_dec);

        ctx.compute(t_work, sigma);
        arb_set(e_now, ctx.abs_z);

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

        Decimal cap_dec(1, pow10_z(work_frac_digits));
        if (work_frac_digits == 0) cap_dec = Decimal(1, 10);
        ctx.arb_set_decimal_safe(cap_arb, cap_dec);

        arb_abs(abs_raw_jump, raw_jump);
        bool is_cap = arb_gt(abs_raw_jump, cap_arb);

        if (is_cap && !lock_found) {
            lock_found = true;
            dir_effective = "auto";
        }

        std::string status = is_cap ? "[CAP]" : "[FIT]";
        int next_work_frac = work_frac_digits;

        if (is_cap) {
            if (work_frac_digits == 0) next_work_frac = 1;
            else next_work_frac = work_frac_digits;

            if (dir_effective == "down") {
                arb_set(safe_jump, cap_arb);      // t_next = t - cap
            } else if (dir_effective == "up") {
                arb_neg(safe_jump, cap_arb);      // t_next = t + cap
            } else {
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
                    arb_neg(safe_jump, cap_arb); // move up
                } else {
                    arb_set(safe_jump, cap_arb); // move down
                }

                arb_clear(e_cap_plus);
                arb_clear(e_cap_minus);
            }
        } else {
            next_work_frac = target_frac;
            arb_set(safe_jump, raw_jump);

            if (!lock_found) {
                if (dir_effective == "down") {
                    arb_abs(safe_jump, raw_jump);      // t decreases
                } else if (dir_effective == "up") {
                    arb_abs(safe_jump, raw_jump);
                    arb_neg(safe_jump, safe_jump);     // t increases
                }
            }
        }

        // If this iteration will hit depth (or beyond), cache FULL slope(mid) in Decimal
        // at exactly args.depth decimals (fixed, no rounding).
        if (next_work_frac >= args.depth) {
            final_slope_mid_dec = arb_mid_to_decimal_exact_trunc(slope, args.depth);
            have_final_slope = true;
        }

        // ---- APPLY UPDATE IN DECIMAL SPACE ----
        Decimal t_next_dec = t_work;

        if (is_cap) {
            if (arb_is_negative(safe_jump)) {
                t_next_dec = t_work + cap_dec;
            } else {
                t_next_dec = t_work - cap_dec;
            }
        } else {
            // Convert safe_jump EXACTLY via mantissa+2exp (no strings),
            // then truncate the jump (never round) before applying.
            int jump_trunc_decimals = next_work_frac + 256;
            Decimal jump_dec = arb_mid_to_decimal_exact_trunc(safe_jump, jump_trunc_decimals);

            t_next_dec = t_work - jump_dec;
        }

        if (next_work_frac > internal_max_frac) next_work_frac = internal_max_frac;
        t_work = strict_truncate_toward0(t_next_dec, next_work_frac);

        int display_decimals = next_work_frac;
        if (display_decimals > args.depth) display_decimals = args.depth;

        std::string t_addr_print = fixed_trunc_str_decimal_exact(t_work, display_decimals);

        // DISPLAY ONLY (keep short)
        std::string slope_mid = arb_mid_to_string_sig(slope, 50);
        std::string jump_mid  = arb_mid_to_string_sig(safe_jump, 20);

        // Compute energy at updated t_work, and also ratio/drop vs previous e_now
        ctx.compute(t_work, sigma);
        arb_set(e_after, ctx.abs_z);
        std::string e_print = ctx.get_energy_str(45);

        // ratio = E_next / E_now; drop% = (1 - ratio) * 100
        if (arb_is_zero(e_now)) {
            arb_zero(e_ratio);
            arb_zero(e_drop_pct);
        } else {
            arb_div(e_ratio, e_after, e_now, prec_bits);
            arb_sub(tmp_arb, one_arb, e_ratio, prec_bits);
            arb_mul(e_drop_pct, tmp_arb, hundred_arb, prec_bits);
        }

        // Cache "final" ratio/drop for wrap-up (fixed Decimal, no rounding).
        final_energy_ratio_dec = arb_mid_to_decimal_exact_trunc(e_ratio, args.depth);
        final_energy_drop_pct_dec = arb_mid_to_decimal_exact_trunc(e_drop_pct, args.depth);
        have_final_energy_ratio = true;

        // DISPLAY ONLY (keep short)
        std::string ratio_mid = arb_mid_to_string_sig(e_ratio, 18);
        std::string drop_mid  = arb_mid_to_string_sig(e_drop_pct, 18);

        std::cout << std::left
                  << std::setw(4)  << work_frac_digits << " | "
                  << std::setw(6)  << step_size_eff << " | "
                  << std::setw(6)  << grow_mode << " | "
                  << std::setw(8)  << status << " | "
                  << std::setw(28) << slope_mid << " | "
                  << std::setw(28) << jump_mid  << " | "
                  << std::setw(18) << ratio_mid << " | "
                  << std::setw(18) << drop_mid  << " | "
                  << t_addr_print << " | "
                  << e_print
                  << "\n";

        work_frac_digits = next_work_frac;

        // Termination: jump then truncate to depth, stop.
        if (work_frac_digits >= args.depth) {
            Decimal t_final = t_work;
            if (work_frac_digits > args.depth) {
                t_final = strict_truncate_toward0(t_work, args.depth);
            }

            ctx.compute(t_final, sigma);
            // ctx.abs_z now holds FINAL SINK ENERGY (as arb_t) at t_final.

            std::cout << ">>> REACHED TARGET DEPTH (SETTLED) <<<\n";
            std::cout << std::string(220, '-') << "\n";

            // Scientific readout (mid +/- rad) immediately above fixed decimal, per your request.
            print_arb_bracket("FINAL SINK ENERGY (SCI, MID+/-RAD)", ctx.abs_z, 50);

            // FIX #1: FINAL SINK ENERGY should not display as 0 due solely to args.depth truncation.
            // DISPLAY ONLY: adaptive fixed decimal (no scientific notation), then trim trailing zeros.
            std::string final_energy_fixed = arb_mid_to_fixed_decimal_adaptive_trim(ctx.abs_z, 120, args.depth + 256);
            std::cout << "FINAL SINK ENERGY (MID, FIXED): " << final_energy_fixed << "\n";

            // Keep strict fixed-width for t_final (your paste discipline)
            std::cout << "FINAL T VALUE:     " << fixed_trunc_str_decimal_exact(t_final, args.depth) << "\n";

            // full slope(mid) at full fixed decimal precision (args.depth decimals), no rounding.
            if (!have_final_slope) {
                final_slope_mid_dec = arb_mid_to_decimal_exact_trunc(slope, args.depth);
                have_final_slope = true;
            }
            std::cout << "FINAL_SLOPE(mid):  " << fixed_trunc_str_decimal_exact(final_slope_mid_dec, args.depth) << "\n";

            // FIX #2: trim trailing zeros for final ratio/drop (DISPLAY ONLY).
            if (!have_final_energy_ratio) {
                final_energy_ratio_dec = arb_mid_to_decimal_exact_trunc(e_ratio, args.depth);
                final_energy_drop_pct_dec = arb_mid_to_decimal_exact_trunc(e_drop_pct, args.depth);
                have_final_energy_ratio = true;
            }

            std::string ratio_out = fixed_trunc_str_decimal_exact(final_energy_ratio_dec, args.depth);
            std::string drop_out  = fixed_trunc_str_decimal_exact(final_energy_drop_pct_dec, args.depth);
            ratio_out = trim_trailing_zeros_fixed(ratio_out);
            drop_out  = trim_trailing_zeros_fixed(drop_out);

            std::cout << "FINAL_ENERGY_RATIO:    " << ratio_out << "\n";
            std::cout << "FINAL_ENERGY_DROP_PCT: " << drop_out << "\n";

            arb_clear(e_now);
            arb_clear(e_plus);
            arb_clear(diff);
            arb_clear(slope);
            arb_clear(raw_jump);
            arb_clear(abs_raw_jump);
            arb_clear(h_arb);
            arb_clear(cap_arb);
            arb_clear(safe_jump);

            arb_clear(e_after);
            arb_clear(e_ratio);
            arb_clear(e_drop_pct);
            arb_clear(one_arb);
            arb_clear(hundred_arb);
            arb_clear(tmp_arb);

            return 0;
        }
    }

    std::cout << std::string(220, '-') << "\n";

    // ---------------------------------------------------------
    // FORCE RAW OUTPUT (The "Unsafe" Read)
    // ---------------------------------------------------------
    ctx.compute(t_work, sigma);

    // RAW: might be scientific; this is your unsafe read. Keep it as-is, but also provide fixed adaptive.
    std::string raw_energy = arb_mid_to_string_sig(ctx.abs_z, 3000);

    // Also provide a fixed adaptive version (DISPLAY ONLY) for readability and to avoid "0" artifacts.
    std::string raw_energy_fixed = arb_mid_to_fixed_decimal_adaptive_trim(ctx.abs_z, 120, args.depth);

    std::cout << "FINAL SINK ENERGY (RAW):       " << raw_energy << "\n";
    std::cout << "FINAL SINK ENERGY (RAW,FIXED): " << raw_energy_fixed << "\n";
    std::cout << "FINAL T VALUE:                 " << fixed_trunc_str_decimal_exact(t_work, work_frac_digits) << "\n";

    // Also print slope at whatever depth we ended with (strict fixed).
    Decimal end_slope_dec = arb_mid_to_decimal_exact_trunc(slope, work_frac_digits);
    std::cout << "FINAL_SLOPE(mid):              " << fixed_trunc_str_decimal_exact(end_slope_dec, work_frac_digits) << "\n";

    // Also print energy ratio/drop at whatever depth we ended with (fixed to depth here), trimmed.
    if (!have_final_energy_ratio) {
        final_energy_ratio_dec = arb_mid_to_decimal_exact_trunc(e_ratio, args.depth);
        final_energy_drop_pct_dec = arb_mid_to_decimal_exact_trunc(e_drop_pct, args.depth);
        have_final_energy_ratio = true;
    }

    std::string ratio_out = fixed_trunc_str_decimal_exact(final_energy_ratio_dec, args.depth);
    std::string drop_out  = fixed_trunc_str_decimal_exact(final_energy_drop_pct_dec, args.depth);
    ratio_out = trim_trailing_zeros_fixed(ratio_out);
    drop_out  = trim_trailing_zeros_fixed(drop_out);

    std::cout << "FINAL_ENERGY_RATIO:            " << ratio_out << "\n";
    std::cout << "FINAL_ENERGY_DROP_PCT:         " << drop_out << "\n";

    arb_clear(e_now);
    arb_clear(e_plus);
    arb_clear(diff);
    arb_clear(slope);
    arb_clear(raw_jump);
    arb_clear(abs_raw_jump);
    arb_clear(h_arb);
    arb_clear(cap_arb);
    arb_clear(safe_jump);

    arb_clear(e_after);
    arb_clear(e_ratio);
    arb_clear(e_drop_pct);
    arb_clear(one_arb);
    arb_clear(hundred_arb);
    arb_clear(tmp_arb);

    return 0;
}
