// nsniperv3.cpp
//
// Closure-stable Johansson-digit matching version (diagnostics stripped).
//
// Your request (implemented):
// - REMOVE the extra zeta call per step that is only used for diagnostics/progress.
//   Specifically, we NO LONGER do ctx.compute(t_work) AFTER updating t_work on each step.
//   Per iteration we now do ONLY:
//     1) E_now   = |zeta(0.5 + i*t_work)|
//     2) E_plus  = |zeta(0.5 + i*(t_work + h))|
//   These two are required to compute slope and jump.
//
// - Any other expensive calls used only for diagnostics were removed:
//     * No per-step ctx.compute(t_next) for progress energy
//     * No per-step ratio/drop computation
//   (Final ratio/drop is still produced, using the final sink energy vs the last pre-step energy.)
//
// - Still show a simple per-step progress line (cheap):
//     * prints one short line per step
//     * uses magnitude estimate of PRE-UPDATE energy E_now (already computed)
//     * uses magnitude estimate of jump (already computed)
//     * prints slope sign only
//     * prints frac after update (next_work_frac) so you see digit growth
//
// - Final report remains full true output:
//     * FINAL sink energy [mid +/- rad] at final t
//     * FINAL sink energy fixed adaptive
//     * FINAL T VALUE strict fixed-width at args.depth decimals (truncate only)
//     * FINAL_SLOPE(mid) strict fixed-width at args.depth decimals (truncate only)
//     * FINAL_ENERGY_RATIO / FINAL_ENERGY_DROP_PCT computed from final sink energy vs last pre-step energy
//
// NEW (your latest request):
// - CAP logic remains unchanged.
// - When CAP happens, set a flag "safety_after_cap".
// - From that point forward (for this run / sink), if decimals permit, reduce FIT digit growth by 1
//   additional digit: step_size_eff = step_size_capped - 2 (instead of -1), clamped to >= 1.
//   This makes jumps slightly less aggressive after the first CAP.
//
// Build:
//   g++ -O2 -std=c++17 nsniperv3.cpp -lflint -lgmp -lmpfr -o nsniperv3

#include <iostream>
#include <string>
#include <iomanip>
#include <cstdlib>
#include <cctype>
#include <cstdio>
#include <climits>

#include "decimal.hpp"
#include "arb_interface.hpp"

#include <flint/arf.h>
#include <flint/fmpz.h>

// -----------------------------
// Helpers (state-safe)
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
// DISPLAY helpers (final-report only)
// -----------------------------

// Trim trailing zeros after decimal point; remove '.' if it becomes the last char.
static std::string trim_trailing_zeros_fixed(std::string s) {
    if (s.find('e') != std::string::npos || s.find('E') != std::string::npos) return s;

    size_t dot = s.find('.');
    if (dot == std::string::npos) return s;

    while (!s.empty() && s.back() == '0') s.pop_back();
    if (!s.empty() && s.back() == '.') s.pop_back();

    if (s == "-0") s = "0";
    return s;
}

// Convert scientific string like "-3.82e-40" to fixed decimal (no exponent).
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

// Prints: LABEL: [mid +/- rad]
static void print_arb_bracket(const char* label, const arb_t x, slong digits) {
    std::cout << label << ": ";
    arb_printn(x, digits, 0);
    std::cout << "\n";
}

// Arb midpoint -> fixed decimal, adaptive enough to avoid printing as 0.
// DISPLAY ONLY.
static std::string arb_mid_to_fixed_decimal_adaptive_trim(const arb_t x, int sig_digits, int max_decimals) {
    const arf_struct *mid = arb_midref(x);

    if (arf_is_zero(mid)) return "0";

    // Attempt render
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

    // Estimate exponent, allow more decimals
    fmpz_t m_f;
    fmpz_init(m_f);
    slong e2 = 0;
    arf_get_fmpz_2exp(m_f, &e2, mid);

    slong m_bits = fmpz_bits(m_f);
    if (m_bits < 1) m_bits = 1;

    const long double LOG10_2 = 0.30102999566398119521L;
    long double log10_m_lo = (long double)(m_bits - 1) * LOG10_2;
    long double log10_mid  = log10_m_lo + (long double)e2 * LOG10_2;

    long long exp10 = (long long)(log10_mid >= 0 ? (long long)log10_mid : (long long)(log10_mid - 1.0L));
    long long needed_decimals = 0;
    if (exp10 < 0) needed_decimals = (-exp10) + 8;

    fmpz_clear(m_f);

    long long use_decimals_ll = (long long)max_decimals;
    if (needed_decimals > use_decimals_ll) use_decimals_ll = needed_decimals;

    const long long HARD_CAP = (long long)max_decimals + 1000;
    if (use_decimals_ll > HARD_CAP) use_decimals_ll = HARD_CAP;

    int use_decimals = (use_decimals_ll < 0) ? 0 : (int)use_decimals_ll;

    char *c2 = arf_get_str(mid, sig_digits * 2);
    if (!c2) {
        std::cerr << "ERROR: arf_get_str failed\n";
        std::exit(1);
    }
    std::string sci2(c2);
    flint_free(c2);

    std::string fixed2 = sci_to_fixed_noexp_trunc(sci2, use_decimals);
    fixed2 = trim_trailing_zeros_fixed(fixed2);

    if (fixed2 == "0") return sci2;
    return fixed2;
}

// -----------------------------
// Hot-loop progress helpers (NO decimal expansion / no extra zeta calls)
// -----------------------------

// Approx floor(log10(|mid|)) using (bitlen(m)-1 + e2) * log10(2).
static long long approx_floor_log10_abs_arf_mid(const arf_struct *mid) {
    if (arf_is_zero(mid)) return LLONG_MIN;

    fmpz_t m_f;
    fmpz_init(m_f);
    slong e2 = 0;
    arf_get_fmpz_2exp(m_f, &e2, mid);

    slong bits = fmpz_bits(m_f);
    if (bits < 1) bits = 1;
    fmpz_clear(m_f);

    const long double LOG10_2 = 0.30102999566398119521L;
    long double log10_est = ((long double)(bits - 1) + (long double)e2) * LOG10_2;

    long long k = (long long)log10_est;
    if (log10_est < 0 && (long double)k != log10_est) k -= 1;
    return k;
}

static long long approx_floor_log10_abs_arb_mid(const arb_t x) {
    return approx_floor_log10_abs_arf_mid(arb_midref(x));
}

static char slope_sign_char(const arb_t slope) {
    const arf_struct *m = arb_midref(slope);
    if (arf_is_zero(m)) return '0';
    return arf_sgn(m) < 0 ? '-' : '+';
}

static const char* dir_short(const std::string& dir_effective) {
    if (dir_effective == "down") return "D";
    if (dir_effective == "up") return "U";
    return "A";
}

// One short line per step. Energy token is PRE-UPDATE energy (E_now), already computed.
static void print_step_line(
    int step_idx,
    int frac_after_update,
    int step_size_eff,
    int probe_depth,
    bool is_cap,
    const std::string& dir_effective,
    const arb_t e_now,
    const arb_t safe_jump,
    const arb_t slope
) {
    char buf[256];

    const char *mode = is_cap ? "CAP" : "FIT";
    const char *dch  = dir_short(dir_effective);

    long long logE = approx_floor_log10_abs_arb_mid(e_now);
    long long logJ = approx_floor_log10_abs_arb_mid(safe_jump);
    char sgn = slope_sign_char(slope);

    char eTok[32], jTok[32];
    if (logE == LLONG_MIN) std::snprintf(eTok, sizeof(eTok), "0");
    else std::snprintf(eTok, sizeof(eTok), "1e%lld", logE);

    if (logJ == LLONG_MIN) std::snprintf(jTok, sizeof(jTok), "0");
    else std::snprintf(jTok, sizeof(jTok), "1e%lld", logJ);

    int n = std::snprintf(
        buf, sizeof(buf),
        "%4d frac=%-6d step=%-5d probe=%-5d %-3s dir=%s  E(pre)~%-10s  J~%-10s  S=%c\n",
        step_idx, frac_after_update, step_size_eff, probe_depth, mode, dch, eTok, jTok, sgn
    );
    if (n < 0) return;
    if ((size_t)n >= sizeof(buf)) n = (int)sizeof(buf) - 1;

    std::fwrite(buf, 1, (size_t)n, stdout);
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
    bool quiet = false;
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
        } else if (k == "--quiet") {
            a.quiet = true;
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
    Decimal t_work = decimal_from_string(t_string);

    int work_frac_digits = count_decimals_str(t_string);

    const int internal_max_frac = args.depth + args.overshoot_frac;

    int start_int_digits = count_int_digits(t_work);
    int compute_dps = (start_int_digits + internal_max_frac) + args.max_step + 256;
    slong prec_bits = dps_to_bits(compute_dps);

    ArbZetaContext ctx(prec_bits);
    const Decimal sigma = Decimal(1, 2);

    std::string dir_effective = args.dir;
    bool lock_found = false;

    // NEW: Once CAP triggers, permanently play it safer on digit growth for this run.
    bool safety_after_cap = false;

    // Cache final slope(mid) as arb (avoid Decimal conversion in hot loop).
    arb_t final_slope_arb;
    arb_init(final_slope_arb);
    bool have_final_slope = false;

    // Cache last pre-step energy (E_now) for final ratio/drop.
    arb_t last_e_pre;
    arb_init(last_e_pre);
    bool have_last_e_pre = false;

    std::cout << "DEPTH (exact frac digits): " << args.depth << "\n";
    std::cout << "OVERSHOOT_FRAC:            " << args.overshoot_frac << "\n";
    std::cout << "INTERNAL_MAX_FRAC:         " << internal_max_frac << " (depth + overshoot)\n";
    std::cout << "MAX_STEP:                  " << args.max_step << "\n";
    std::cout << "START_INT_DIGITS:          " << start_int_digits << "\n";
    std::cout << "COMPUTE_DPS:               " << compute_dps << "\n";
    std::cout << "ARB PREC (bits):           " << prec_bits << "\n\n";

    if (!args.quiet) {
        std::cout << "STEP TRACE (cheap): step, frac(after), step_size, probe_depth, CAP/FIT, dir, log10(|E_pre|), log10(|J|), slope sign\n";
    }

    int total_steps = 0;

    arb_t e_now, e_plus, diff, slope, raw_jump, abs_raw_jump;
    arb_t h_arb, cap_arb, safe_jump;
    arb_t e_ratio, e_drop_pct, one_arb, hundred_arb, tmp_arb; // final only
    arb_init(e_now);
    arb_init(e_plus);
    arb_init(diff);
    arb_init(slope);
    arb_init(raw_jump);
    arb_init(abs_raw_jump);
    arb_init(h_arb);
    arb_init(cap_arb);
    arb_init(safe_jump);

    arb_init(e_ratio);
    arb_init(e_drop_pct);
    arb_init(one_arb);
    arb_init(hundred_arb);
    arb_init(tmp_arb);
    arb_one(one_arb);
    arb_set_si(hundred_arb, 100);

    while (total_steps < args.max_steps) {
        total_steps += 1;

        int step_size_raw = (work_frac_digits < 5) ? 1 : work_frac_digits;

        int step_size_capped = step_size_raw;
        if (step_size_capped > args.max_step) step_size_capped = args.max_step;

        int step_size_eff;
        if (work_frac_digits < 5) {
            step_size_eff = 1;
        } else {
            int subtract = 1;

            // After CAP triggers once, play it safer: subtract 2 if it makes sense.
            // "if decimals available" -> require step_size_capped >= 3 so (capped - 2) stays >= 1.
            if (safety_after_cap && step_size_capped >= 3) {
                subtract = 2;
            }

            step_size_eff = step_size_capped - subtract;
            if (step_size_eff < 1) step_size_eff = 1;
        }

        int target_frac = work_frac_digits + step_size_eff;
        if (target_frac > internal_max_frac) target_frac = internal_max_frac;

        int probe_depth;
        if (step_size_capped == args.max_step) probe_depth = args.depth + 64;
        else probe_depth = target_frac + 5;

        if (probe_depth > internal_max_frac) probe_depth = internal_max_frac;
        if (probe_depth < 1) probe_depth = 1;

        Decimal h_dec(1, pow10_z(probe_depth));
        ctx.arb_set_decimal_safe(h_arb, h_dec);

        // 1) E_now at current t_work
        ctx.compute(t_work, sigma);
        arb_set(e_now, ctx.abs_z);

        // cache last pre-step energy for final ratio/drop
        arb_set(last_e_pre, e_now);
        have_last_e_pre = true;

        // 2) E_plus at t_work + h
        Decimal t_plus = t_work + h_dec;
        ctx.compute(t_plus, sigma);
        arb_set(e_plus, ctx.abs_z);

        // slope = (E_plus - E_now)/h
        arb_sub(diff, e_plus, e_now, prec_bits);
        if (arb_is_zero(diff)) {
            arb_set_str(slope, "1e-300", prec_bits);
        } else {
            arb_div(slope, diff, h_arb, prec_bits);
        }

        // raw_jump = E_now / slope
        arb_div(raw_jump, e_now, slope, prec_bits);

        // cap = 10^-work_frac_digits (or 0.1 at integer start)
        Decimal cap_dec(1, pow10_z(work_frac_digits));
        if (work_frac_digits == 0) cap_dec = Decimal(1, 10);
        ctx.arb_set_decimal_safe(cap_arb, cap_dec);

        arb_abs(abs_raw_jump, raw_jump);
        bool is_cap = arb_gt(abs_raw_jump, cap_arb);

        // NEW: CAP still functions normally, but sets the safety flag once triggered.
        if (is_cap) {
            safety_after_cap = true;
        }

        if (is_cap && !lock_found) {
            lock_found = true;
            dir_effective = "auto";
        }

        int next_work_frac = work_frac_digits;

        if (is_cap) {
            if (work_frac_digits == 0) next_work_frac = 1;
            else next_work_frac = work_frac_digits;

            if (dir_effective == "down") {
                arb_set(safe_jump, cap_arb);
            } else if (dir_effective == "up") {
                arb_neg(safe_jump, cap_arb);
            } else {
                // CAP auto chooses direction via 2 extra zeta calls (this IS used to decide next update).
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
                    arb_neg(safe_jump, cap_arb);
                } else {
                    arb_set(safe_jump, cap_arb);
                }

                arb_clear(e_cap_plus);
                arb_clear(e_cap_minus);
            }
        } else {
            next_work_frac = target_frac;
            arb_set(safe_jump, raw_jump);

            if (!lock_found) {
                if (dir_effective == "down") {
                    arb_abs(safe_jump, raw_jump);
                } else if (dir_effective == "up") {
                    arb_abs(safe_jump, raw_jump);
                    arb_neg(safe_jump, safe_jump);
                }
            }
        }

        // cache slope(mid) arb when we first hit depth-range
        if (!have_final_slope && next_work_frac >= args.depth) {
            arb_set(final_slope_arb, slope);
            have_final_slope = true;
        }

        // ---- APPLY UPDATE IN DECIMAL SPACE ----
        Decimal t_next_dec = t_work;

        if (is_cap) {
            if (arb_is_negative(safe_jump)) t_next_dec = t_work + cap_dec;
            else t_next_dec = t_work - cap_dec;
        } else {
            int jump_trunc_decimals = next_work_frac + 256;
            Decimal jump_dec = arb_mid_to_decimal_exact_trunc(safe_jump, jump_trunc_decimals);
            t_next_dec = t_work - jump_dec;
        }

        if (next_work_frac > internal_max_frac) next_work_frac = internal_max_frac;
        t_work = strict_truncate_toward0(t_next_dec, next_work_frac);

        // ---- cheap per-step trace (uses E_now, no extra zeta compute) ----
        if (!args.quiet) {
            print_step_line(total_steps, next_work_frac, step_size_eff, probe_depth, is_cap, dir_effective, e_now, safe_jump, slope);
        }

        work_frac_digits = next_work_frac;

        // Termination: when reached depth, compute final energy once at final t (required).
        if (work_frac_digits >= args.depth) {
            Decimal t_final = t_work;
            if (work_frac_digits > args.depth) {
                t_final = strict_truncate_toward0(t_work, args.depth);
            }

            ctx.compute(t_final, sigma); // FINAL sink energy at t_final in ctx.abs_z

            // Final ratio/drop computed vs last pre-step energy (E_pre)
            if (!have_last_e_pre || arb_is_zero(last_e_pre)) {
                arb_zero(e_ratio);
                arb_zero(e_drop_pct);
            } else {
                arb_div(e_ratio, ctx.abs_z, last_e_pre, prec_bits);
                arb_sub(tmp_arb, one_arb, e_ratio, prec_bits);
                arb_mul(e_drop_pct, tmp_arb, hundred_arb, prec_bits);
            }

            std::cout << ">>> REACHED TARGET DEPTH (SETTLED) <<<\n";
            std::cout << std::string(220, '-') << "\n";

            print_arb_bracket("FINAL SINK ENERGY (SCI, MID+/-RAD)", ctx.abs_z, 50);

            std::string final_energy_fixed = arb_mid_to_fixed_decimal_adaptive_trim(ctx.abs_z, 120, args.depth + 256);
            std::cout << "FINAL SINK ENERGY (MID, FIXED): " << final_energy_fixed << "\n";

            std::cout << "FINAL T VALUE:     " << fixed_trunc_str_decimal_exact(t_final, args.depth) << "\n";

            if (!have_final_slope) {
                arb_set(final_slope_arb, slope);
                have_final_slope = true;
            }
            Decimal final_slope_mid_dec = arb_mid_to_decimal_exact_trunc(final_slope_arb, args.depth);
            std::cout << "FINAL_SLOPE(mid):  " << fixed_trunc_str_decimal_exact(final_slope_mid_dec, args.depth) << "\n";

            Decimal final_energy_ratio_dec = arb_mid_to_decimal_exact_trunc(e_ratio, args.depth);
            Decimal final_energy_drop_pct_dec = arb_mid_to_decimal_exact_trunc(e_drop_pct, args.depth);

            std::string ratio_out = fixed_trunc_str_decimal_exact(final_energy_ratio_dec, args.depth);
            std::string drop_out  = fixed_trunc_str_decimal_exact(final_energy_drop_pct_dec, args.depth);
            ratio_out = trim_trailing_zeros_fixed(ratio_out);
            drop_out  = trim_trailing_zeros_fixed(drop_out);

            std::cout << "FINAL_ENERGY_RATIO:    " << ratio_out << "\n";
            std::cout << "FINAL_ENERGY_DROP_PCT: " << drop_out << "\n";

            // clear
            arb_clear(e_now);
            arb_clear(e_plus);
            arb_clear(diff);
            arb_clear(slope);
            arb_clear(raw_jump);
            arb_clear(abs_raw_jump);
            arb_clear(h_arb);
            arb_clear(cap_arb);
            arb_clear(safe_jump);

            arb_clear(e_ratio);
            arb_clear(e_drop_pct);
            arb_clear(one_arb);
            arb_clear(hundred_arb);
            arb_clear(tmp_arb);

            arb_clear(final_slope_arb);
            arb_clear(last_e_pre);

            return 0;
        }
    }

    std::cout << std::string(220, '-') << "\n";

    // ---------------------------------------------------------
    // If max_steps reached without depth, do a final compute once.
    // ---------------------------------------------------------
    ctx.compute(t_work, sigma);

    std::cout << "FINAL SINK ENERGY (RAW):       ";
    arb_printn(ctx.abs_z, 3000, 0);
    std::cout << "\n";

    std::string raw_energy_fixed = arb_mid_to_fixed_decimal_adaptive_trim(ctx.abs_z, 120, args.depth);
    std::cout << "FINAL SINK ENERGY (RAW,FIXED): " << raw_energy_fixed << "\n";
    std::cout << "FINAL T VALUE:                 " << fixed_trunc_str_decimal_exact(t_work, work_frac_digits) << "\n";

    Decimal end_slope_dec = arb_mid_to_decimal_exact_trunc(slope, work_frac_digits);
    std::cout << "FINAL_SLOPE(mid):              " << fixed_trunc_str_decimal_exact(end_slope_dec, work_frac_digits) << "\n";

    if (!have_last_e_pre || arb_is_zero(last_e_pre)) {
        arb_zero(e_ratio);
        arb_zero(e_drop_pct);
    } else {
        arb_div(e_ratio, ctx.abs_z, last_e_pre, prec_bits);
        arb_sub(tmp_arb, one_arb, e_ratio, prec_bits);
        arb_mul(e_drop_pct, tmp_arb, hundred_arb, prec_bits);
    }

    Decimal final_energy_ratio_dec = arb_mid_to_decimal_exact_trunc(e_ratio, args.depth);
    Decimal final_energy_drop_pct_dec = arb_mid_to_decimal_exact_trunc(e_drop_pct, args.depth);

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

    arb_clear(e_ratio);
    arb_clear(e_drop_pct);
    arb_clear(one_arb);
    arb_clear(hundred_arb);
    arb_clear(tmp_arb);

    arb_clear(final_slope_arb);
    arb_clear(last_e_pre);

    return 0;
}
