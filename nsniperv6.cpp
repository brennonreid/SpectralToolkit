// nsniperv5.cpp
//
// Single-call-per-iteration HardyZ sniper.
//
// Hot-loop:
// - Exactly ONE HardyZ(+deriv) evaluation per iteration.
// - No probe evaluation.
// - Slope uses derivative: d|Z|/dt = sign(Z)*Z'(t).
//
// Outputs:
// - If --out is provided: write 3 files (out, out.steps.txt, out.t.txt unless overridden).
// - If --out is omitted: NO files are written (stdout only).
//
// T-only file content:
// - Exactly the final T value digits, followed by '\n' only.

#include <iostream>
#include <string>
#include <cstdlib>
#include <cctype>
#include <cstdio>
#include <climits>
#include <cstdint>
#include <unordered_map>
#include <cstring>

#include "decimal.hpp"
#include "arb_interface.hpp"

#include <flint/arf.h>
#include <flint/fmpz.h>

// -----------------------------
// pow10 cache (mpz)
// -----------------------------
struct Pow10Cache {
    std::unordered_map<int, mpz_class> cache;

    const mpz_class& get(int n) {
        if (n < 0) {
            std::cerr << "ERROR: pow10 exponent must be >= 0\n";
            std::exit(1);
        }
        auto it = cache.find(n);
        if (it != cache.end()) return it->second;

        mpz_class v = pow10_z(n);
        auto res = cache.emplace(n, v);
        return res.first->second;
    }
};

// -----------------------------
// Decimal helpers
// -----------------------------
static Decimal strict_truncate_toward0_cached(const Decimal& x, int decimals, Pow10Cache& p10) {
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

    const mpz_class& scale = p10.get(decimals);
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

static inline Decimal t_can_cached(const Decimal& t, int frac_digits, Pow10Cache& p10) {
    return strict_truncate_toward0_cached(t, frac_digits, p10);
}

// Fixed-point string with exactly "decimals" digits after '.', truncating toward zero.
static std::string fixed_trunc_str_decimal_exact_cached(const Decimal& val, int decimals, Pow10Cache& p10) {
    if (decimals < 0) {
        std::cerr << "ERROR: decimals must be >= 0\n";
        std::exit(1);
    }

    Decimal x = val;
    bool neg = (x < Decimal(0));
    if (neg) x = -x;

    if (decimals == 0) {
        Decimal t = strict_truncate_toward0_cached(x, 0, p10);
        std::string s = decimal_to_string(t, 0);
        if (neg && s != "0") return "-" + s;
        return s;
    }

    const mpz_class& scale = p10.get(decimals);

    Decimal scaled = strict_truncate_toward0_cached(x * Decimal(scale), 0, p10);
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

static int count_int_digits_cached(const Decimal& x, Pow10Cache& p10) {
    Decimal t_int = strict_truncate_toward0_cached(x, 0, p10);
    std::string is = decimal_to_string(t_int, 0);
    if (!is.empty() && is[0] == '-') is.erase(0, 1);
    int d = (int)is.size();
    if (d < 1) d = 1;
    return d;
}

static slong dps_to_bits(int dps) {
    long double bits = (long double)dps * 3.32192809488736234787L;
    bits += 64.0L;
    if (bits < 128.0L) bits = 128.0L;
    return (slong)bits;
}

static slong clamp_step_bits(slong b) {
    if (b < 256) b = 256;
    return b;
}

// Convert arf midpoint -> Decimal, then truncate.
static Decimal arf_mid_to_decimal_exact_trunc_pow10_cached(const arf_struct* mid, int decimals, Pow10Cache& p10) {
    if (decimals < 0) {
        std::cerr << "ERROR: decimals must be >= 0\n";
        std::exit(1);
    }
    if (arf_is_zero(mid)) return Decimal(0);

    fmpz_t m_f;
    fmpz_init(m_f);

    slong e2 = 0;
    arf_get_fmpz_2exp(m_f, &e2, mid);

    char* m_str = fmpz_get_str(NULL, 10, m_f);
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

    return strict_truncate_toward0_cached(r, decimals, p10);
}

static Decimal arb_mid_to_decimal_exact_trunc_cached(const arb_t x, int decimals, Pow10Cache& p10) {
    return arf_mid_to_decimal_exact_trunc_pow10_cached(arb_midref(x), decimals, p10);
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

    bool out_enabled = false;
    std::string out_path = "";
    std::string out_steps_path = "";
    std::string out_t_path = "";
};

static Args parse_args(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; i++) {
        std::string k = argv[i];

        auto require_value = [&](const std::string& flag) {
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
        } else if (k == "--out") {
            require_value(k);
            a.out_enabled = true;
            a.out_path = argv[++i];
            if (a.out_path.empty()) {
                std::cerr << "ERROR: --out path must be non-empty\n";
                std::exit(1);
            }
        } else if (k == "--out_steps") {
            require_value(k);
            a.out_steps_path = argv[++i];
            if (a.out_steps_path.empty()) {
                std::cerr << "ERROR: --out_steps path must be non-empty\n";
                std::exit(1);
            }
        } else if (k == "--out_t") {
            require_value(k);
            a.out_t_path = argv[++i];
            if (a.out_t_path.empty()) {
                std::cerr << "ERROR: --out_t path must be non-empty\n";
                std::exit(1);
            }
        } else {
            std::cerr << "ERROR: Unknown argument: " << k << "\n";
            std::exit(1);
        }
    }

    if (!a.out_enabled) {
        if (!a.out_steps_path.empty() || !a.out_t_path.empty()) {
            std::cerr << "ERROR: --out_steps/--out_t require --out\n";
            std::exit(1);
        }
    } else {
        if (a.out_steps_path.empty()) a.out_steps_path = a.out_path + ".steps.txt";
        if (a.out_t_path.empty()) a.out_t_path = a.out_path + ".t.txt";
    }

    return a;
}

static inline std::string trim_local(std::string s) {
    while (!s.empty() && std::isspace((unsigned char)s.front())) s.erase(s.begin());
    while (!s.empty() && std::isspace((unsigned char)s.back())) s.pop_back();
    return s;
}

// log10(|mid(x)|) -> out
static void arb_log10_abs_mid(arb_t out, const arb_t x, const arb_t ln10, slong prec_bits) {
    arb_t t;
    arb_init(t);
    arb_abs(t, x);
    arb_log(t, t, prec_bits);
    arb_div(out, t, ln10, prec_bits);
    arb_clear(t);
}

static void print_arb_sci(FILE* f_or_null, const arb_t x, slong digits) {
    if (f_or_null) arb_fprintn(f_or_null, x, digits, 0);
    else arb_printn(x, digits, 0);
}

static void emit_divider(FILE* f_or_null) {
    const char* line = "----------------------------------------------------------------\n";
    if (f_or_null) std::fwrite(line, 1, std::strlen(line), f_or_null);
    else std::fwrite(line, 1, std::strlen(line), stdout);
}

static void emit_kv_int(FILE* f_or_null, const char* key, long long v) {
    if (f_or_null) std::fprintf(f_or_null, "%s: %lld\n", key, v);
    else std::printf("%s: %lld\n", key, v);
}

static void emit_kv_str(FILE* f_or_null, const char* key, const std::string& v) {
    if (f_or_null) std::fprintf(f_or_null, "%s: %s\n", key, v.c_str());
    else std::printf("%s: %s\n", key, v.c_str());
}

static void emit_kv_arb(FILE* f_or_null, const char* key, const arb_t x, slong digits) {
    if (f_or_null) std::fprintf(f_or_null, "%s: ", key);
    else std::printf("%s: ", key);
    print_arb_sci(f_or_null, x, digits);
    if (f_or_null) std::fprintf(f_or_null, "\n");
    else std::printf("\n");
}

static void emit_step_block(
    FILE* f_or_null,
    int iter,
    int frac_digits,
    const char* type,
    int step_size_eff,
    const std::string& gap_str,
    const arb_t jump_mid,
    const arb_t j_log10,
    const arb_t slope,
    const arb_t energy,
    const arb_t e_log10,
    const std::string& t_str,
    slong digits
) {
    emit_divider(f_or_null);
    emit_kv_int(f_or_null, "ITER", iter);
    emit_kv_int(f_or_null, "FRAC", frac_digits);
    emit_kv_str(f_or_null, "TYPE", type);
    emit_kv_int(f_or_null, "STEP_SIZE", step_size_eff);
    emit_kv_str(f_or_null, "GAP", gap_str);
    emit_kv_arb(f_or_null, "JUMP(mid)", jump_mid, digits);
    emit_kv_arb(f_or_null, "J_LOG10", j_log10, digits);
    emit_kv_arb(f_or_null, "SLOPE(mid)", slope, digits);
    emit_kv_arb(f_or_null, "ENERGY", energy, digits);
    emit_kv_arb(f_or_null, "E_LOG10", e_log10, digits);
    emit_kv_str(f_or_null, "T", t_str);
}

static void emit_done_block(
    FILE* f_or_null,
    int iter,
    int frac_digits,
    const char* type,
    int step_size_eff,
    const std::string& gap_str,
    const arb_t jump_mid,
    const arb_t j_log10,
    const arb_t slope,
    const arb_t energy,
    const arb_t e_log10,
    const std::string& t_str,
    const std::string& final_gap_str,
    slong digits
) {
    emit_divider(f_or_null);
    emit_kv_str(f_or_null, "DONE", "1");
    emit_kv_int(f_or_null, "ITER", iter);
    emit_kv_int(f_or_null, "FRAC", frac_digits);
    emit_kv_str(f_or_null, "TYPE", type);
    emit_kv_int(f_or_null, "STEP_SIZE", step_size_eff);
    emit_kv_str(f_or_null, "GAP", gap_str);
    emit_kv_arb(f_or_null, "JUMP(mid)", jump_mid, digits);
    emit_kv_arb(f_or_null, "J_LOG10", j_log10, digits);
    emit_kv_arb(f_or_null, "SLOPE(mid)", slope, digits);
    emit_kv_arb(f_or_null, "ENERGY", energy, digits);
    emit_kv_arb(f_or_null, "E_LOG10", e_log10, digits);
    emit_kv_str(f_or_null, "T", t_str);
    emit_kv_str(f_or_null, "FINAL_GAP", final_gap_str);
    emit_divider(f_or_null);
}

// -----------------------------
// Main
// -----------------------------
int main(int argc, char** argv) {
    Args args = parse_args(argc, argv);
    Pow10Cache p10;

    FILE* fout = NULL;
    FILE* fsteps = NULL;
    FILE* ft = NULL;

    if (args.out_enabled) {
        fout = std::fopen(args.out_path.c_str(), "wb");
        if (!fout) { std::perror("ERROR: fopen(--out) failed"); std::exit(1); }

        fsteps = std::fopen(args.out_steps_path.c_str(), "wb");
        if (!fsteps) { std::perror("ERROR: fopen(--out_steps) failed"); std::fclose(fout); std::exit(1); }

        ft = std::fopen(args.out_t_path.c_str(), "wb");
        if (!ft) { std::perror("ERROR: fopen(--out_t) failed"); std::fclose(fsteps); std::fclose(fout); std::exit(1); }
    }

    std::string t_string = trim_local(args.t);
    Decimal t_work = decimal_from_string(t_string);

    int work_frac_digits = count_decimals_str(t_string);
    t_work = t_can_cached(t_work, work_frac_digits, p10);

    const int internal_max_frac = args.depth + args.overshoot_frac;

    int start_int_digits = count_int_digits_cached(t_work, p10);
    int compute_dps = (start_int_digits + internal_max_frac) + args.max_step + 256;
    slong prec_bits = dps_to_bits(compute_dps);

    ArbZetaContext ctx(prec_bits);

    std::string dir_effective = args.dir;
    bool safety_after_cap = false;

    if (args.out_enabled) {
        std::fprintf(fout, "nsniperv5 HardyZ single-call hot-loop log\n");
        std::fprintf(fout, "t=%s\n", args.t.c_str());
        std::fprintf(fout, "depth=%d\n", args.depth);
        std::fprintf(fout, "max_steps=%d\n", args.max_steps);
        std::fprintf(fout, "overshoot_frac=%d\n", args.overshoot_frac);
        std::fprintf(fout, "max_step=%d\n", args.max_step);
        std::fprintf(fout, "dir=%s\n", args.dir.c_str());
        std::fprintf(fout, "compute_dps=%d\n", compute_dps);
        std::fprintf(fout, "arb_prec_bits=%lld\n", (long long)prec_bits);
        std::fprintf(fout, "----------------------------------------\n");
        std::fflush(fout);
    }

    int total_steps = 0;

    arb_t e_now, slope, raw_jump, abs_raw_jump;
    arb_t cap_arb, safe_jump;

    arb_init(e_now);
    arb_init(slope);
    arb_init(raw_jump);
    arb_init(abs_raw_jump);
    arb_init(cap_arb);
    arb_init(safe_jump);

    // ln(10) helper
    arb_t ln10;
    arb_init(ln10);
    arb_log_ui(ln10, 10, prec_bits);

    // GAP this iteration (Decimal, as actually applied)
    Decimal gap_dec = Decimal(0);

    while (total_steps < args.max_steps) {
        total_steps += 1;

        // Step size schedule (doubling regime) and exact-landing clamp.
        int step_size_raw = (work_frac_digits < 5) ? 1 : work_frac_digits;

        int step_size_capped = step_size_raw;
        if (step_size_capped > args.max_step) step_size_capped = args.max_step;

        int step_size_eff;
        if (work_frac_digits < 5) {
            step_size_eff = 1;
        } else {
            int subtract = 1;
            if (safety_after_cap && step_size_capped >= 3) subtract = 2;
            step_size_eff = step_size_capped - subtract;
            if (step_size_eff < 1) step_size_eff = 1;
        }

        if (work_frac_digits < args.depth) {
            int remaining = args.depth - work_frac_digits;
            if (remaining < 1) remaining = 1;
            if (step_size_eff > remaining) step_size_eff = remaining;
        } else {
            step_size_eff = 1;
        }

        int target_frac = work_frac_digits + step_size_eff;
        if (target_frac > args.depth) target_frac = args.depth;

        // Precision driver (not a probe).
        int prec_frac = work_frac_digits;
        if (prec_frac < 6) prec_frac = 6;
        prec_frac += 5;
        if (prec_frac > internal_max_frac) prec_frac = internal_max_frac;
        if (prec_frac < 1) prec_frac = 1;

        int base_guard = 256;
        int cancel_guard = prec_frac + 64;
        int step_dps = (start_int_digits + prec_frac) + base_guard + cancel_guard;
        slong prec_step_bits = clamp_step_bits(dps_to_bits(step_dps));

        // Reporting digits (reporting-only)
        slong report_digits = (slong)step_dps;
        if (report_digits < 80) report_digits = 80;
        report_digits += 8;

        // Single HardyZ+derivative evaluation.
        ctx.compute_hardy_z_and_deriv_prec(t_work, prec_step_bits);

        // Energy = |Z(t)|.
        arb_set(e_now, ctx.abs_z);

        // Slope = d|Z|/dt = sign(Z)*Z'(t).
        arb_set(slope, ctx.hardy_zp);
        if (arb_is_negative(ctx.hardy_z)) {
            arb_neg(slope, slope);
        }
        if (arb_is_zero(slope)) {
            arb_set_str(slope, "1e-300", prec_step_bits);
        }

        // Newton jump.
        arb_div(raw_jump, e_now, slope, prec_step_bits);

        // CAP threshold.
        Decimal cap_dec;
        if (work_frac_digits == 0) cap_dec = Decimal(1, 10);
        else cap_dec = Decimal(1, p10.get(work_frac_digits));
        ctx.arb_set_decimal_exact_prec(cap_arb, cap_dec, prec_step_bits);

        arb_abs(abs_raw_jump, raw_jump);
        bool is_cap = arb_gt(abs_raw_jump, cap_arb);

        if (is_cap && work_frac_digits > 5) safety_after_cap = true;

        int next_work_frac = work_frac_digits;

        if (is_cap) {
            if (work_frac_digits == 0) next_work_frac = 1;
            else next_work_frac = work_frac_digits;

            if (dir_effective == "down") {
                arb_set(safe_jump, cap_arb);
            } else if (dir_effective == "up") {
                arb_neg(safe_jump, cap_arb);
            } else {
                // Auto CAP direction without extra eval: follow Newton direction from raw_jump sign.
                if (arb_is_negative(raw_jump)) arb_neg(safe_jump, cap_arb);
                else arb_set(safe_jump, cap_arb);
            }

            gap_dec = cap_dec;
        } else {
            next_work_frac = target_frac;
            arb_set(safe_jump, raw_jump);

            if (dir_effective == "down") {
                arb_abs(safe_jump, raw_jump);
            } else if (dir_effective == "up") {
                arb_abs(safe_jump, raw_jump);
                arb_neg(safe_jump, safe_jump);
            }
        }

        // Apply update (truncation enforces exact landing behavior).
        Decimal t_next_dec = t_work;

        if (is_cap) {
            if (arb_is_negative(safe_jump)) t_next_dec = t_work + cap_dec;
            else t_next_dec = t_work - cap_dec;
        } else {
            int jump_trunc_decimals = next_work_frac + 256;
            Decimal jump_dec = arb_mid_to_decimal_exact_trunc_cached(safe_jump, jump_trunc_decimals, p10);

            // GAP for FIT is abs(jump_dec) as actually applied
            gap_dec = (jump_dec < Decimal(0)) ? (-jump_dec) : jump_dec;

            t_next_dec = t_work - jump_dec;
        }

        if (next_work_frac > internal_max_frac) next_work_frac = internal_max_frac;
        if (next_work_frac > args.depth) next_work_frac = args.depth;

        t_work = strict_truncate_toward0_cached(t_next_dec, next_work_frac, p10);

        // Reporting fields for this iteration
        arb_t e_log10, j_log10;
        arb_init(e_log10);
        arb_init(j_log10);

        arb_log10_abs_mid(e_log10, e_now, ln10, prec_step_bits);
        arb_log10_abs_mid(j_log10, safe_jump, ln10, prec_step_bits);

        std::string t_str = fixed_trunc_str_decimal_exact_cached(t_work, next_work_frac, p10);
        std::string gap_str = fixed_trunc_str_decimal_exact_cached(gap_dec, next_work_frac, p10);
        const char* type = is_cap ? "CAP" : "FIT";

        // Finalization: if this iteration reaches depth, emit DONE once.
        if (next_work_frac >= args.depth) {
            std::string final_gap_str = fixed_trunc_str_decimal_exact_cached(gap_dec, args.depth, p10);

            if (!args.quiet) {
                emit_done_block(NULL, total_steps, next_work_frac, type, step_size_eff,
                                gap_str, safe_jump, j_log10, slope, e_now, e_log10, t_str,
                                final_gap_str, report_digits);
            }
            if (args.out_enabled) {
                emit_done_block(fsteps, total_steps, next_work_frac, type, step_size_eff,
                                gap_str, safe_jump, j_log10, slope, e_now, e_log10, t_str,
                                final_gap_str, report_digits);
                std::fflush(fsteps);

                // Optional human log file (still not a table)
                std::fprintf(fout, ">>> DONE <<<\n");
                std::fprintf(fout, "T: %s\n", t_str.c_str());
                std::fprintf(fout, "FINAL_GAP: %s\n", final_gap_str.c_str());
                std::fprintf(fout, "ENERGY: ");
                arb_fprintn(fout, e_now, report_digits, 0);
                std::fprintf(fout, "\n");
                std::fprintf(fout, "SLOPE(mid): ");
                arb_fprintn(fout, slope, report_digits, 0);
                std::fprintf(fout, "\n");
                std::fprintf(fout, "JUMP(mid): ");
                arb_fprintn(fout, safe_jump, report_digits, 0);
                std::fprintf(fout, "\n");
                std::fprintf(fout, "E_LOG10: ");
                arb_fprintn(fout, e_log10, report_digits, 0);
                std::fprintf(fout, "\n");
                std::fprintf(fout, "J_LOG10: ");
                arb_fprintn(fout, j_log10, report_digits, 0);
                std::fprintf(fout, "\n");
                std::fflush(fout);

                // T-only file: exactly digits + '\n'
                std::fwrite(t_str.data(), 1, t_str.size(), ft);
                std::fwrite("\n", 1, 1, ft);
                std::fflush(ft);
            }

            arb_clear(e_log10);
            arb_clear(j_log10);

            arb_clear(ln10);

            arb_clear(e_now);
            arb_clear(slope);
            arb_clear(raw_jump);
            arb_clear(abs_raw_jump);
            arb_clear(cap_arb);
            arb_clear(safe_jump);

            if (args.out_enabled) {
                std::fclose(ft);
                std::fclose(fsteps);
                std::fclose(fout);
            }
            return 0;
        }

        // Normal per-step block
        if (!args.quiet) {
            emit_step_block(NULL, total_steps, next_work_frac, type, step_size_eff,
                            gap_str, safe_jump, j_log10, slope, e_now, e_log10, t_str,
                            report_digits);
        }
        if (args.out_enabled) {
            emit_step_block(fsteps, total_steps, next_work_frac, type, step_size_eff,
                            gap_str, safe_jump, j_log10, slope, e_now, e_log10, t_str,
                            report_digits);
            if ((total_steps % 25) == 0) std::fflush(fsteps);
        }

        arb_clear(e_log10);
        arb_clear(j_log10);

        work_frac_digits = next_work_frac;
    }

    // Max steps reached: one call here is unavoidable.
    ctx.compute_hardy_z_and_deriv_prec(t_work, prec_bits);

    arb_t final_slope;
    arb_init(final_slope);
    arb_set(final_slope, ctx.hardy_zp);
    if (arb_is_negative(ctx.hardy_z)) arb_neg(final_slope, final_slope);

    std::string t_str = fixed_trunc_str_decimal_exact_cached(t_work, work_frac_digits, p10);
    std::string gap_str = fixed_trunc_str_decimal_exact_cached(gap_dec, work_frac_digits, p10);

    arb_t e_log10, j_log10;
    arb_init(e_log10);
    arb_init(j_log10);
    arb_log10_abs_mid(e_log10, ctx.abs_z, ln10, prec_bits);
    arb_log10_abs_mid(j_log10, safe_jump, ln10, prec_bits);

    slong report_digits = (slong)compute_dps + 16;
    if (report_digits < 80) report_digits = 80;

    if (!args.quiet) {
        emit_done_block(NULL, args.max_steps, work_frac_digits, "MAX_STEPS", 0,
                        gap_str, safe_jump, j_log10, final_slope, ctx.abs_z, e_log10, t_str,
                        gap_str, report_digits);
    }
    if (args.out_enabled) {
        emit_done_block(fsteps, args.max_steps, work_frac_digits, "MAX_STEPS", 0,
                        gap_str, safe_jump, j_log10, final_slope, ctx.abs_z, e_log10, t_str,
                        gap_str, report_digits);
        std::fflush(fsteps);

        std::fwrite(t_str.data(), 1, t_str.size(), ft);
        std::fwrite("\n", 1, 1, ft);
        std::fflush(ft);
    }

    arb_clear(e_log10);
    arb_clear(j_log10);

    arb_clear(final_slope);

    arb_clear(ln10);

    arb_clear(e_now);
    arb_clear(slope);
    arb_clear(raw_jump);
    arb_clear(abs_raw_jump);
    arb_clear(cap_arb);
    arb_clear(safe_jump);

    if (args.out_enabled) {
        std::fclose(ft);
        std::fclose(fsteps);
        std::fclose(fout);
    }

    return 0;
}
