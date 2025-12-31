// nsniperv4.cpp
//
// Minimal, exact-landing v4:
//
// - NO energy caching.
// - NO secant slope.
// - Keeps probe sizing logic unchanged.
// - Keeps CAP logic unchanged.
// - Adds "exact landing" step scheduling: doubling regime, but clamps the effective step
//   so we never overshoot args.depth.
//
// Output files:
//   --out <path>        : main log file (final results + key run metadata)
//   --out_steps <path>  : per-step trace file (one line per step)
//   --out_t <path>      : ONLY final T value (raw digits, one line)
//
// Defaults:
//   --out        defaults to "nsniperv4_output.txt"
//   --out_steps  defaults to "<out>.steps.txt"
//   --out_t      defaults to "<out>.t.txt"
//
// Hard-fail if any output file cannot be opened.
//
// Build:
//   g++ -O2 -std=c++17 nsniperv4.cpp -lflint -lgmp -lmpfr -o nsniperv4

#include <iostream>
#include <string>
#include <iomanip>
#include <cstdlib>
#include <cctype>
#include <cstdio>
#include <climits>
#include <chrono>
#include <cstdint>
#include <unordered_map>

#include "decimal.hpp"
#include "arb_interface.hpp"

#include <flint/arf.h>
#include <flint/fmpz.h>

// -----------------------------
// pow10 cache (mpz) - HOT
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
// Helpers (state-safe)
// -----------------------------
static Decimal strict_truncate_toward0_cached(const Decimal &x, int decimals, Pow10Cache& p10) {
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

static inline Decimal t_can_cached(const Decimal &t, int frac_digits, Pow10Cache& p10) {
    return strict_truncate_toward0_cached(t, frac_digits, p10);
}

static std::string fixed_trunc_str_decimal_exact_cached(const Decimal &val, int decimals, Pow10Cache& p10) {
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

static int count_int_digits_cached(const Decimal &x, Pow10Cache& p10) {
    Decimal t_int = strict_truncate_toward0_cached(x, 0, p10);
    std::string is = decimal_to_string(t_int, 0);
    if (!is.empty() && is[0] == '-') is.erase(0, 1);
    int d = (int)is.size();
    if (d < 1) d = 1;
    return d;
}

static Decimal arf_mid_to_decimal_exact_trunc_pow10_cached(const arf_struct *mid, int decimals, Pow10Cache& p10) {
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

    return strict_truncate_toward0_cached(r, decimals, p10);
}

static Decimal arb_mid_to_decimal_exact_trunc_cached(const arb_t x, int decimals, Pow10Cache& p10) {
    const arf_struct *mid = arb_midref(x);
    return arf_mid_to_decimal_exact_trunc_pow10_cached(mid, decimals, p10);
}

// -----------------------------
// DISPLAY helpers (final-report only)
// -----------------------------
static std::string trim_trailing_zeros_fixed(std::string s) {
    if (s.find('e') != std::string::npos || s.find('E') != std::string::npos) return s;

    size_t dot = s.find('.');
    if (dot == std::string::npos) return s;

    while (!s.empty() && s.back() == '0') s.pop_back();
    if (!s.empty() && s.back() == '.') s.pop_back();

    if (s == "-0") s = "0";
    return s;
}

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

static void print_arb_bracket(const char* label, const arb_t x, slong digits) {
    std::cout << label << ": ";
    arb_printn(x, digits, 0);
    std::cout << "\n";
}

static void fprint_arb_bracket(FILE* f, const char* label, const arb_t x, slong digits) {
    std::fprintf(f, "%s: ", label);
    arb_fprintn(f, x, digits, 0);
    std::fprintf(f, "\n");
}

static std::string arb_mid_to_fixed_decimal_adaptive_trim(const arb_t x, int sig_digits, int max_decimals) {
    const arf_struct *mid = arb_midref(x);

    if (arf_is_zero(mid)) return "0";

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

    // fallback
    char *c2 = arf_get_str(mid, sig_digits * 2);
    if (!c2) {
        std::cerr << "ERROR: arf_get_str failed\n";
        std::exit(1);
    }
    std::string sci2(c2);
    flint_free(c2);

    std::string fixed2 = sci_to_fixed_noexp_trunc(sci2, max_decimals);
    fixed2 = trim_trailing_zeros_fixed(fixed2);

    if (fixed2 == "0") return sci2;
    return fixed2;
}

// -----------------------------
// Hot-loop progress helpers (NO decimal expansion / no extra zeta calls)
// -----------------------------
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

// One short line per step (stdout).
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

// Same line, but to file (steps log).
static void fprint_step_line(
    FILE* f,
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

    std::fwrite(buf, 1, (size_t)n, f);
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

    // kept for compatibility but unused in this minimal build:
    bool profile = false;
    int profile_every = 1;
    int secant_min_frac = 6;
    int secant_stable_needed = 2;
    int secant_lock_prefix = 60;

    // Output files
    std::string out_path = "nsniperv4_output.txt";
    std::string out_steps_path = "";
    std::string out_t_path = "";
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
        } else if (k == "--profile") {
            // ignored (compat)
            a.profile = true;
        } else if (k == "--profile_every") {
            require_value(k);
            a.profile_every = std::stoi(argv[++i]);
            if (a.profile_every < 1) a.profile_every = 1;
        } else if (k == "--secant_min_frac") {
            require_value(k);
            a.secant_min_frac = std::stoi(argv[++i]);
            if (a.secant_min_frac < 0) a.secant_min_frac = 0;
        } else if (k == "--secant_stable") {
            require_value(k);
            a.secant_stable_needed = std::stoi(argv[++i]);
            if (a.secant_stable_needed < 0) a.secant_stable_needed = 0;
        } else if (k == "--secant_lock_prefix") {
            require_value(k);
            a.secant_lock_prefix = std::stoi(argv[++i]);
            if (a.secant_lock_prefix < 1) a.secant_lock_prefix = 1;
            if (a.secant_lock_prefix > 400) a.secant_lock_prefix = 400;
        } else if (k == "--out") {
            require_value(k);
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
    return a;
}

static inline std::string trim_local(std::string s) {
    while (!s.empty() && std::isspace((unsigned char)s.front())) s.erase(s.begin());
    while (!s.empty() && std::isspace((unsigned char)s.back())) s.pop_back();
    return s;
}

static std::string default_steps_path_from_out(const std::string& out) {
    return out + ".steps.txt";
}

static std::string default_t_path_from_out(const std::string& out) {
    return out + ".t.txt";
}

// -----------------------------
// Main
// -----------------------------
int main(int argc, char **argv) {
    Args args = parse_args(argc, argv);

    Pow10Cache p10;

    if (args.out_steps_path.empty()) {
        args.out_steps_path = default_steps_path_from_out(args.out_path);
    }
    if (args.out_t_path.empty()) {
        args.out_t_path = default_t_path_from_out(args.out_path);
    }

    FILE* fout = std::fopen(args.out_path.c_str(), "wb");
    if (!fout) {
        std::perror("ERROR: fopen(--out) failed");
        std::exit(1);
    }

    FILE* fsteps = std::fopen(args.out_steps_path.c_str(), "wb");
    if (!fsteps) {
        std::perror("ERROR: fopen(--out_steps) failed");
        std::fclose(fout);
        std::exit(1);
    }

    FILE* ft = std::fopen(args.out_t_path.c_str(), "wb");
    if (!ft) {
        std::perror("ERROR: fopen(--out_t) failed");
        std::fclose(fsteps);
        std::fclose(fout);
        std::exit(1);
    }

    std::cout << "OUT:       " << args.out_path << "\n";
    std::cout << "OUT_STEPS: " << args.out_steps_path << "\n";
    std::cout << "OUT_T:     " << args.out_t_path << "\n";

    std::string t_string = trim_local(args.t);
    Decimal t_work = decimal_from_string(t_string);

    int work_frac_digits = count_decimals_str(t_string);
    t_work = t_can_cached(t_work, work_frac_digits, p10);

    const int internal_max_frac = args.depth + args.overshoot_frac;

    int start_int_digits = count_int_digits_cached(t_work, p10);
    int compute_dps = (start_int_digits + internal_max_frac) + args.max_step + 256;
    slong prec_bits = dps_to_bits(compute_dps);

    ArbZetaContext ctx(prec_bits);
    const Decimal sigma = Decimal(1, 2);

    std::string dir_effective = args.dir;
    bool lock_found = false;
    bool safety_after_cap = false;

    arb_t final_slope_arb;
    arb_init(final_slope_arb);
    bool have_final_slope = false;

    arb_t last_e_pre;
    arb_init(last_e_pre);
    bool have_last_e_pre = false;

    // Header to files
    std::fprintf(fout, "nsniperv4 minimal exact-landing run log\n");
    std::fprintf(fout, "t=%s\n", args.t.c_str());
    std::fprintf(fout, "depth=%d\n", args.depth);
    std::fprintf(fout, "max_steps=%d\n", args.max_steps);
    std::fprintf(fout, "overshoot_frac=%d\n", args.overshoot_frac);
    std::fprintf(fout, "max_step=%d\n", args.max_step);
    std::fprintf(fout, "dir=%s\n", args.dir.c_str());
    std::fprintf(fout, "compute_dps=%d\n", compute_dps);
    std::fprintf(fout, "arb_prec_bits=%lld\n", (long long)prec_bits);
    std::fprintf(fout, "out_steps=%s\n", args.out_steps_path.c_str());
    std::fprintf(fout, "out_t=%s\n", args.out_t_path.c_str());
    std::fprintf(fout, "----------------------------------------\n");
    std::fflush(fout);

    std::fprintf(fsteps, "STEP TRACE (cheap): step, frac(after), step_size, probe_depth, CAP/FIT, dir, log10(|E_pre|), log10(|J|), slope sign\n");
    std::fflush(fsteps);

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
    arb_t e_ratio, e_drop_pct, one_arb, hundred_arb, tmp_arb;

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

        // ----- EXACT-LANDING STEP SCHEDULER (doubling regime, safe clamp) -----
        // Base doubling rule (unchanged):
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

        // Exact landing clamp (NEW):
        // Never overshoot args.depth. If we're within range, shrink the effective step.
        if (work_frac_digits < args.depth) {
            int remaining = args.depth - work_frac_digits;
            if (remaining < 1) remaining = 1;
            if (step_size_eff > remaining) step_size_eff = remaining;
        } else {
            // already at/over depth, force termination path below
            step_size_eff = 1;
        }

        int target_frac = work_frac_digits + step_size_eff;
        if (target_frac > args.depth) target_frac = args.depth; // hard safety

        // Probe depth logic: DO NOT CHANGE (as requested)
        int probe_depth = work_frac_digits;
        if (probe_depth < 6) probe_depth = 6;
        probe_depth += 5;
        if (probe_depth > internal_max_frac) probe_depth = internal_max_frac;
        if (probe_depth < 1) probe_depth = 1;

        // Prec policy (same pattern as your current file)
        int base_guard = 256;
        int cancel_guard = probe_depth + 64;
        int step_dps = (start_int_digits + probe_depth) + base_guard + cancel_guard;
        slong prec_step_bits = clamp_step_bits(dps_to_bits(step_dps));

        // Compute E_now (no cache)
        ctx.compute(t_work, sigma, prec_step_bits);
        arb_set(e_now, ctx.abs_z);

        arb_set(last_e_pre, e_now);
        have_last_e_pre = true;

        // Forward-diff slope only (minimal)
        mpz_class h_den = p10.get(probe_depth);
        Decimal h_dec = Decimal(1, h_den);
        ctx.arb_set_decimal_exact_prec(h_arb, h_dec, prec_step_bits);

        Decimal t_plus_can = t_can_cached(t_work + h_dec, probe_depth, p10);

        ctx.compute(t_plus_can, sigma, prec_step_bits);
        arb_set(e_plus, ctx.abs_z);

        arb_sub(diff, e_plus, e_now, prec_step_bits);
        if (arb_is_zero(diff)) {
            arb_set_str(slope, "1e-300", prec_step_bits);
        } else {
            arb_div(slope, diff, h_arb, prec_step_bits);
        }
        arb_div(raw_jump, e_now, slope, prec_step_bits);

        // CAP threshold
        Decimal cap_dec;
        if (work_frac_digits == 0) cap_dec = Decimal(1, 10);
        else cap_dec = Decimal(1, p10.get(work_frac_digits));
        ctx.arb_set_decimal_exact_prec(cap_arb, cap_dec, prec_step_bits);

        arb_abs(abs_raw_jump, raw_jump);
        bool is_cap = arb_gt(abs_raw_jump, cap_arb);

        if (is_cap && work_frac_digits > 5) safety_after_cap = true;

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
                arb_t e_cap_plus, e_cap_minus;
                arb_init(e_cap_plus);
                arb_init(e_cap_minus);

                Decimal t_cap_plus_can  = t_can_cached(t_work + cap_dec, work_frac_digits, p10);
                Decimal t_cap_minus_can = t_can_cached(t_work - cap_dec, work_frac_digits, p10);

                ctx.compute(t_cap_plus_can, sigma, prec_step_bits);
                arb_set(e_cap_plus, ctx.abs_z);
                ctx.compute(t_cap_minus_can, sigma, prec_step_bits);
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
            // exact landing: we already computed target_frac clamped to args.depth
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

        if (!have_final_slope && next_work_frac >= args.depth) {
            arb_set(final_slope_arb, slope);
            have_final_slope = true;
        }

        Decimal t_next_dec = t_work;

        if (is_cap) {
            if (arb_is_negative(safe_jump)) t_next_dec = t_work + cap_dec;
            else t_next_dec = t_work - cap_dec;
        } else {
            int jump_trunc_decimals = next_work_frac + 256;
            Decimal jump_dec = arb_mid_to_decimal_exact_trunc_cached(safe_jump, jump_trunc_decimals, p10);
            t_next_dec = t_work - jump_dec;
        }

        // Maintain original internal_max_frac safety, but never store beyond args.depth
        if (next_work_frac > internal_max_frac) next_work_frac = internal_max_frac;
        if (next_work_frac > args.depth) next_work_frac = args.depth;

        t_work = strict_truncate_toward0_cached(t_next_dec, next_work_frac, p10);

        if (!args.quiet) {
            print_step_line(total_steps, next_work_frac, step_size_eff, probe_depth, is_cap, dir_effective, e_now, safe_jump, slope);
        }
        fprint_step_line(fsteps, total_steps, next_work_frac, step_size_eff, probe_depth, is_cap, dir_effective, e_now, safe_jump, slope);
        if ((total_steps % 50) == 0) std::fflush(fsteps);

        work_frac_digits = next_work_frac;

        if (work_frac_digits >= args.depth) {
            Decimal t_final = t_work;
            if (work_frac_digits > args.depth) {
                t_final = strict_truncate_toward0_cached(t_work, args.depth, p10);
            }

            ctx.compute(t_final, sigma, prec_bits);

            if (!have_last_e_pre || arb_is_zero(last_e_pre)) {
                arb_zero(e_ratio);
                arb_zero(e_drop_pct);
            } else {
                arb_div(e_ratio, ctx.abs_z, last_e_pre, prec_bits);
                arb_sub(tmp_arb, one_arb, e_ratio, prec_bits);
                arb_mul(e_drop_pct, tmp_arb, hundred_arb, prec_bits);
            }

            std::cout << ">>> REACHED TARGET DEPTH (SETTLED) <<<\n";
            print_arb_bracket("FINAL SINK ENERGY (SCI, MID+/-RAD)", ctx.abs_z, 50);
            std::cout << "WROTE FULL OUTPUT TO: " << args.out_path << "\n";
            std::cout << "WROTE STEP TRACE TO:  " << args.out_steps_path << "\n";
            std::cout << "WROTE T-ONLY TO:      " << args.out_t_path << "\n";

            std::fprintf(fout, ">>> REACHED TARGET DEPTH (SETTLED) <<<\n");
            std::fprintf(fout, "----------------------------------------\n");
            fprint_arb_bracket(fout, "FINAL SINK ENERGY (SCI, MID+/-RAD)", ctx.abs_z, 80);

            std::string final_energy_fixed = arb_mid_to_fixed_decimal_adaptive_trim(ctx.abs_z, 120, args.depth + 256);
            std::fprintf(fout, "FINAL SINK ENERGY (MID, FIXED): %s\n", final_energy_fixed.c_str());

            std::string t_final_str = fixed_trunc_str_decimal_exact_cached(t_final, args.depth, p10);
            std::fprintf(fout, "FINAL T VALUE:     %s\n", t_final_str.c_str());

            // T-only file: exactly digits + '\n'
            std::fwrite(t_final_str.data(), 1, t_final_str.size(), ft);
            std::fflush(ft);

            if (!have_final_slope) {
                arb_set(final_slope_arb, slope);
                have_final_slope = true;
            }
            Decimal final_slope_mid_dec = arb_mid_to_decimal_exact_trunc_cached(final_slope_arb, args.depth, p10);
            std::string slope_str = fixed_trunc_str_decimal_exact_cached(final_slope_mid_dec, args.depth, p10);
            std::fprintf(fout, "FINAL_SLOPE(mid):  %s\n", slope_str.c_str());

            Decimal final_energy_ratio_dec = arb_mid_to_decimal_exact_trunc_cached(e_ratio, args.depth, p10);
            Decimal final_energy_drop_pct_dec = arb_mid_to_decimal_exact_trunc_cached(e_drop_pct, args.depth, p10);

            std::string ratio_out = fixed_trunc_str_decimal_exact_cached(final_energy_ratio_dec, args.depth, p10);
            std::string drop_out  = fixed_trunc_str_decimal_exact_cached(final_energy_drop_pct_dec, args.depth, p10);
            ratio_out = trim_trailing_zeros_fixed(ratio_out);
            drop_out  = trim_trailing_zeros_fixed(drop_out);

            std::fprintf(fout, "FINAL_ENERGY_RATIO:    %s\n", ratio_out.c_str());
            std::fprintf(fout, "FINAL_ENERGY_DROP_PCT: %s\n", drop_out.c_str());

            std::fflush(fout);
            std::fflush(fsteps);

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

            std::fclose(ft);
            std::fclose(fsteps);
            std::fclose(fout);
            return 0;
        }
    }

    // max_steps reached without depth: do one final compute and log to file
    ctx.compute(t_work, sigma, prec_bits);

    std::cout << ">>> MAX_STEPS REACHED (FINAL COMPUTE) <<<\n";
    print_arb_bracket("FINAL SINK ENERGY (SCI, MID+/-RAD)", ctx.abs_z, 50);
    std::cout << "WROTE FULL OUTPUT TO: " << args.out_path << "\n";
    std::cout << "WROTE STEP TRACE TO:  " << args.out_steps_path << "\n";
    std::cout << "WROTE T-ONLY TO:      " << args.out_t_path << "\n";

    std::fprintf(fout, ">>> MAX_STEPS REACHED (FINAL COMPUTE) <<<\n");
    std::fprintf(fout, "----------------------------------------\n");
    fprint_arb_bracket(fout, "FINAL SINK ENERGY (SCI, MID+/-RAD)", ctx.abs_z, 80);

    std::string raw_energy_fixed = arb_mid_to_fixed_decimal_adaptive_trim(ctx.abs_z, 120, args.depth);
    std::fprintf(fout, "FINAL SINK ENERGY (RAW,FIXED): %s\n", raw_energy_fixed.c_str());

    std::string t_str = fixed_trunc_str_decimal_exact_cached(t_work, work_frac_digits, p10);
    std::fprintf(fout, "FINAL T VALUE:                 %s\n", t_str.c_str());

    // T-only file
    std::fwrite(t_str.data(), 1, t_str.size(), ft);
    std::fflush(ft);

    std::fflush(fout);
    std::fflush(fsteps);

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

    std::fclose(ft);
    std::fclose(fsteps);
    std::fclose(fout);
    return 0;
}
