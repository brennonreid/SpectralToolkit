// heatsonar.cpp
//
// Logical heatmap + seed extraction + jet-spider refinement for Hardy Z zeros,
// using SINGLE-CALL (jet) evaluations.
//
// Each evaluation uses exactly one call to:
//   ctx.compute_hardy_z_and_deriv_prec(T, prec_bits);
//
// Primary use-cases:
//  1) Range scan over integer interval [int_start, int_end] with step (e.g., 0.01 or 0.001)
//     -> emit "hits" (local minima) and write candidates to a seeds file.
//  2) Local spider around a center T0 using heatmap-derived seeds.
//
// STRICT POLICY:
// - If your ArbZetaContext does not provide the required methods, compilation MUST fail.
// - No fallbacks.
//
// REQUIRED ArbZetaContext API (from your arb_interface.hpp):
//   void compute_hardy_z_and_deriv_prec(const Decimal& T, slong prec_bits);
//   arb_t hardy_z;   // Z(T) stored as arb_t
//   arb_t hardy_zp;  // Z'(T) stored as arb_t
//   std::string get_energy_str(int digits) const; // returns |Z| string from last compute
//
// Build example (MSYS2 UCRT64):
//   g++ -O2 -std=c++17 heatsonar.cpp C:/CToolkitV1/zeta_mpc.cpp -o heatsonar.exe \
//     -IC:/msys64/ucrt64/include -LC:/msys64/ucrt64/lib -lmpc -lmpfr -lgmp -lgmpxx -lflint
//
// Example scan 10..100 step 0.001, keep ALL minima hits:
//   ./heatsonar.exe --mode range --int_start 10 --int_end 100 --T_step 0.001 --seeds seeds.txt --csv hits.csv
//
// Example scan 10..100 step 0.001, keep minima with |Z| <= 1e-6:
//   ./heatsonar.exe --mode range --int_start 10 --int_end 100 --T_step 0.001 --target_energy 1e-6 --seeds seeds.txt --csv hits.csv
//
// DEFAULT OUTPUTS (minimal-change extension):
// - By default, the program writes TWO files unless you override them:
//     1) heatmap file:     --heatmap <path>   (default: heatmap.csv)
//     2) seeds file:       --seeds   <path>   (default: seeds.txt)
// - Existing outputs remain intact:
//     --csv still writes "hits" CSV in range mode,
//     and writes spider report CSV in spider mode.
//
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <limits>
#include <cstdint>

#include "decimal.hpp"
#include "arb_interface.hpp"

#include <flint/flint.h>
#include <flint/arf.h>
#include <flint/arb.h>

// ---------------------------
// Small helpers
// ---------------------------
static inline long double ld_abs(long double x) { return (x < 0) ? -x : x; }

static inline std::string ld_to_sci(long double x, int prec = 18) {
    std::ostringstream oss;
    oss.setf(std::ios::scientific);
    oss << std::setprecision(prec) << (double)x;
    return oss.str();
}

// Convert arb_t midpoint to long double (via double)
static inline long double arb_mid_to_ld(const arb_t x) {
    const arf_struct *mid = arb_midref(x);
    return (long double)arf_get_d(mid, ARF_RND_NEAR);
}

// Render arb_t to string (no exponent), flint allocates
static inline std::string arb_to_str_noexp(const arb_t x, slong digits) {
    char *ss = arb_get_str(x, digits, ARF_STR_NOEXPONENT);
    std::string res(ss);
    flint_free(ss);
    return res;
}

// Count decimals in a Decimal by string inspection.
static int count_decimals(const Decimal &T) {
    std::string s = decimal_to_string(T, 80);
    return count_decimals_str(s);
}

static Decimal truncate_T(const Decimal &T, int digits) {
    return truncate_decimal(T, digits);
}

// ---------------------------
// Jet evaluation (single call)
// ---------------------------
struct JetEval {
    Decimal T;
    long double Z;       // Hardy Z(T)
    long double Zp;      // d/dT Hardy Z
    long double E;       // |Z|
    long double slope;   // d|Z|/dT = sign(Z)*Zp
    int slope_sign;      // sign(slope)
    std::string Z_str;   // pretty Z(T)
    std::string E_str;   // pretty |Z(T)|
};

static JetEval eval_jet(ArbZetaContext &ctx, const Decimal &T, slong prec_bits, int pretty_digits) {
    ctx.compute_hardy_z_and_deriv_prec(T, prec_bits);

    JetEval j;
    j.T = T;

    // ArbZetaContext stores hardy_z and hardy_zp as arb_t
    j.Z  = arb_mid_to_ld(ctx.hardy_z);
    j.Zp = arb_mid_to_ld(ctx.hardy_zp);

    j.E = ld_abs(j.Z);

    long double signZ = (j.Z >= 0) ? 1.0L : -1.0L;
    j.slope = signZ * j.Zp;
    j.slope_sign = (j.slope > 0) ? 1 : (j.slope < 0 ? -1 : 0);

    // Pretty strings for reporting
    j.Z_str = arb_to_str_noexp(ctx.hardy_z, pretty_digits);
    j.E_str = ctx.get_energy_str(pretty_digits);

    return j;
}

// ---------------------------
// Logical heatmap cell
// ---------------------------
struct HeatCell {
    Decimal T;
    long double log10E;
    long double log10Jump;   // log10(|E/slope|) or +inf if slope==0
    int dir;                 // +1 means increasing T reduces |Z|, -1 means decreasing T reduces |Z|, 0 unknown
    bool cap;
    long double E;
    long double slope;
};

static HeatCell make_cell(const JetEval &j, long double max_step_ld) {
    HeatCell c;
    c.T = j.T;
    c.E = j.E;
    c.slope = j.slope;

    if (j.E <= 0) c.log10E = -std::numeric_limits<long double>::infinity();
    else c.log10E = std::log10((double)j.E);

    if (j.slope == 0.0L) {
        c.log10Jump = std::numeric_limits<long double>::infinity();
        c.dir = 0;
        c.cap = true;
        return c;
    }

    long double raw_jump = j.E / j.slope; // Newton step in T units

    // Update is T_next = T - raw_jump.
    int raw_sign = (raw_jump > 0) ? 1 : (raw_jump < 0 ? -1 : 0);
    c.dir = (raw_sign == 0) ? 0 : -raw_sign;

    long double aj = ld_abs(raw_jump);
    c.log10Jump = (aj <= 0) ? -std::numeric_limits<long double>::infinity() : std::log10((double)aj);

    // Logical cap: compare jump magnitude vs max_step
    c.cap = (aj > max_step_ld);

    return c;
}

// Local minima indices in log10E (3-point).
static std::vector<size_t> find_local_minima(const std::vector<HeatCell> &cells) {
    std::vector<size_t> idx;
    if (cells.size() < 3) return idx;
    for (size_t i = 1; i + 1 < cells.size(); ++i) {
        if (cells[i].log10E <= cells[i-1].log10E && cells[i].log10E <= cells[i+1].log10E) {
            idx.push_back(i);
        }
    }
    return idx;
}

// Write hits CSV
static void write_cells_csv(const std::string &path, const std::vector<HeatCell> &cells) {
    std::ofstream f(path);
    if (!f) {
        std::cerr << "ERROR: cannot open csv: " << path << std::endl;
        std::exit(1);
    }
    f << "T,log10E,log10Jump,dir,cap,E,slope\n";
    for (const auto &c : cells) {
        f << decimal_to_string(c.T, 40) << ","
          << ld_to_sci(c.log10E, 18) << ","
          << ld_to_sci(c.log10Jump, 18) << ","
          << c.dir << ","
          << (c.cap ? 1 : 0) << ","
          << ld_to_sci(c.E, 18) << ","
          << ld_to_sci(c.slope, 18) << "\n";
    }
}

// Write FULL heatmap CSV (every sampled cell).
static void write_heatmap_csv(const std::string &path, const std::vector<HeatCell> &cells) {
    std::ofstream f(path);
    if (!f) {
        std::cerr << "ERROR: cannot open heatmap: " << path << std::endl;
        std::exit(1);
    }
    f << "idx,T,log10E,log10Jump,dir,cap,E,slope\n";
    for (size_t i = 0; i < cells.size(); ++i) {
        const auto &c = cells[i];
        f << i << ","
          << decimal_to_string(c.T, 40) << ","
          << ld_to_sci(c.log10E, 18) << ","
          << ld_to_sci(c.log10Jump, 18) << ","
          << c.dir << ","
          << (c.cap ? 1 : 0) << ","
          << ld_to_sci(c.E, 18) << ","
          << ld_to_sci(c.slope, 18) << "\n";
    }
}

// Write seed list (one T per line).
static void write_seeds_txt(const std::string &path, const std::vector<Decimal> &seeds, int digits) {
    std::ofstream f(path);
    if (!f) {
        std::cerr << "ERROR: cannot open seeds file: " << path << std::endl;
        std::exit(1);
    }
    for (const auto &t : seeds) {
        f << decimal_to_string(t, digits) << "\n";
    }
}

// ---------------------------
// Jet spider refinement (single-call per iteration)
// ---------------------------
struct SpiderArgs {
    Decimal T0;
    Decimal band_radius;   // +- radius for heatmap band
    Decimal T_step;        // heatmap sampling step
    int heat_points_max;   // guard
    int seed_count;        // number of minima to spider
    int depth;             // truncate T to this many fractional digits each step
    int max_steps;         // per-seed steps
    Decimal max_step;      // cap magnitude in T units
    long double target_energy; // stop if E <= target_energy, negative disables
    int dps;               // precision driver (decimal digits -> bits)
    int pretty_digits;     // digits for string reporting
};

struct SpiderResult {
    Decimal T;
    long double E;
    std::string E_str;
    int steps;
};

static SpiderResult spider_one_seed(ArbZetaContext &ctx, const Decimal &seedT, const SpiderArgs &a) {
    Decimal T = seedT;

    SpiderResult out;
    out.T = T;
    out.E = std::numeric_limits<long double>::infinity();
    out.steps = 0;

    long double max_step_ld = (long double)std::stod(decimal_to_string(a.max_step, 40));

    for (int k = 0; k < a.max_steps; ++k) {
        slong prec_bits = (slong)(a.dps * 3.33 + 20);

        JetEval j = eval_jet(ctx, T, prec_bits, a.pretty_digits);

        out.E = j.E;
        out.E_str = j.E_str;
        out.T = T;
        out.steps = k;

        if (a.target_energy > 0.0L && j.E <= a.target_energy) break;
        if (j.slope == 0.0L) break;

        long double raw_jump_ld = (j.E / j.slope);
        Decimal raw_jump = Decimal((double)raw_jump_ld);

        long double aj_ld = ld_abs(raw_jump_ld);
        if (aj_ld > max_step_ld) {
            raw_jump = (raw_jump < Decimal(0)) ? -a.max_step : a.max_step;
        }

        T = T - raw_jump;

        if (a.depth >= 0) {
            T = truncate_T(T, a.depth);
        }
    }

    return out;
}

// ---------------------------
// MODE: range
// ---------------------------
static void run_range_mode(ArbZetaContext &ctx,
                          int int_start,
                          int int_end,
                          const Decimal &T_step_in,
                          int dps,
                          int pretty_digits,
                          const Decimal &max_step,
                          long double target_energy,
                          const std::string &heatmap_path,
                          const std::string &csv_path,
                          const std::string &seeds_path,
                          const std::string &order) {
    if (int_end < int_start) {
        std::cerr << "ERROR: int_end must be >= int_start\n";
        std::exit(1);
    }
    if (T_step_in == Decimal(0)) {
        std::cerr << "ERROR: T_step must be non-zero\n";
        std::exit(1);
    }

    Decimal T_step = decimal_abs(T_step_in);
    long double max_step_ld = (long double)std::stod(decimal_to_string(max_step, 40));

    std::vector<HeatCell> cells;
    cells.reserve((size_t)((int_end - int_start + 1) * 1000));

    Decimal T = Decimal(int_start);
    Decimal T_end = Decimal(int_end);

    std::cout << "[RANGE] int_start=" << int_start
              << " int_end=" << int_end
              << " step=" << decimal_to_string(T_step, 20)
              << " dps=" << dps
              << std::endl;

    std::int64_t safety = 0;
    const std::int64_t safety_limit = 400000000; // hard guard

    while (T <= T_end) {
        slong prec_bits = (slong)(dps * 3.33 + 20);
        JetEval j = eval_jet(ctx, T, prec_bits, pretty_digits);
        cells.push_back(make_cell(j, max_step_ld));

        T = T + T_step;

        if (++safety > safety_limit) {
            std::cerr << "ERROR: range sweep safety limit hit\n";
            std::exit(1);
        }
    }

    std::cout << "[RANGE] samples=" << cells.size() << std::endl;

    // Always write full heatmap file (defaulted by main unless overridden).
    if (!heatmap_path.empty()) {
        write_heatmap_csv(heatmap_path, cells);
        std::cout << "[RANGE] wrote heatmap: " << heatmap_path << std::endl;
    }

    auto mins = find_local_minima(cells);
    std::cout << "[RANGE] local_minima=" << mins.size() << std::endl;

    std::vector<HeatCell> hits;
    hits.reserve(mins.size());
    for (auto i : mins) hits.push_back(cells[i]);

    // Filter by threshold (if provided). Otherwise keep ALL minima hits.
    if (target_energy > 0.0L) {
        std::vector<HeatCell> filtered;
        filtered.reserve(hits.size());
        for (const auto &h : hits) {
            if (h.E <= target_energy) filtered.push_back(h);
        }
        hits.swap(filtered);
        std::cout << "[RANGE] hits_after_target_energy=" << hits.size()
                  << " (target_energy=" << ld_to_sci(target_energy, 6) << ")"
                  << std::endl;
    }

    // Output ordering
    if (order == "energy") {
        std::sort(hits.begin(), hits.end(), [](const HeatCell &a, const HeatCell &b){
            return a.E < b.E;
        });
    } else if (order == "T" || order == "t") {
        std::sort(hits.begin(), hits.end(), [](const HeatCell &a, const HeatCell &b){
            return a.T < b.T;
        });
    } else {
        std::cerr << "ERROR: --order must be T or energy\n";
        std::exit(1);
    }

    int show = (int)std::min<size_t>(20, hits.size());
    for (int k = 0; k < show; ++k) {
        const auto &c = hits[(size_t)k];
        std::cout << "  [HIT " << k << "] T=" << decimal_to_string(c.T, 25)
                  << " log10E=" << ld_to_sci(c.log10E, 18)
                  << " dir=" << c.dir
                  << " log10Jump=" << ld_to_sci(c.log10Jump, 18)
                  << (c.cap ? " [CAP]" : "")
                  << std::endl;
    }

    // Existing: --csv writes HITS (minima) CSV in range mode.
    if (!csv_path.empty()) {
        write_cells_csv(csv_path, hits);
        std::cout << "[RANGE] wrote hits csv: " << csv_path << std::endl;
    }

    // Existing: --seeds writes candidate T list; now defaulted to seeds.txt unless overridden.
    if (!seeds_path.empty()) {
        std::vector<Decimal> seeds;
        seeds.reserve(hits.size());
        for (const auto &h : hits) seeds.push_back(h.T);

        write_seeds_txt(seeds_path, seeds, 80);
        std::cout << "[RANGE] wrote seeds: " << seeds_path << " (count=" << seeds.size() << ")\n";
    }

    if (hits.empty()) {
        std::cout << "[RANGE] WARNING: no hits matched your filter.\n";
        if (target_energy > 0.0L) {
            std::cout << "         If you used --target_energy, try a larger threshold.\n";
        }
    }
}

// ---------------------------
// MODE: spider (local band)
// ---------------------------
static void run_spider_mode(ArbZetaContext &ctx,
                            const SpiderArgs &a,
                            const std::string &heatmap_path,
                            const std::string &csv_path,
                            const std::string &seeds_path) {
    Decimal left = a.T0 - a.band_radius;
    Decimal right = a.T0 + a.band_radius;

    long double max_step_ld = (long double)std::stod(decimal_to_string(a.max_step, 40));

    std::vector<HeatCell> cells;
    cells.reserve((size_t)a.heat_points_max);

    Decimal T = left;
    Decimal step = a.T_step;
    if (step == Decimal(0)) {
        std::cerr << "ERROR: T_step must be non-zero" << std::endl;
        std::exit(1);
    }

    while (T <= right) {
        if ((int)cells.size() >= a.heat_points_max) break;

        slong prec_bits = (slong)(a.dps * 3.33 + 20);
        JetEval j = eval_jet(ctx, T, prec_bits, a.pretty_digits);
        cells.push_back(make_cell(j, max_step_ld));

        T = T + step;
    }

    std::cout << "[SPIDER] band=[" << decimal_to_string(left, 20) << ", " << decimal_to_string(right, 20) << "]"
              << " step=" << decimal_to_string(step, 20)
              << " samples=" << cells.size() << std::endl;

    // Always write full heatmap file (defaulted by main unless overridden).
    if (!heatmap_path.empty()) {
        write_heatmap_csv(heatmap_path, cells);
        std::cout << "[SPIDER] wrote heatmap: " << heatmap_path << std::endl;
    }

    auto mins = find_local_minima(cells);
    if (mins.empty()) {
        std::cerr << "[SPIDER] no local minima found in heatmap band" << std::endl;
        std::exit(2);
    }

    struct Seed { Decimal T; long double log10E; };
    std::vector<Seed> seeds;
    seeds.reserve(mins.size());
    for (auto idx : mins) seeds.push_back({cells[idx].T, cells[idx].log10E});
    std::sort(seeds.begin(), seeds.end(), [&](const Seed &aa, const Seed &bb){
        return aa.log10E < bb.log10E;
    });

    int nseed = std::min(a.seed_count, (int)seeds.size());
    std::cout << "[SPIDER] seeds=" << nseed << " (from " << seeds.size() << " minima)" << std::endl;

    std::vector<SpiderResult> out;
    out.reserve((size_t)nseed);

    for (int i = 0; i < nseed; ++i) {
        Decimal seedT = seeds[i].T;
        std::cout << "  [SEED " << i << "] T=" << decimal_to_string(seedT, 40)
                  << " log10E=" << ld_to_sci(seeds[i].log10E, 18) << std::endl;

        SpiderResult r = spider_one_seed(ctx, seedT, a);
        out.push_back(r);

        std::cout << "    [DONE] steps=" << r.steps
                  << " T=" << decimal_to_string(r.T, 60)
                  << " E=" << r.E_str
                  << std::endl;
    }

    std::sort(out.begin(), out.end(), [&](const SpiderResult &aa, const SpiderResult &bb){
        return aa.E < bb.E;
    });

    std::cout << "\n[SPIDER REPORT]\n";
    for (size_t i = 0; i < out.size(); ++i) {
        std::cout << "  #" << i
                  << " steps=" << out[i].steps
                  << " T=" << decimal_to_string(out[i].T, 60)
                  << " E=" << out[i].E_str
                  << "\n";
    }

    // Existing: --csv writes spider REPORT csv (ranked refined outputs).
    if (!csv_path.empty()) {
        std::ofstream f(csv_path);
        if (!f) { std::cerr << "ERROR: cannot open csv: " << csv_path << "\n"; std::exit(1); }
        f << "rank,steps,T,E_str\n";
        for (size_t i = 0; i < out.size(); ++i) {
            f << i << "," << out[i].steps << "," << decimal_to_string(out[i].T, 80) << "," << out[i].E_str << "\n";
        }
        std::cout << "[SPIDER] wrote csv: " << csv_path << std::endl;
    }

    // Existing: --seeds writes final refined candidates; now defaulted to seeds.txt unless overridden.
    if (!seeds_path.empty()) {
        std::vector<Decimal> s;
        s.reserve(out.size());
        for (const auto &r : out) s.push_back(r.T);
        write_seeds_txt(seeds_path, s, 80);
        std::cout << "[SPIDER] wrote seeds: " << seeds_path << " (count=" << s.size() << ")\n";
    }
}

// ---------------------------
// CLI
// ---------------------------
static void print_help() {
    std::cout <<
    "Usage: heatsonar --mode range|spider [options]\n"
    "\n"
    "Modes:\n"
    "  range : scan [int_start..int_end] with --T_step, extract local minima hits, write seeds\n"
    "  spider: local band heatmap around --T, pick seeds, refine by jet-spider\n"
    "\n"
    "Default outputs (2 files unless overridden):\n"
    "  --heatmap <path>         Full heatmap CSV of every sampled cell (default: heatmap.csv)\n"
    "  --seeds <path>           Seeds/candidates (one T per line) (default: seeds.txt)\n"
    "\n"
    "Common options:\n"
    "  --dps <int>              Decimal digits for Arb precision (default 80)\n"
    "  --pretty_digits <int>    Digits for string reporting (default 25)\n"
    "  --max_step <num>         Cap magnitude for spider jump in T units (default 200)\n"
    "  --csv <path>             Existing CSV output (range: hits CSV, spider: report CSV)\n"
    "  --seeds <path>           Write seeds file (one T per line)\n"
    "  --heatmap <path>         Write full heatmap CSV (every sampled cell)\n"
    "\n"
    "Range options:\n"
    "  --int_start <int>\n"
    "  --int_end <int>\n"
    "  --T_step <num>           Step size (e.g., 0.01 or 0.001)\n"
    "  --target_energy <dbl>    Keep hits with |Z| <= target_energy (optional)\n"
    "  --order <T|energy>       Ordering of displayed/written hits (default T)\n"
    "\n"
    "Spider options:\n"
    "  --T <num>                Center\n"
    "  --band_radius <num>      Band radius (default 0.2)\n"
    "  --T_step <num>           Heatmap sampling step (default 1e-3)\n"
    "  --heat_points_max <int>  Max samples in band (default 200000)\n"
    "  --seed_count <int>       Number of minima to spider (default 6)\n"
    "  --depth <int>            Truncate T to this many fractional digits (default 120)\n"
    "  --max_steps <int>        Max spider steps per seed (default 500)\n"
    "  --target_energy <dbl>    Stop spider if |Z| <= target_energy (optional)\n";
}

int main(int argc, char** argv) {
    if (argc < 2) { print_help(); return 1; }

    std::string mode = "";
    std::string csv_path = "";
    std::string seeds_path = "";      // seeds/candidates file
    std::string heatmap_path = "";    // full heatmap CSV
    std::string order = "T";          // T (forward) or energy

    // Defaults
    int dps = 80;
    int pretty_digits = 25;
    Decimal max_step = decimal_from_string("200");

    // Range defaults
    int int_start = 0;
    int int_end = -1;
    std::string T_step_str = "";
    long double target_energy = -1.0L;

    // Spider defaults
    std::string T_str = "";
    Decimal band_radius = decimal_from_string("0.2");
    Decimal spider_T_step = decimal_from_string("1e-3");
    int heat_points_max = 200000;
    int seed_count = 6;
    int depth = 120;
    int max_steps = 500;

    // Manual arg parse (strict)
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        auto need = [&](const std::string &name) {
            if (i + 1 >= argc) { std::cerr << "ERROR: missing value for " << name << std::endl; std::exit(1); }
        };

        if (arg == "--mode") { need(arg); mode = argv[++i]; }
        else if (arg == "--csv") { need(arg); csv_path = argv[++i]; }
        else if (arg == "--seeds") { need(arg); seeds_path = argv[++i]; }
        else if (arg == "--heatmap") { need(arg); heatmap_path = argv[++i]; }
        else if (arg == "--dps") { need(arg); dps = std::stoi(argv[++i]); }
        else if (arg == "--pretty_digits") { need(arg); pretty_digits = std::stoi(argv[++i]); }
        else if (arg == "--max_step") { need(arg); max_step = decimal_from_string(argv[++i]); }

        // Range
        else if (arg == "--int_start") { need(arg); int_start = std::stoi(argv[++i]); }
        else if (arg == "--int_end") { need(arg); int_end = std::stoi(argv[++i]); }
        else if (arg == "--T_step") { need(arg); T_step_str = argv[++i]; }
        else if (arg == "--target_energy") { need(arg); target_energy = (long double)std::stod(argv[++i]); }
        else if (arg == "--order") { need(arg); order = argv[++i]; }

        // Spider
        else if (arg == "--T") { need(arg); T_str = argv[++i]; }
        else if (arg == "--band_radius") { need(arg); band_radius = decimal_from_string(argv[++i]); }
        else if (arg == "--heat_points_max") { need(arg); heat_points_max = std::stoi(argv[++i]); }
        else if (arg == "--seed_count") { need(arg); seed_count = std::stoi(argv[++i]); }
        else if (arg == "--depth") { need(arg); depth = std::stoi(argv[++i]); }
        else if (arg == "--max_steps") { need(arg); max_steps = std::stoi(argv[++i]); }

        else {
            std::cerr << "ERROR: unknown arg: " << arg << std::endl;
            std::exit(1);
        }
    }

    if (mode.empty()) { std::cerr << "ERROR: --mode required\n"; return 1; }

    // Defaults: always write 2 files unless overridden.
    if (heatmap_path.empty()) heatmap_path = "heatmap.csv";
    if (seeds_path.empty()) seeds_path = "seeds.txt";

    slong prec_bits = (slong)(dps * 3.33 + 20);
    ArbZetaContext ctx(prec_bits);

    if (mode == "range") {
        if (int_end < int_start) { std::cerr << "ERROR: range requires --int_start and --int_end\n"; return 1; }
        if (T_step_str.empty()) { std::cerr << "ERROR: range requires --T_step\n"; return 1; }
        Decimal T_step = decimal_from_string(T_step_str);

        run_range_mode(ctx, int_start, int_end, T_step, dps, pretty_digits, max_step,
                       target_energy, heatmap_path, csv_path, seeds_path, order);
        return 0;
    }

    if (mode == "spider") {
        if (T_str.empty()) { std::cerr << "ERROR: spider requires --T\n"; return 1; }
        Decimal T0 = decimal_from_string(T_str);

        if (!T_step_str.empty()) spider_T_step = decimal_from_string(T_step_str);

        SpiderArgs a;
        a.T0 = T0;
        a.band_radius = band_radius;
        a.T_step = spider_T_step;
        a.heat_points_max = heat_points_max;
        a.seed_count = seed_count;
        a.depth = depth;
        a.max_steps = max_steps;
        a.max_step = max_step;
        a.target_energy = target_energy;
        a.dps = dps;
        a.pretty_digits = pretty_digits;

        run_spider_mode(ctx, a, heatmap_path, csv_path, seeds_path);
        return 0;
    }

    std::cerr << "ERROR: unknown mode: " << mode << "\n";
    return 1;
}
