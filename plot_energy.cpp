// plot_energy.cpp
//
// C++ energy landscape port using FLINT/Arb (acb_zeta).
//
// Computes:
//   E(sigma) = sum_{k} |zeta(sigma + i*T_k)|
// over a band of sigma centered at 0.5.
//
// Output: CSV with columns "sigma,energy".
//
// Build example (UCRT64 shell):
//   g++ -std=c++17 -g \
//       /c/CToolkitV1/plot_energy.cpp \
//       -o /c/CToolkitV1/plot_energy.exe \
//       -I/ucrt64/include \
//       -L/ucrt64/lib \
//       -lflint -lmpfr -lgmp
//
// Run example:
//   ./plot_energy.exe zeros.txt --zoom 0.01 --resolution 201 \
//       --subset 50 --start_idx 50000 --precision 60 \
//       --output energy_zeros50000.csv

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <limits>
#include <algorithm>
#include <cctype>

#include <flint/flint.h>
#include <flint/arb.h>
#include <flint/arf.h>
#include <flint/acb.h>

// ----------------------- Config / CLI parsing -----------------------

struct Config {
    std::string zeros_file;
    int    resolution   = 201;        // --resolution
    int    subset       = 50;         // --subset
    int    precision    = 60;         // --precision (decimal digits)
    double zoom         = 0.01;       // --zoom (sigma total width)
    long long start_idx = 50000;      // --start_idx
    std::string output  = "energy_landscape.csv"; // --output
};

static void print_usage(const char* prog) {
    std::cerr << "Usage: " << prog << " zeros_file [options]\n\n"
              << "Options:\n"
              << "  --resolution N   Number of sigma samples (default 201)\n"
              << "  --subset N       Number of zeros to use (default 50)\n"
              << "  --precision N    Decimal digits for Arb (default 60)\n"
              << "  --zoom Z         Total sigma width around 0.5 (default 0.01)\n"
              << "  --start_idx K    Starting zero index (default 50000)\n"
              << "  --output FILE    Output CSV filename (default energy_landscape.csv)\n";
}

static bool parse_int(const char* s, int& out) {
    char* end = nullptr;
    long v = std::strtol(s, &end, 10);
    if (!end || *end != '\0') return false;
    out = static_cast<int>(v);
    return true;
}

static bool parse_long_long(const char* s, long long& out) {
    char* end = nullptr;
    long long v = std::strtoll(s, &end, 10);
    if (!end || *end != '\0') return false;
    out = v;
    return true;
}

static bool parse_double(const char* s, double& out) {
    char* end = nullptr;
    double v = std::strtod(s, &end);
    if (!end || *end != '\0') return false;
    out = v;
    return true;
}

static bool parse_args(int argc, char** argv, Config& cfg) {
    if (argc < 2) {
        print_usage(argv[0]);
        return false;
    }

    cfg.zeros_file = argv[1];

    int i = 2;
    while (i < argc) {
        std::string arg = argv[i];

        auto need_val = [&](const char* name) -> const char* {
            if (i + 1 >= argc) {
                std::cerr << "Missing value after " << name << "\n";
                std::exit(1);
            }
            return argv[++i];
        };

        if (arg == "--resolution") {
            const char* v = need_val("--resolution");
            if (!parse_int(v, cfg.resolution)) {
                std::cerr << "Invalid --resolution: " << v << "\n";
                return false;
            }
        } else if (arg == "--subset") {
            const char* v = need_val("--subset");
            if (!parse_int(v, cfg.subset)) {
                std::cerr << "Invalid --subset: " << v << "\n";
                return false;
            }
        } else if (arg == "--precision") {
            const char* v = need_val("--precision");
            if (!parse_int(v, cfg.precision)) {
                std::cerr << "Invalid --precision: " << v << "\n";
                return false;
            }
        } else if (arg == "--zoom") {
            const char* v = need_val("--zoom");
            if (!parse_double(v, cfg.zoom)) {
                std::cerr << "Invalid --zoom: " << v << "\n";
                return false;
            }
        } else if (arg == "--start_idx") {
            const char* v = need_val("--start_idx");
            if (!parse_long_long(v, cfg.start_idx)) {
                std::cerr << "Invalid --start_idx: " << v << "\n";
                return false;
            }
        } else if (arg == "--output") {
            const char* v = need_val("--output");
            cfg.output = v;
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            return false;
        }

        ++i;
    }

    return true;
}

// ----------------------- Zero loading -----------------------
//
// We store zeros as strings so we can feed them directly into Arb
// with arb_set_str() at full precision.

static std::vector<std::string> load_zeros_strings(const std::string& path) {
    std::vector<std::string> zeros;
    std::ifstream in(path);
    if (!in) {
        std::cerr << "[FATAL] Could not open zeros file: " << path << "\n";
        return zeros;
    }

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;

        std::string token;
        std::string temp;
        for (char c : line) {
            if (std::isspace(static_cast<unsigned char>(c))) {
                if (!temp.empty()) {
                    token = temp;
                    temp.clear();
                }
            } else {
                temp.push_back(c);
            }
        }
        if (!temp.empty()) {
            token = temp;
        }

        if (token.empty()) continue;

        // Quick sanity check that token is numeric:
        char* end = nullptr;
        (void)std::strtod(token.c_str(), &end);
        if (end && *end == '\0') {
            zeros.push_back(token);
        }
    }

    // Assume file is already sorted; no numeric sort here to avoid double.
    return zeros;
}

// ----------------------- Main -----------------------

int main(int argc, char** argv) {
    Config cfg;
    if (!parse_args(argc, argv, cfg)) {
        return 1;
    }

    // Convert decimal precision to bits for Arb
    const double LOG2_10 = 3.32192809488736234787; // log2(10)
    slong bits_prec = static_cast<slong>(std::ceil(cfg.precision * LOG2_10));
    if (bits_prec < 64) bits_prec = 64;

    std::cout << "[INFO] Loading zeros from: " << cfg.zeros_file << "\n";
    auto all_zero_strs = load_zeros_strings(cfg.zeros_file);
    if (all_zero_strs.empty()) {
        std::cerr << "[FATAL] No zeros loaded.\n";
        return 1;
    }

    long long start_index = cfg.start_idx;
    if (start_index < 0 ||
        start_index >= static_cast<long long>(all_zero_strs.size())) {
        std::cerr << "[FATAL] start_idx out of range.\n";
        return 1;
    }

    long long end_index = start_index + cfg.subset;
    if (end_index > static_cast<long long>(all_zero_strs.size())) {
        end_index = static_cast<long long>(all_zero_strs.size());
    }

    std::vector<std::string> sample_zero_strs;
    sample_zero_strs.reserve(end_index - start_index);
    for (long long k = start_index; k < end_index; ++k) {
        sample_zero_strs.push_back(all_zero_strs[static_cast<std::size_t>(k)]);
    }

    if (sample_zero_strs.empty()) {
        std::cerr << "[FATAL] Zero subset is empty; check --start_idx / --subset.\n";
        return 1;
    }

    double center_t = std::strtod(sample_zero_strs.front().c_str(), nullptr);
    std::cout << "[INFO] Using " << sample_zero_strs.size()
              << " zeros; first T ~ " << center_t << "\n";
    std::cout << "[INFO] Precision (decimal): " << cfg.precision
              << " -> bits: " << bits_prec << "\n";

    // Sigma grid centered at 0.5
    if (cfg.resolution < 2) cfg.resolution = 2;

    int    N    = cfg.resolution;
    double step = cfg.zoom / (N - 1); // total width = zoom
    int    mid  = N / 2;              // sigma[mid] = 0.5 exactly

    std::vector<double> sigmas;
    sigmas.reserve(N);
    for (int i = 0; i < N; ++i) {
        sigmas.push_back(0.5 + (static_cast<double>(i - mid) * step));
    }

    std::vector<double> energies(sigmas.size(), 0.0);

    std::cout << "[INFO] Computing energy landscape...\n";

    double min_energy = std::numeric_limits<double>::infinity();
    double min_sigma  = 0.0;

    const std::size_t num_sigmas = sigmas.size();
    const std::size_t num_zeros  = sample_zero_strs.size();

    // Arb temporaries reused across loops
    acb_t s, z;
    arb_t sigma_arb, T_arb, abs_z, energy_arb;
    acb_init(s);
    acb_init(z);
    arb_init(sigma_arb);
    arb_init(T_arb);
    arb_init(abs_z);
    arb_init(energy_arb);

    for (std::size_t i = 0; i < num_sigmas; ++i) {
        double sigma = sigmas[i];

        // energy_arb = 0
        arb_zero(energy_arb);

        // sigma_arb = sigma (from double; grid control, not sacred data)
        arb_set_d(sigma_arb, sigma);

        for (std::size_t k = 0; k < num_zeros; ++k) {
            const std::string& T_str = sample_zero_strs[k];

            // T_arb = T_k from full-precision string
            if (arb_set_str(T_arb, T_str.c_str(), bits_prec) != 0) {
                std::cerr << "[WARN] Failed to parse zero '" << T_str
                          << "' as arb; skipping.\n";
                continue;
            }

            // s = sigma + i*T
            acb_set_arb_arb(s, sigma_arb, T_arb);

            // z = zeta(s)
            acb_zeta(z, s, bits_prec);

            // abs_z = |z|
            acb_abs(abs_z, z, bits_prec);

            // energy_arb += abs_z
            arb_add(energy_arb, energy_arb, abs_z, bits_prec);
        }

        // Convert energy_arb's midpoint to double for logging + CSV.
        double energy_d = arf_get_d(arb_midref(energy_arb), ARF_RND_NEAR);
        energies[i] = energy_d;

        if (energy_d < min_energy) {
            min_energy = energy_d;
            min_sigma  = sigma;
        }

        if (num_sigmas >= 10 && (i % (num_sigmas / 10) == 0)) {
            std::cout << "  [progress] " << i << " / " << num_sigmas
                      << "  sigma=" << sigma
                      << "  E~" << energy_d << "\n";
        }
    }

    acb_clear(s);
    acb_clear(z);
    arb_clear(sigma_arb);
    arb_clear(T_arb);
    arb_clear(abs_z);
    arb_clear(energy_arb);
    flint_cleanup();

    // Write CSV
    std::ofstream out(cfg.output);
    if (!out) {
        std::cerr << "[FATAL] Could not open output: " << cfg.output << "\n";
        return 1;
    }

    out.setf(std::ios::fixed);
    out.precision(15);
    out << "# sigma,energy\n";
    for (std::size_t i = 0; i < sigmas.size(); ++i) {
        out << sigmas[i] << "," << energies[i] << "\n";
    }
    out.close();

    std::cout << "[INFO] CSV saved to: " << cfg.output << "\n";
    std::cout << "[INFO] Minimum energy (double) at sigma = "
              << min_sigma << "  E_min ~ " << min_energy << "\n";

    return 0;
}
