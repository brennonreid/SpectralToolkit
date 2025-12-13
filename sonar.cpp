/**
 * sonar.cpp
 * Updated for Multi-Target "Hunter-Killer" Logic
 * Fixed: Decimal("1e-6") crash in deduplication
 */

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <future>
#include <thread>
#include <mutex>
#include <limits>
#include <getopt.h>

#include "decimal.hpp"
#include "arb_interface.hpp"

// ------------------------------------------------------------
// Data Structures
// ------------------------------------------------------------

struct Result {
    Decimal T;
    Decimal sigma_min;
    long double E_min_dbl;
    std::string E_min_str; 
    std::string classification; 
    
    // For carrying candidates up the stack
    std::vector<Result> multi_candidates; 
};

// ------------------------------------------------------------
// Helpers
// ------------------------------------------------------------

int count_decimals(const Decimal &T) {
    std::string s = decimal_to_string(T, 50); 
    return count_decimals_str(s);
}

Decimal truncate_T(const Decimal &T, int digits) {
    return truncate_decimal(T, digits);
}

// Deduplicate results that have converged to the same Zero
std::vector<Result> deduplicate_results(std::vector<Result> &raw) {
    if (raw.empty()) return {};

    // Sort by T
    std::sort(raw.begin(), raw.end(), [](const Result &a, const Result &b) {
        return a.T < b.T;
    });

    std::vector<Result> unique;
    unique.push_back(raw[0]);

    // Safety threshold
    Decimal safe_threshold = decimal_from_string("1e-6");

    for (size_t i = 1; i < raw.size(); ++i) {
        // If T is very close (within 1e-6), assume same zero and keep the one with lower Energy
        Decimal diff = raw[i].T - unique.back().T;
        if (diff < 0) diff = -diff;

        if (diff < safe_threshold) {
            // Merge: keep the one with better energy
            if (raw[i].E_min_dbl < unique.back().E_min_dbl) {
                unique.back() = raw[i];
            }
        } else {
            unique.push_back(raw[i]);
        }
    }
    return unique;
}

// ------------------------------------------------------------
// Core Logic
// ------------------------------------------------------------

Result scan_single_T(ArbZetaContext &ctx, 
                     const Decimal &T, 
                     const Decimal &sigma_left, 
                     const Decimal &sigma_right, 
                     int sigma_points) {
    
    if (sigma_points < 2) sigma_points = 2;
    Decimal step = (sigma_right - sigma_left) / Decimal(sigma_points - 1);
    
    Decimal sigma = sigma_left;
    Decimal best_sigma = sigma_left;
    long double min_E = std::numeric_limits<long double>::max();
    
    for (int i = 0; i < sigma_points; ++i) {
        long double E = ctx.energy_at_sigma(T, sigma);
        if (E < min_E) {
            min_E = E;
            best_sigma = sigma;
        }
        sigma += step;
    }
    
    // Compute high-precision string for potential reporting
    ctx.compute(T, best_sigma);
    std::string e_str = ctx.get_energy_str(25);

    Result res;
    res.T = T;
    res.sigma_min = best_sigma;
    res.E_min_dbl = min_E;
    res.E_min_str = e_str;
    res.classification = "NONE";
    return res;
}

std::string classify_result(const Result &res, 
                            const Decimal &sigma_target, 
                            const Decimal &sigma_tol, 
                            const Decimal &energy_tol) {
    
    Decimal sig_diff = res.sigma_min - sigma_target;
    if (sig_diff < 0) sig_diff = -sig_diff;
    
    bool sigma_ok = (sig_diff <= sigma_tol);
    double e_tol_d = std::stod(decimal_to_string(energy_tol, 10));
    bool energy_ok = (res.E_min_dbl <= e_tol_d);

    if (sigma_ok && energy_ok) return "LOCK";
    return "NONE";
}

// ------------------------------------------------------------
// Spider Logic
// ------------------------------------------------------------

Result spider_scan_T(ArbZetaContext &ctx,
                     const Decimal &T,
                     Decimal sigma_left,
                     Decimal sigma_right,
                     int sigma_points,
                     int max_levels,
                     long double target_energy = -1.0) {
    
    Result best_res;
    best_res.E_min_dbl = std::numeric_limits<long double>::max();

    Decimal cur_left = sigma_left;
    Decimal cur_right = sigma_right;
    
    Decimal limit = decimal_from_string("1e-20");

    for (int level = 0; level < max_levels; ++level) {
        Result res = scan_single_T(ctx, T, cur_left, cur_right, sigma_points);

        if (res.E_min_dbl < best_res.E_min_dbl) {
            best_res = res;
        }

        if (target_energy > 0.0 && res.E_min_dbl <= target_energy) {
            break;
        }

        Decimal width = cur_right - cur_left;
        if (width <= limit) break;

        Decimal new_half = width / 4;
        if (new_half <= limit) break;

        Decimal new_left = res.sigma_min - new_half;
        Decimal new_right = res.sigma_min + new_half;

        if (new_left < sigma_left) new_left = sigma_left;
        if (new_right > sigma_right) new_right = sigma_right;

        if ((new_right - new_left) <= limit) break;

        cur_left = new_left;
        cur_right = new_right;
    }
    return best_res;
}

std::pair<Decimal, Result> refine_T_band(
    ArbZetaContext &ctx,
    const Decimal &T_center,
    int digit_depth,
    const Decimal &sigma_left,
    const Decimal &sigma_right,
    int sigma_points,
    int sigma_levels,
    long double target_energy,
    const Decimal &T_step_min,
    int prefix_decimals,
    const Decimal &T_floor
) {
    int d = digit_depth;
    Decimal prefix_d = truncate_T(T_center, d);

    mpz_class p10 = pow10_z(d + 2);
    Decimal fine_step(1, p10);

    if (fine_step < T_step_min) {
        fine_step = T_step_min;
    }

    Decimal band_left, band_right;
    if (d == 0) {
        band_left = prefix_d - fine_step;
        band_right = prefix_d + 1;
    } else {
        Decimal coarse_unit(1, pow10_z(d));
        band_left = prefix_d - fine_step;
        band_right = prefix_d + coarse_unit;
    }

    std::cout << "[DIGIT-BAND] T_center=" << decimal_to_string(T_center, 40)
              << " depth=" << d 
              << " band=[" << decimal_to_string(band_left, 10) << ", " << decimal_to_string(band_right, 10) << "]"
              << std::endl;

    Result best_res;
    best_res.E_min_dbl = std::numeric_limits<long double>::max();
    Decimal best_T_raw = T_center;
    bool found_valid = false;

    struct Sample {
        Decimal T;
        Decimal sigma;
        long double E;
        bool prefix_ok;
    };
    std::vector<Sample> samples;

    Decimal T = band_left;
    Decimal stop_T = band_right + (fine_step / 2);

    while (T <= stop_T) {
        Result res_T = spider_scan_T(ctx, T, sigma_left, sigma_right, sigma_points, sigma_levels, target_energy);
        
        Decimal T_prefix = truncate_T(T, (prefix_decimals > 0 ? prefix_decimals : 0));
        bool prefix_ok = (T_prefix == T_floor);

        // Optional: reduce verbosity if needed
        std::cout << "    [digit-sample] T=" << decimal_to_string(T, 40)
                  << " sig=" << decimal_to_string(res_T.sigma_min, 10)
                  << " E=" << res_T.E_min_dbl << (prefix_ok ? "" : " [SKIP]") << std::endl;

        samples.push_back({T, res_T.sigma_min, res_T.E_min_dbl, prefix_ok});

        if (prefix_ok) {
            if (!found_valid || res_T.E_min_dbl < best_res.E_min_dbl) {
                best_res = res_T;
                best_T_raw = T;
                found_valid = true;
            }
        }

        T += fine_step;
    }

    if (!found_valid) {
        Result dummy;
        dummy.T = T_center;
        dummy.E_min_dbl = 1e9;
        dummy.classification = "NONE";
        return {T_center, dummy};
    }

    Decimal best_T_coarse = truncate_T(best_T_raw, d + 1);

    std::vector<Result> multi_candidates;
    for (size_t i = 1; i < samples.size() - 1; ++i) {
        if (!samples[i].prefix_ok) continue;

        double E_prev = samples[i-1].E;
        double E_curr = samples[i].E;
        double E_next = samples[i+1].E;

        // Detect local minimum
        if (E_curr <= E_prev && E_curr <= E_next) {
            Result cand_res;
            cand_res.T = truncate_T(samples[i].T, d + 1);
            cand_res.sigma_min = samples[i].sigma;
            cand_res.E_min_dbl = E_curr;
            multi_candidates.push_back(cand_res);
        }
    }

    if (!multi_candidates.empty()) {
        std::cout << "    [digit-candidates] found " << multi_candidates.size() << std::endl;
    }

    best_res.T = best_T_coarse;
    best_res.multi_candidates = multi_candidates;

    return {best_T_coarse, best_res};
}

// ------------------------------------------------------------
// Multi-Target Hunter-Killer Spider
// ------------------------------------------------------------
Result spider_search_2d(ArbZetaContext &ctx,
                        const Decimal &T0,
                        const Decimal &sigma_left,
                        const Decimal &sigma_right,
                        int sigma_points,
                        int sigma_levels,
                        long double target_energy,
                        int vertical_levels,
                        Decimal T_step_min) {
    
    int prefix_decimals = count_decimals(T0);
    Decimal T_floor = truncate_T(T0, prefix_decimals);
    int digit_depth = prefix_decimals;

    // ACTIVE LIST: Nodes still being refined
    std::vector<Decimal> active_centers;
    active_centers.push_back(T0);

    // COMPLETED LIST: Nodes that hit target or timed out
    std::vector<Result> final_locks;

    for (int v_level = 0; v_level < vertical_levels; ++v_level) {
        if (active_centers.empty()) break;

        mpz_class p10 = pow10_z(digit_depth + 2);
        Decimal fine_step0(1, p10);
        
        if (fine_step0 < T_step_min) {
            std::cout << "[SPIDER] Stopping: fine step too small." << std::endl;
            // Capture current status of active nodes before exiting
            for (const auto &t : active_centers) {
                Result r = spider_scan_T(ctx, t, sigma_left, sigma_right, sigma_points, sigma_levels, target_energy);
                final_locks.push_back(r);
            }
            break;
        }

        std::cout << "[SPIDER] v_level=" << v_level << " digit_depth=" << digit_depth 
                  << " active_nodes=" << active_centers.size() 
                  << " finished_locks=" << final_locks.size() << std::endl;

        std::vector<Decimal> next_active_raw;

        for (size_t idx = 0; idx < active_centers.size(); ++idx) {
            Decimal tc = active_centers[idx];
            
            auto pair = refine_T_band(ctx, tc, digit_depth, sigma_left, sigma_right, 
                                      sigma_points, sigma_levels, target_energy, 
                                      T_step_min, prefix_decimals, T_floor);
            
            Decimal primary_T = pair.first;
            Result res_T = pair.second;

            // CHECK PRIMARY
            if (target_energy > 0.0 && res_T.E_min_dbl <= target_energy) {
                std::cout << "    [TARGET HIT] T=" << decimal_to_string(primary_T, 40) << " E=" << res_T.E_min_dbl << std::endl;
                final_locks.push_back(res_T);
            } else {
                next_active_raw.push_back(primary_T);
            }

            // CHECK CANDIDATES
            for (const auto &cand : res_T.multi_candidates) {
                if (target_energy > 0.0 && cand.E_min_dbl <= target_energy) {
                     std::cout << "    [TARGET HIT-CANDIDATE] T=" << decimal_to_string(cand.T, 40) << " E=" << cand.E_min_dbl << std::endl;
                     final_locks.push_back(cand);
                } else {
                    next_active_raw.push_back(cand.T);
                }
            }
        }

        // Deduplicate active centers for next round
        std::sort(next_active_raw.begin(), next_active_raw.end());
        next_active_raw.erase(std::unique(next_active_raw.begin(), next_active_raw.end()), next_active_raw.end());
        active_centers = next_active_raw;

        digit_depth++;
    }

    // Any remaining active centers that didn't hit target? Save them now.
    for (const auto &t : active_centers) {
        Result r = spider_scan_T(ctx, t, sigma_left, sigma_right, sigma_points, sigma_levels, target_energy);
        final_locks.push_back(r);
    }

    // Clean up Final Locks (deduplicate)
    final_locks = deduplicate_results(final_locks);

    // Sort to find the "best" one to return as primary
    if (final_locks.empty()) {
        Result dummy; dummy.T = T0; dummy.E_min_dbl = 1e9;
        return dummy;
    }

    std::sort(final_locks.begin(), final_locks.end(), [](const Result &a, const Result &b){
        return a.E_min_dbl < b.E_min_dbl;
    });

    Result primary = final_locks[0];
    
    // Put the rest in multi_candidates
    for (size_t i = 1; i < final_locks.size(); ++i) {
        primary.multi_candidates.push_back(final_locks[i]);
    }

    return primary;
}


// ------------------------------------------------------------
// Main
// ------------------------------------------------------------

void print_help() {
    std::cout << "Usage: sonar --mode [single|sweep|spider] --T <val> [options]\n";
}

int main(int argc, char** argv) {
    if (argc < 2) {
        print_help();
        return 1;
    }

    // Defaults
    std::string mode = "single";
    std::string T_str = "";
    std::string T_start_str, T_end_str, T_step_str = "0.1";
    int dps = 40;
    std::string sigma_left_str = "0.4", sigma_right_str = "0.6";
    int sigma_points = 41;
    std::string sigma_tol_str = "0.01", energy_tol_str = "0.1";
    std::string target_energy_str = "";
    int workers = 1;
    std::string csv_path;
    
    int spider_levels = 6;
    int vertical_levels = 12;
    std::string T_step_min_str = "1e-12";

    // Manual Arg Parsing
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--mode" && i+1 < argc) mode = argv[++i];
        else if (arg == "--T" && i+1 < argc) T_str = argv[++i];
        else if (arg == "--T_start" && i+1 < argc) T_start_str = argv[++i];
        else if (arg == "--T_end" && i+1 < argc) T_end_str = argv[++i];
        else if (arg == "--T_step" && i+1 < argc) T_step_str = argv[++i];
        else if (arg == "--dps" && i+1 < argc) dps = std::stoi(argv[++i]);
        else if (arg == "--sigma_left" && i+1 < argc) sigma_left_str = argv[++i];
        else if (arg == "--sigma_right" && i+1 < argc) sigma_right_str = argv[++i];
        else if (arg == "--sigma_points" && i+1 < argc) sigma_points = std::stoi(argv[++i]);
        else if (arg == "--sigma_tol" && i+1 < argc) sigma_tol_str = argv[++i];
        else if (arg == "--energy_tol" && i+1 < argc) energy_tol_str = argv[++i];
        else if (arg == "--target_energy" && i+1 < argc) target_energy_str = argv[++i];
        else if (arg == "--workers" && i+1 < argc) workers = std::stoi(argv[++i]);
        else if (arg == "--csv" && i+1 < argc) csv_path = argv[++i];
        else if (arg == "--spider_levels" && i+1 < argc) spider_levels = std::stoi(argv[++i]);
        else if (arg == "--vertical_levels" && i+1 < argc) vertical_levels = std::stoi(argv[++i]);
        else if (arg == "--T_step_min" && i+1 < argc) T_step_min_str = argv[++i];
    }

    slong prec_bits = (slong)(dps * 3.33 + 20);
    ArbZetaContext ctx(prec_bits);

    Decimal sigma_left = decimal_from_string(sigma_left_str);
    Decimal sigma_right = decimal_from_string(sigma_right_str);
    Decimal sigma_tol = decimal_from_string(sigma_tol_str);
    Decimal energy_tol = decimal_from_string(energy_tol_str);
    
    long double target_energy = -1.0;
    if (!target_energy_str.empty()) {
        target_energy = std::stod(target_energy_str);
    }

    std::vector<Result> results;
    std::mutex results_mutex;

    if (mode == "single") {
        if (T_str.empty()) { std::cerr << "--T required" << std::endl; return 1; }
        Decimal T = decimal_from_string(T_str);
        Result res = scan_single_T(ctx, T, sigma_left, sigma_right, sigma_points);
        res.classification = classify_result(res, decimal_from_string("0.5"), sigma_tol, energy_tol);
        results.push_back(res);

        std::cout << "[SINGLE] T = " << decimal_to_string(res.T, 20) << "\n"
                  << "  sigma_min = " << decimal_to_string(res.sigma_min, 20) << "\n"
                  << "  E_min     = " << res.E_min_str << "\n"
                  << "  class     = " << res.classification << std::endl;
    }
    else if (mode == "spider") {
        if (T_str.empty()) { std::cerr << "--T required" << std::endl; return 1; }
        Decimal T0 = decimal_from_string(T_str);
        Decimal T_step_min = decimal_from_string(T_step_min_str);
        
        Result res = spider_search_2d(ctx, T0, sigma_left, sigma_right, sigma_points, 
                                      spider_levels, target_energy, vertical_levels, T_step_min);
        
        res.classification = classify_result(res, decimal_from_string("0.5"), sigma_tol, energy_tol);
        results.push_back(res);

        // Sort secondary candidates by T for clean reporting
        std::sort(res.multi_candidates.begin(), res.multi_candidates.end(), [](const Result &a, const Result &b){
            return a.T < b.T;
        });

        std::cout << "\n========================================\n"
                  << " SPIDER REPORT\n"
                  << "========================================\n";
        
        // Print the winner (lowest energy)
        std::cout << "PRIMARY LOCK:\n"
                  << "  T         = " << decimal_to_string(res.T, 40) << "\n"
                  << "  E_min     = " << res.E_min_str << "\n"
                  << "  Class     = " << res.classification << "\n\n";

        // Print others
        if (!res.multi_candidates.empty()) {
            std::cout << "SECONDARY LOCKS (" << res.multi_candidates.size() << " found):\n";
            for (auto &cand : res.multi_candidates) {
                std::string cls = classify_result(cand, decimal_from_string("0.5"), sigma_tol, energy_tol);
                
                // If we don't have the high-prec string computed yet, compute it now
                if (cand.E_min_str.empty()) {
                     ctx.compute(cand.T, cand.sigma_min);
                     cand.E_min_str = ctx.get_energy_str(25);
                }

                std::cout << "  [Target]\n"
                          << "    T     = " << decimal_to_string(cand.T, 40) << "\n"
                          << "    E_min = " << cand.E_min_str << "\n" // Use high prec
                          << "    Class = " << cls << "\n";
            }
        }
        std::cout << "========================================\n";
    }
else if (mode == "sweep") {
    if (T_start_str.empty() || T_end_str.empty()) {
        std::cerr << "T_start/end required" << std::endl;
        return 1;
    }
    if (T_step_str.empty()) {
        std::cerr << "T_step required for sweep mode" << std::endl;
        return 1;
    }

    Decimal start = decimal_from_string(T_start_str);
    Decimal end   = decimal_from_string(T_end_str);
    Decimal raw_step = decimal_from_string(T_step_str);
    Decimal zero(0);

    if (raw_step == zero) {
        std::cerr << "T_step must be non-zero" << std::endl;
        return 1;
    }

    // Use the *magnitude* of the step as the stride.
    // Direction is derived from start vs end so we can sweep either way.
    Decimal step_mag = decimal_abs(raw_step);
    Decimal dir = (end >= start) ? Decimal(1) : Decimal(-1);
    Decimal step = dir * step_mag;

    // Build T_values explicitly: start, start+step, ..., until we cross end.
    std::vector<Decimal> T_values;
    T_values.reserve(1024); // optional, will grow as needed

    Decimal T = start;
    std::size_t safety_counter = 0;
    const std::size_t safety_limit = 100000000; // hard guard against infinite loop

    if (dir > zero) {
        // Increasing sweep: start <= end, step > 0
        while (T <= end) {
            T_values.push_back(T);
            T = T + step;
            if (++safety_counter > safety_limit) {
                std::cerr << "ERROR: sweep loop exceeded safety limit (increasing)" << std::endl;
                return 1;
            }
        }
    } else {
        // Decreasing sweep: start >= end, step < 0
        while (T >= end) {
            T_values.push_back(T);
            T = T + step;
            if (++safety_counter > safety_limit) {
                std::cerr << "ERROR: sweep loop exceeded safety limit (decreasing)" << std::endl;
                return 1;
            }
        }
    }

    std::cout << "Sweep: " << T_values.size() << " points." << std::endl;

    if (workers < 1) workers = 1;
    
    auto worker_func = [&](size_t idx_start, size_t idx_end) {
        ArbZetaContext local_ctx(prec_bits);
        
        std::vector<Result> local_results;
        local_results.reserve(idx_end - idx_start);

        for (size_t i = idx_start; i < idx_end; ++i) {
            Result r = scan_single_T(local_ctx, T_values[i], sigma_left, sigma_right, sigma_points);
            r.classification = classify_result(r, decimal_from_string("0.5"), sigma_tol, energy_tol);
            local_results.push_back(r);
        }
        
        std::lock_guard<std::mutex> lock(results_mutex);
        results.insert(results.end(), local_results.begin(), local_results.end());
    };

    std::vector<std::future<void>> futures;
    size_t total = T_values.size();
    size_t chunk_size = (total + workers - 1) / workers;

    for (int w = 0; w < workers; ++w) {
        size_t idx_start = w * chunk_size;
        size_t idx_end   = std::min(idx_start + chunk_size, total);
        if (idx_start >= total) break;
        
        futures.push_back(std::async(std::launch::async, worker_func, idx_start, idx_end));
    }

    for (auto &f : futures) {
        f.get();
    }
    
    // Always sort by T so output is ordered, regardless of sweep direction.
    std::sort(results.begin(), results.end(), [](const Result &a, const Result &b){
        return a.T < b.T;
    });

    for (const auto &r : results) {
        std::cout << "T = " << decimal_to_string(r.T, 10) 
                  << "  sigma_min = " << decimal_to_string(r.sigma_min, 5)
                  << "  E_min = " << r.E_min_dbl
                  << "  class = " << r.classification << std::endl;
    }
}


    if (!csv_path.empty() && !results.empty()) {
        std::ofstream f(csv_path);
        f << "T,sigma_min,E_min,class\n";
        for (const auto &r : results) {
            f << decimal_to_string(r.T, 40) << ","
              << decimal_to_string(r.sigma_min, 20) << ","
              << r.E_min_str << "," 
              << r.classification << "\n";
        }
    }

    return 0;
}