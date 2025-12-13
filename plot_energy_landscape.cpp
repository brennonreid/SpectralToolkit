#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <chrono> // Standard C++ timing

#include "decimal.hpp"      // Provides Decimal (mpq_class) and utilities
#include "arb_interface.hpp" // Provides ArbZetaContext and FLINT/Arb interface

// --- PLACEHOLDER: Load Zeros from file (CRITICAL USER SECTION) ---
// NOTE: This function must correctly load your zeros.txt file. 
std::vector<Decimal> load_zeros_robust(const std::string& path, int start_idx, int subset) {
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "[FATAL] Could not open zeros file: " << path << std::endl;
        return {};
    }

    std::vector<Decimal> all_zeros;
    std::string line;
    while (std::getline(file, line)) {
        std::string t_str;
        std::size_t last_space = line.find_last_of(" \t");
        if (last_space != std::string::npos) {
            t_str = line.substr(last_space + 1);
        } else {
            t_str = line;
        }
        
        try {
            Decimal t_val = decimal_from_string(t_str);
            all_zeros.push_back(t_val);
        } catch (...) {
            // Ignore junk lines
        }
    }

    if (start_idx >= (int)all_zeros.size()) {
        std::cerr << "[FATAL] Start index out of bounds (" << start_idx << " >= " << all_zeros.size() << ")\n";
        return {};
    }
    int end_idx = std::min(start_idx + subset, (int)all_zeros.size());
    
    std::vector<Decimal> sample_zeros(all_zeros.begin() + start_idx, all_zeros.begin() + end_idx);
    
    return sample_zeros;
}


// --- CORE SINGLE-THREADED COMPUTATION ---
void compute_energy_landscape(
    const std::vector<Decimal>& sigmas,
    const std::vector<Decimal>& zeros_T,
    slong prec,
    std::vector<long double>& energies
) {
    energies.assign(sigmas.size(), 0.0L);
    size_t num_sigmas = sigmas.size();

    // Context is created once and reused (safe in single-threaded mode)
    ArbZetaContext context(prec); 

    // Loop is now strictly serial (single-threaded)
    for (size_t i = 0; i < num_sigmas; ++i) {
        long double total_energy = 0.0L;

        // Sum the energy contribution from all T values (zeros)
        for (const auto& T : zeros_T) {
            context.compute(T, sigmas[i]);
            total_energy += context.get_energy_double();
        }
        
        energies[i] = total_energy;
    }
}


// --- MAIN EXECUTION ---
int main(int argc, char* argv[]) {
    // --- Configuration Defaults ---
    const int PREC_BITS = 128;
    std::string output_file = "zeta_landscape_cpp_data.txt";
    std::string zeros_file = "";
    double zoom_width = 0.01;
    int resolution = 201;
    int start_idx = 50000;
    int subset = 50;
    slong precision_bits = PREC_BITS;

    // --- Argument Parsing ---
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <zeros_file> [options]\n";
        return 1;
    }
    
    zeros_file = argv[1];

    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--zoom" && i + 1 < argc) zoom_width = std::stod(argv[++i]);
        else if (arg == "--res" && i + 1 < argc) resolution = std::stoi(argv[++i]);
        else if (arg == "--start" && i + 1 < argc) start_idx = std::stoi(argv[++i]);
        else if (arg == "--subset" && i + 1 < argc) subset = std::stoi(argv[++i]);
        else if (arg == "--prec" && i + 1 < argc) precision_bits = std::stol(argv[++i]);
        else if (arg == "--output" && i + 1 < argc) output_file = argv[++i];
    }
    
    // --- 1. Load Zeros and Setup ---
    std::vector<Decimal> sample_zeros = load_zeros_robust(zeros_file, start_idx, subset);
    if (sample_zeros.empty()) return 1;

    std::cout << "[INFO] Running in **SINGLE-THREADED MODE**.\n";
    std::cout << "[INFO] Using " << sample_zeros.size() << " zeros; precision: " << precision_bits << " bits." << std::endl;
    std::cout << "[INFO] Computing landscape on " << resolution << " sigma points..." << std::endl;

    // Sigma grid
    Decimal start_s(0.5 - (zoom_width / 2.0));
    Decimal end_s(0.5 + (zoom_width / 2.0));
    Decimal step = (end_s - start_s) / Decimal(resolution - 1);
    
    std::vector<Decimal> sigmas = frange_decimal(start_s, end_s, step);
    if (sigmas.size() > (size_t)resolution) sigmas.resize(resolution);

    // --- 2. Compute Energy ---
    std::vector<long double> energies;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    compute_energy_landscape(sigmas, sample_zeros, precision_bits, energies);
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

    std::cout << "[INFO] Finished computation in " << elapsed << " seconds." << std::endl;

    // --- 3. Find Minimum (for reporting) ---
    auto min_it = std::min_element(energies.begin(), energies.end());
    size_t min_idx = std::distance(energies.begin(), min_it);
    
    std::string min_sigma_str = decimal_to_string(sigmas[min_idx], 10);
    
    std::cout << "[INFO] Minimum energy (long double) = " << *min_it << std::endl;
    std::cout << "[INFO] True minimum found at sigma = " << min_sigma_str << std::endl;

    // --- 4. Output Data to File ---
    std::ofstream outfile(output_file);
    if (!outfile.is_open()) {
        std::cerr << "[FATAL] Could not write output file: " << output_file << std::endl;
        return 1;
    }

    outfile << "# Sigma \t Energy_Magnitude\n";
    outfile << std::fixed << std::setprecision(15);
    
    for (size_t i = 0; i < sigmas.size(); ++i) {
        outfile << decimal_to_string(sigmas[i], 15) << "\t" << energies[i] << "\n";
    }
    outfile.close();
    std::cout << "[INFO] Data saved to: " << output_file << std::endl;

    return 0;
}