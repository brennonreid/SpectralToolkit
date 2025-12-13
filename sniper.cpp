/**
 * sniper.cpp
 * Robust "Ballistic Cone Descent"
 * FIX: Sanitizes Scientific Notation from Arb to prevent "Energy: 0" parsing errors.
 */

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <clocale>
#include <algorithm> // for std::max
#include <cctype>
#include "decimal.hpp"
#include "arb_interface.hpp"

// ------------------------------------------------------------
// Small helpers for Decimal
// ------------------------------------------------------------

static Decimal decimal_abs(const Decimal &x) {
    Decimal zero(0);
    return (x < zero) ? -x : x;
}

// ------------------------------------------------------------
// String sanitizers for Arb "[center +/- radius]" output
// ------------------------------------------------------------

// Trim whitespace in-place
static void trim_inplace(std::string &s) {
    size_t start = 0;
    while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start]))) {
        ++start;
    }
    size_t end = s.size();
    while (end > start && std::isspace(static_cast<unsigned char>(s[end - 1]))) {
        --end;
    }
    if (start == 0 && end == s.size()) return;
    s = s.substr(start, end - start);
}

// Convert a simple numeric string with optional scientific notation
// (e.g., "1.23e-4") into plain decimal. Assumes no brackets or "+/-".
static std::string sanitize_sci_number(const std::string &input) {
    std::string s = input;
    trim_inplace(s);
    if (s.empty()) return "0";

    // If there is no exponent, just return trimmed.
    size_t e_pos = s.find('e');
    if (e_pos == std::string::npos) e_pos = s.find('E');
    if (e_pos == std::string::npos) {
        return s;
    }

    // Parse sign
    bool neg = false;
    size_t i = 0;
    if (s[i] == '+' || s[i] == '-') {
        neg = (s[i] == '-');
        ++i;
    }

    // Parse mantissa digits and optional '.'
    std::string mantissa;
    int dec_pos = -1; // index in mantissa where '.' would be
    for (; i < s.size(); ++i) {
        char c = s[i];
        if (std::isdigit(static_cast<unsigned char>(c))) {
            mantissa.push_back(c);
        } else if (c == '.') {
            dec_pos = static_cast<int>(mantissa.size());
        } else {
            break;
        }
    }

    if (mantissa.empty()) {
        return "0";
    }
    if (dec_pos == -1) {
        dec_pos = static_cast<int>(mantissa.size());
    }

    // Parse exponent
    int exp_val = 0;
    if (i < s.size() && (s[i] == 'e' || s[i] == 'E')) {
        ++i;
        bool exp_neg = false;
        if (i < s.size() && (s[i] == '+' || s[i] == '-')) {
            exp_neg = (s[i] == '-');
            ++i;
        }
        int e = 0;
        for (; i < s.size(); ++i) {
            char c = s[i];
            if (!std::isdigit(static_cast<unsigned char>(c))) break;
            e = e * 10 + (c - '0');
        }
        exp_val = exp_neg ? -e : e;
    }

    int new_dec = dec_pos + exp_val;

    // Decimal point moves to the left of all digits
    if (new_dec <= 0) {
        std::string result;
        if (neg) result.push_back('-');
        result += "0.";
        result.append(static_cast<std::string::size_type>(-new_dec), '0');
        result += mantissa;

        // Trim trailing zeros
        size_t r_end = result.size();
        while (r_end > 0 && result[r_end - 1] == '0') {
            --r_end;
        }
        if (r_end > 0 && result[r_end - 1] == '.') {
            --r_end;
        }
        return result.substr(0, r_end);
    }

    // Decimal point moves to / stays to the right of all digits
    if (new_dec >= static_cast<int>(mantissa.size())) {
        std::string result;
        if (neg) result.push_back('-');
        result += mantissa;
        result.append(static_cast<std::string::size_type>(new_dec - mantissa.size()), '0');
        return result;
    }

    // Decimal point ends up somewhere inside mantissa
    std::string result;
    if (neg) result.push_back('-');
    result.append(mantissa.substr(0, static_cast<std::string::size_type>(new_dec)));
    result.push_back('.');
    result.append(mantissa.substr(static_cast<std::string::size_type>(new_dec)));

    // Trim trailing zeros
    size_t r_end = result.size();
    while (r_end > 0 && result[r_end - 1] == '0') {
        --r_end;
    }
    if (r_end > 0 && result[r_end - 1] == '.') {
        --r_end;
    }
    return result.substr(0, r_end);
}

// MAIN HELPER: Take an Arb "[center +/- radius]" string and return a
// plain decimal for the CENTER only, with scientific notation expanded.
std::string sanitize_arb_string(const std::string &raw) {
    if (raw.empty()) return "0";

    std::string inner = raw;

    // Strip outer brackets if present.
    size_t l = inner.find('[');
    size_t r = inner.rfind(']');
    if (l != std::string::npos && r != std::string::npos && r > l + 1) {
        inner = inner.substr(l + 1, r - l - 1);
    }

    // If there's a "+/-", keep only the center part to the left.
    size_t pm = inner.find("+/-");
    if (pm != std::string::npos) {
        inner = inner.substr(0, pm);
    }

    // Trim whitespace
    trim_inplace(inner);
    if (inner.empty()) return "0";

    // If the center itself uses scientific notation, expand it.
    size_t e_pos = inner.find('e');
    if (e_pos == std::string::npos) e_pos = inner.find('E');
    if (e_pos == std::string::npos) {
        // No exponent in the center; return as-is.
        return inner;
    }

    return sanitize_sci_number(inner);
}

// ------------------------------------------------------------
// Energy extraction
// ------------------------------------------------------------

Decimal get_energy_as_decimal(ArbZetaContext &ctx, int digits_prec) {
    // Ask Arb for a bit more than we need so the center is not rounded too aggressively.
    std::string raw = ctx.get_energy_str(digits_prec + 10);

    std::string sanitized = sanitize_arb_string(raw);
    if (sanitized.empty()) sanitized = "0";

    // DEBUG: only print when Arb uses scientific notation (usually near very small values)
    if (raw.find('e') != std::string::npos || raw.find('E') != std::string::npos) {
        std::cout << " [DEBUG] Arb energy (sci): [" << raw
                  << "] -> sanitized: " << sanitized << "\n";
    }

    return decimal_from_string(sanitized);
}

// ------------------------------------------------------------
// ROBUST Dual-Scale Slope Calculation
// ------------------------------------------------------------

Decimal get_robust_slope(ArbZetaContext &ctx, const Decimal &T, const Decimal &h_base, int prec_digits) {
    Decimal sigma = decimal_from_string("0.5");
    Decimal zero(0);

    // 1. Measure Base Energy
    ctx.compute(T, sigma);
    Decimal E0 = get_energy_as_decimal(ctx, prec_digits);

    // 2. Calculate Narrow Slope (the "Sniper"): central difference
    Decimal h_narrow = h_base;
    if (h_narrow == zero) {
        h_narrow = decimal_from_string("1e-6");
    }

    ctx.compute(T + h_narrow, sigma);
    Decimal E_p_n = get_energy_as_decimal(ctx, prec_digits);

    ctx.compute(T - h_narrow, sigma);
    Decimal E_m_n = get_energy_as_decimal(ctx, prec_digits);

    Decimal slope_narrow = (E_p_n - E_m_n) / (Decimal(2) * h_narrow);

    // 3. Calculate Wide Slope (the "Cone"): central difference at larger h
    Decimal h_wide = h_narrow * Decimal(10);
    Decimal h_wide_min = decimal_from_string("1e-4");
    if (decimal_abs(h_wide) < h_wide_min) {
        h_wide = (h_wide < zero) ? -h_wide_min : h_wide_min;
    }

    ctx.compute(T + h_wide, sigma);
    Decimal E_p_w = get_energy_as_decimal(ctx, prec_digits);

    ctx.compute(T - h_wide, sigma);
    Decimal E_m_w = get_energy_as_decimal(ctx, prec_digits);

    Decimal slope_wide = (E_p_w - E_m_w) / (Decimal(2) * h_wide);

    // 4. CONSENSUS LOGIC

    // Case A: Trapped in noise? (Narrow is ~0, Wide is moving)
    if (slope_narrow == zero && slope_wide != zero) {
        std::cout << " [AUTO-CORRECT] Escaping precision trap using Wide Slope.\n";
        return slope_wide;
    }

    // Case B: Overshoot risk? (Opposite signs -> take the safer narrow slope)
    if ((slope_wide > zero && slope_narrow < zero) ||
        (slope_wide < zero && slope_narrow > zero)) {
        return slope_narrow;
    }

    // Default: trust narrow slope
    return slope_narrow;
}

// ------------------------------------------------------------
// Locked digits: count shared leading digits between successive T
// (ignoring the decimal point).
// ------------------------------------------------------------

int count_locked_digits(const Decimal &prev_T, const Decimal &curr_T, int target_dps) {
    std::string s_prev = decimal_to_string(prev_T, target_dps);
    std::string s_curr = decimal_to_string(curr_T, target_dps);

    int locked = 0;
    size_t i_prev = 0;
    size_t i_curr = 0;
    while (i_prev < s_prev.size() && i_curr < s_curr.size()) {
        char c_prev = s_prev[i_prev];
        char c_curr = s_curr[i_curr];

        // Skip decimal points in both strings
        if (c_prev == '.') {
            ++i_prev;
            continue;
        }
        if (c_curr == '.') {
            ++i_curr;
            continue;
        }

        if (c_prev != c_curr) break;

        ++locked;
        ++i_prev;
        ++i_curr;
    }
    return locked;
}

// ------------------------------------------------------------
// Main sniper loop
// ------------------------------------------------------------

void run_sniper(const std::string &start_T_str, int target_dps) {
    std::setlocale(LC_ALL, "C");

    int internal_bits = static_cast<int>(target_dps * 3.33) + 256;
    ArbZetaContext ctx(internal_bits);

    // Parse input string directly (no double -> string round-trip).
    Decimal T = decimal_from_string(start_T_str);
    Decimal sigma = decimal_from_string("0.5");

    std::cout << "==============================================================\n";
    std::cout << " SNIPER: Ballistic Cone Descent (Value Preserved)\n";
    std::cout << " Target Precision: " << target_dps << " digits\n";
    std::cout << " Start T: " << decimal_to_string(T, target_dps) << "\n";
    std::cout << "==============================================================\n";

    Decimal zero(0);
    Decimal energy_floor = decimal_from_string("1e-60");
    Decimal min_step = decimal_from_string("1e-8");
    Decimal max_step = decimal_from_string("0.1");

    Decimal prev_T = T;
    int locked = 0;

    int step = 0;
    const int max_steps = 100;

    while (step < max_steps) {
        ++step;

        // 1. Current energy
        ctx.compute(T, sigma);
        Decimal E_val = get_energy_as_decimal(ctx, target_dps);
        std::string e_str = decimal_to_string(E_val, 20);

        if (E_val == zero) {
            std::string raw_s = ctx.get_energy_str(target_dps + 10);
            std::cout << " [DEBUG] Step " << step
                      << " Energy is 0. Raw Arb string: [" << raw_s << "]\n";
        }

        // Termination check: tiny energy (but avoid early bail on step 1).
        if (E_val < energy_floor && step > 1) {
            std::cout << "\n[ZERO LOCK] Energy floor reached (1e-60).\n";
            break;
        }

        // 2. Step size h (base for narrow slope)
        Decimal h_base;
        if (E_val == zero) {
            h_base = min_step;
        } else {
            h_base = E_val / Decimal(1000);
            if (h_base < min_step) h_base = min_step;
            if (h_base > max_step) h_base = max_step;
        }

        // 3. Compute robust slope
        Decimal slope = get_robust_slope(ctx, T, h_base, target_dps);

        // Slope vanish safeguard
        if (slope == zero) {
            std::cout << "[WARNING] Slope vanish. Stuck.\n";
            std::cout << " [AUTO-KICK] Attempting to nudge +0.001\n";
            T = T + decimal_from_string("0.001");
            prev_T = T;
            step = 0; // restart counting from next iteration
            continue;
        }

        // 4. Compute shift (ballistic cone step).
        // Use normalized slope sign to avoid insane jumps.
        Decimal sign = (slope > zero) ? Decimal(1) : Decimal(-1);
        Decimal mag = E_val;
        if (mag == zero) {
            mag = min_step;
        }
        Decimal shift = -sign * mag;

        // Clamp shift to sane bounds
        if (decimal_abs(shift) < min_step) {
            shift = (shift < zero ? -min_step : min_step);
        }
        if (decimal_abs(shift) > max_step) {
            shift = (shift < zero ? -max_step : max_step);
        }

        Decimal T_new = T + shift;

        // 5. Locked digits (relative to previous T)
        locked = count_locked_digits(prev_T, T_new, target_dps);

        // 6. Report step
        std::cout << " Step " << step << ":\n";
        std::cout << "   Energy : " << e_str << "\n";
        std::cout << "   Slope  : " << decimal_to_string(slope, 10) << "\n";
        std::cout << "   Shift  : " << decimal_to_string(shift, 20) << "\n";
        std::cout << "   Locked : " << locked << " digits\n";

        prev_T = T;
        T = T_new;

        if (locked >= target_dps) {
            std::cout << "\n[TARGET PRECISION REACHED]\n";
            break;
        }
    }

    std::cout << "\n--------------------------------------------------------------\n";
    std::cout << " FINAL SNIPE RESULT:\n";
    std::cout << " " << decimal_to_string(T, target_dps) << "\n";
    std::cout << "--------------------------------------------------------------\n";
}

// ------------------------------------------------------------
// CLI entry point
// ------------------------------------------------------------

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: sniper <T_approx> [precision=100]\n";
        return 1;
    }
    std::string start_T = argv[1];
    int precision = 100;
    if (argc > 2) precision = std::stoi(argv[2]);
    try {
        run_sniper(start_T, precision);
    } catch (const std::exception& e) {
        std::cerr << "CRITICAL ERROR: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
