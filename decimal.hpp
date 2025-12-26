#ifndef DECIMAL_HPP
#define DECIMAL_HPP

#include <gmpxx.h>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cctype>

// Rational decimal type
using Decimal = mpq_class;

// 10^n as big integer
inline mpz_class pow10_z(long n)
{
    mpz_class r = 1;
    if (n == 0) return 1;
    if (n < 0) return 0; 

    // Use efficient GMP pow_ui
    mpz_ui_pow_ui(r.get_mpz_t(), 10, static_cast<unsigned long>(n));
    return r;
}

// Trim whitespace
inline std::string trim(const std::string &s_in)
{
    std::string s = s_in;
    auto not_space = [](unsigned char c) { return !std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
    s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
    return s;
}

static Decimal decimal_abs(const Decimal &x) {
    Decimal zero(0);
    return (x < zero) ? -x : x;
}


// Helper to parse "[-]ddd.ddd" (NO exponent)
inline Decimal decimal_from_simple_string(const std::string &s)
{
    if (s.empty()) return Decimal(0);

    int sign = 1;
    std::size_t pos = 0;
    if (s[pos] == '+') {
        ++pos;
    } else if (s[pos] == '-') {
        sign = -1;
        ++pos;
    }

    std::size_t dot = s.find('.', pos);
    std::string int_part, frac_part;
    if (dot == std::string::npos) {
        int_part = s.substr(pos);
    } else {
        int_part = s.substr(pos, dot - pos);
        frac_part = s.substr(dot + 1);
    }

    // Strip non-digits (safeguard)
    auto strip = [](std::string &t) {
        t.erase(std::remove_if(t.begin(), t.end(), 
            [](unsigned char c){ return !isdigit(c); }), t.end());
    };

    strip(int_part);
    strip(frac_part);
    
    if (int_part.empty()) int_part = "0";

    // FIX: Force Base 10 to prevent Octal detection on leading zeros
    mpz_class num(int_part, 10);
    mpz_class den = 1;

    if (!frac_part.empty()) {
        // FIX: Force Base 10 here too
        mpz_class frac(frac_part, 10);
        mpz_class scale = pow10_z(static_cast<long>(frac_part.size()));
        num = num * scale + frac;
        den = scale;
    }

    if (sign < 0) num = -num;

    return Decimal(num, den);
}

// Main parser: Handles "123", "12.34", "1.2e-5", "1E10"
inline Decimal decimal_from_string(const std::string &s_in)
{
    std::string s = trim(s_in);
    if (s.empty()) return Decimal(0);

    // Check for exponent
    std::size_t e_pos = s.find('e');
    if (e_pos == std::string::npos) e_pos = s.find('E');

    if (e_pos == std::string::npos) {
        return decimal_from_simple_string(s);
    }

    // Split base and exponent
    std::string base_str = s.substr(0, e_pos);
    std::string exp_str = s.substr(e_pos + 1);

    Decimal base = decimal_from_simple_string(base_str);
    
    if (exp_str.empty()) return base;

    long exp_val = 0;
    try {
        exp_val = std::stol(exp_str);
    } catch (...) {
        return base; 
    }

    if (exp_val == 0) return base;

    if (exp_val > 0) {
        mpz_class scale = pow10_z(exp_val);
        return base * Decimal(scale);
    } else {
        mpz_class scale = pow10_z(-exp_val);
        return base / Decimal(scale);
    }
}

// Truncate x to "digits" decimal places (toward -infinity)
inline Decimal truncate_decimal(const Decimal &x, int digits)
{
    mpz_class scale;
    if (digits >= 0) scale = pow10_z(digits);
    else return Decimal(0); 

    Decimal y = x * Decimal(scale);
    
    mpz_class n = y.get_num();
    mpz_class d = y.get_den();
    mpz_class q;

    if (n >= 0) q = n / d;
    else q = -((-n + d - 1) / d); 

    return Decimal(q, scale);
}

// Convert Decimal to string
inline std::string decimal_to_string(const Decimal &x, int digits)
{
    Decimal q = x;
    bool neg = (q < 0);
    if (neg) q = -q;

    mpz_class n = q.get_num();
    mpz_class d = q.get_den();

    mpz_class integer = n / d;
    mpz_class rem = n % d;

    std::string res = integer.get_str();

    if (digits > 0) {
        std::string frac;
        frac.reserve(digits);
        for (int i = 0; i < digits; ++i) {
            rem *= 10;
            mpz_class digit = rem / d;
            rem %= d;
            frac.push_back(static_cast<char>('0' + digit.get_si()));
        }
        while (!frac.empty() && frac.back() == '0') frac.pop_back();
        if (!frac.empty()) res += "." + frac;
    }

    if (neg) res = "-" + res;
    return res;
}

// Count significant decimal digits in string
inline int count_decimals_str(const std::string &s_in)
{
    std::string s = trim(s_in);
    std::size_t dot = s.find('.');
    if (dot == std::string::npos) return 0;
    
    std::size_t e_pos = s.find('e');
    if (e_pos == std::string::npos) e_pos = s.find('E');
    
    std::string frac;
    if (e_pos == std::string::npos) {
        frac = s.substr(dot + 1);
    } else {
        if (e_pos > dot) frac = s.substr(dot + 1, e_pos - dot - 1);
    }

    while (!frac.empty() && frac.back() == '0') {
        frac.pop_back();
    }
    return static_cast<int>(frac.size());
}

// Range generation
inline std::vector<Decimal> frange_decimal(const Decimal &start,
                                           const Decimal &end,
                                           const Decimal &step)
{
    std::vector<Decimal> v;
    if (step <= Decimal(0)) return v;
    Decimal x = start;
    int safety = 0;
    const int MAX_STEPS = 1000000; 

    while (x <= end && safety < MAX_STEPS) {
        v.push_back(x);
        x += step;
        ++safety;
    }
    return v;
}



#endif // DECIMAL_HPP