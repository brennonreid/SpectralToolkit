#ifndef ARB_INTERFACE_HPP
#define ARB_INTERFACE_HPP

#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/arf.h>
#include "decimal.hpp"
#include <string>

// COMPATIBILITY FIX:
#ifndef ARF_STR_NOEXPONENT
#define ARF_STR_NOEXPONENT 0
#endif

struct ArbZetaContext {
    slong prec;
    acb_t s;
    acb_t z;
    arb_t abs_z;

    explicit ArbZetaContext(slong bits_prec) : prec(bits_prec) {
        acb_init(s);
        acb_init(z);
        arb_init(abs_z);
    }

    ~ArbZetaContext() {
        acb_clear(s);
        acb_clear(z);
        arb_clear(abs_z);
    }

    // Set arb_t from Decimal safely via String
    // This bypasses memory layout mismatches between GMP and FLINT
    void arb_set_decimal_safe(arb_t out, const Decimal &q) {
        // Convert to string with high precision (40% of bits is roughly enough digits)
        std::string s_val = decimal_to_string(q, (int)(prec * 0.4)); 
        arb_set_str(out, s_val.c_str(), prec);
    }

    // Compute |zeta(s)|
    void compute(const Decimal &T, const Decimal &sigma) {
        arb_t re, im;
        arb_init(re);
        arb_init(im);

        // USE SAFE STRING SETTERS
        arb_set_decimal_safe(re, sigma);
        arb_set_decimal_safe(im, T);

        arb_set(acb_realref(s), re);
        arb_set(acb_imagref(s), im);

        acb_zeta(z, s, prec);
        acb_abs(abs_z, z, prec);

        arb_clear(re);
        arb_clear(im);
    }

    // Return energy as high-precision string
    std::string get_energy_str(slong digits) {
        char *s = arb_get_str(abs_z, digits, ARF_STR_NOEXPONENT);
        std::string res(s);
        flint_free(s);
        return res;
    }
    
    // Compatibility helper
    long double energy_at_sigma(const Decimal &T, const Decimal &sigma) {
        compute(T, sigma);
        const arf_struct *mid = arb_midref(abs_z);
        return arf_get_d(mid, ARF_RND_NEAR);
    }
};

#endif // ARB_INTERFACE_HPP