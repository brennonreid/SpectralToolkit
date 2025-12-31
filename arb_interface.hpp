// arb_interface.hpp
#ifndef ARB_INTERFACE_HPP
#define ARB_INTERFACE_HPP

#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/arf.h>
#include <flint/fmpz.h>
#include <flint/acb_dirichlet.h>   // acb_dirichlet_hardy_z
#include <flint/dirichlet.h>

#include "decimal.hpp"

#include <string>
#include <iostream>
#include <cstdlib>
#include <gmpxx.h>

// COMPATIBILITY FIX:
#ifndef ARF_STR_NOEXPONENT
#define ARF_STR_NOEXPONENT 0
#endif

struct ArbZetaContext {
    slong prec;
    acb_t s;
    acb_t z;
    arb_t abs_z;

    // Hardy Z jet storage: res[0] = Z(t), res[1] = Z'(t)
    acb_struct hardy_jet[2];

    arb_t hardy_z;   // real part of Z(t)
    arb_t hardy_zp;  // real part of Z'(t)

    // Reuse temporaries (avoid per-call init/clear at huge prec)
    arb_t re;
    arb_t im;

    explicit ArbZetaContext(slong bits_prec) : prec(bits_prec) {
        acb_init(s);
        acb_init(z);
        arb_init(abs_z);

        acb_init(hardy_jet + 0);
        acb_init(hardy_jet + 1);


        arb_init(hardy_z);
        arb_init(hardy_zp);

        arb_init(re);
        arb_init(im);
    }

    ~ArbZetaContext() {
        acb_clear(s);
        acb_clear(z);
        arb_clear(abs_z);

        acb_clear(hardy_jet + 0);
        acb_clear(hardy_jet + 1);

        arb_clear(hardy_z);
        arb_clear(hardy_zp);

        arb_clear(re);
        arb_clear(im);
    }

    // Optional: allow changing default precision for subsequent calls.
    void set_prec(slong bits_prec) { prec = bits_prec; }

    // ------------------------------------------------------------------
    // EXACT Decimal -> arb_t (NO STRINGS)
    // out = q.num / q.den exactly, then rounded to Arb precision.
    // ------------------------------------------------------------------
    void arb_set_decimal_exact(arb_t out, const Decimal &q) {
        arb_set_decimal_exact_prec(out, q, prec);
    }

    // NEW (additive): caller-chosen precision
    void arb_set_decimal_exact_prec(arb_t out, const Decimal &q, slong prec_bits) {
        mpz_class n = q.get_num();
        mpz_class d = q.get_den();

        if (d == 0) {
            std::cerr << "ERROR: Decimal denominator is 0\n";
            std::exit(1);
        }

        fmpz_t fn, fd;
        fmpz_init(fn);
        fmpz_init(fd);

        fmpz_set_mpz(fn, n.get_mpz_t());
        fmpz_set_mpz(fd, d.get_mpz_t());

        arb_set_fmpz(out, fn);
        arb_div_fmpz(out, out, fd, prec_bits);

        fmpz_clear(fn);
        fmpz_clear(fd);
    }

    // ------------------------------------------------------------------
    // BACKWARD COMPAT: Compute |zeta(sigma + iT)|
    // ------------------------------------------------------------------
    void compute(const Decimal &T, const Decimal &sigma) {
        compute_prec(T, sigma, prec);
    }

    // Allow compute() to accept explicit precision
    void compute(const Decimal &T, const Decimal &sigma, slong step_prec) {
        compute_prec(T, sigma, step_prec);
    }

    void compute_prec(const Decimal &T, const Decimal &sigma, slong prec_bits) {
        arb_set_decimal_exact_prec(re, sigma, prec_bits);
        arb_set_decimal_exact_prec(im, T, prec_bits);

        arb_set(acb_realref(s), re);
        arb_set(acb_imagref(s), im);

        acb_zeta(z, s, prec_bits);
        acb_abs(abs_z, z, prec_bits);
    }

    // ------------------------------------------------------------------
    // BACKWARD COMPAT: Compute |HardyZ(T)| (critical line)
    // ------------------------------------------------------------------
    void compute_hardy_abs(const Decimal &T) {
        compute_hardy_abs_prec(T, prec);
    }

    void compute_hardy_abs_prec(const Decimal &T, slong prec_bits) {
        // Hardy Z takes a real argument; pass as acb with imag=0.
        arb_set_decimal_exact_prec(re, T, prec_bits);
        arb_zero(im);

        arb_set(acb_realref(s), re);
        arb_set(acb_imagref(s), im);

        acb_dirichlet_hardy_z(z, s, NULL, NULL, 1, prec_bits);
        acb_abs(abs_z, z, prec_bits);
    }

    // ------------------------------------------------------------------
    // NEW: Compute Hardy Z(T) AND Z'(T) in a single call (len=2).
    // - hardy_z  = Re(Z(T))
    // - hardy_zp = Re(Z'(T))
    // - abs_z    = |hardy_z| (for backwards-compatible magnitude logic)
    // ------------------------------------------------------------------
    void compute_hardy_z_and_deriv(const Decimal &T) {
        compute_hardy_z_and_deriv_prec(T, prec);
    }

    void compute_hardy_z_and_deriv_prec(const Decimal &T, slong prec_bits) {
        // Hardy Z takes a real argument; pass as acb with imag=0.
        arb_set_decimal_exact_prec(re, T, prec_bits);
        arb_zero(im);

        arb_set(acb_realref(s), re);
        arb_set(acb_imagref(s), im);

        // len=2: jet[0] = Z(t), jet[1] = Z'(t)
        acb_dirichlet_hardy_z(hardy_jet, s, NULL, NULL, 2, prec_bits);

        // Extract the real parts (Hardy Z should be real on real input).
        arb_set(hardy_z,  acb_realref(hardy_jet + 0));
        arb_set(hardy_zp, acb_realref(hardy_jet + 1));


        // Keep abs_z as |hardy_z| for existing magnitude-based code paths.
        arb_abs(abs_z, hardy_z);
    }

    // DISPLAY / logging only
    std::string get_energy_str(slong digits) {
        char *ss = arb_get_str(abs_z, digits, ARF_STR_NOEXPONENT);
        std::string res(ss);
        flint_free(ss);
        return res;
    }

    long double energy_at_sigma(const Decimal &T, const Decimal &sigma) {
        compute(T, sigma);
        const arf_struct *mid = arb_midref(abs_z);
        return arf_get_d(mid, ARF_RND_NEAR);
    }
};

#endif // ARB_INTERFACE_HPP
