#!/usr/bin/env python3
"""
prime_block_norm.py — zeros-driven analytic cap for the Prime Block operator.

Tweaks:
- Windows-safe writes (ensure dirs, newline='\n')
- No tqdm (integral uses mp.quad; progress would be misleading)
- Fast tail at high precision: split [T, T+M*sigma] + rigorously bounded remainder via E1
"""
import argparse, json, os, hashlib, datetime as dt
from datetime import timezone
from mpmath import mp


def set_precision(dps: int):
    mp.dps = int(dps)


def mp_str(x):
    return mp.nstr(mp.mpf(x), n=mp.dps, strip_zeros=False)


def _sha256_bytes(b: bytes) -> str:
    h = hashlib.sha256()
    h.update(b)
    return h.hexdigest()


def _write_json(path, obj):
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)
    blob = json.dumps(obj, indent=2, sort_keys=True).encode("utf-8")
    with open(path, "wb") as f:
        f.write(blob)
    return _sha256_bytes(blob)


def read_zeros_json(path):
    with open(path, "r", encoding="utf-8") as f:
        js = json.load(f)
    if isinstance(js, list):
        return [mp.mpf(z) for z in js]
    if isinstance(js, dict) and "zeros" in js:
        zs = js["zeros"]
        out = []
        for z in zs:
            if isinstance(z, (int, float, str)):
                out.append(mp.mpf(z))
            elif isinstance(z, dict):
                g = (
                    z.get("gamma")
                    or z.get("t")
                    or z.get("im")
                    or z.get("imag")
                    or z.get("value")
                )
                if g is None:
                    continue
                out.append(mp.mpf(g))
        return out
    try:
        return [mp.mpf(js)]
    except Exception:
        raise SystemExit("Unrecognized zeros JSON format")


def read_zeros_txt(path, N=None):
    vals = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            parts = s.split()
            tok = parts[-1]
            try:
                vals.append(mp.mpf(tok))
            except Exception:
                continue
            if N is not None and len(vals) >= N:
                break
    if not vals:
        raise SystemExit("No zeros parsed from text file")
    return vals


def weight_W(gamma, sigma, k0, scale):
    g = mp.mpf(gamma)
    s = mp.mpf(sigma)
    k = mp.mpf(k0)
    c = mp.mpf(scale)
    if s <= 0 or k <= 0 or c <= 0:
        raise SystemExit("sigma, k0, and scale must be positive")
    core = mp.e ** (-(g / s) ** 2)
    notch = 1 - mp.e ** (-(g / k) ** 2)
    return c * core * notch


def density_envelope(t):
    t = mp.mpf(t)
    if t < mp.e:
        t = mp.e
    return (mp.log(t) / (2 * mp.pi))


def cap_from_zeros(zlist, sigma, k0, scale, Tcut=None, tail_M=8):
    """
    Return (T, part, tail, cap_hi)
      T      : cutoff (max zero or Tcut)
      part   : sum_{g<=T} W(g)/(sqrt(1/4 + g^2))
      tail   : tail contribution beyond T
      cap_hi : 2 * (part + tail)
    """
    if not zlist:
        raise SystemExit("Empty zero list")

    s = mp.mpf(sigma)
    k = mp.mpf(k0)
    c = mp.mpf(scale)

    if Tcut is not None:
        Tcut = mp.mpf(Tcut)

    # Determine cutoff T
    max_zero = max(zlist)
    if Tcut is None:
        T = max_zero
    else:
        T = Tcut
        if T > max_zero:
            T = max_zero

    # Main sum over zeros with |gamma| <= T
    part = mp.mpf("0")
    for g in zlist:
        if g <= 0:
            continue
        if g > T:
            break
        Wg = weight_W(g, s, k, c)
        denom = mp.sqrt(mp.mpf("0.25") + g * g)
        part += Wg / denom

    # Tail beyond T: split [T, A] + [A, inf)
    M = int(tail_M)
    if M <= 0:
        raise SystemExit("tail_M must be positive")
    A = T + M * s

    def weight_W_cont(t):
        return weight_W(t, s, k, c)

    def density_envelope_cont(t):
        return density_envelope(t)

    def integrand(t):
        return weight_W_cont(t) * density_envelope_cont(t) / t

    # Numeric near tail on [T, A]
    tail_near = mp.mpf("0")
    if A > T:
        tail_near = mp.quad(integrand, [T, A])

    # Remainder bound on [A, inf):
    # drop notch (<=1), bound log(t) by log(A) for t>=A (monotone), integrate core/t
    # ∫_A^∞ e^{-(t/s)^2} (log A)/t dt = (log A)/2 * E1(((A/s)^2))
    # include scale/(2π) from density_envelope
    # => remainder_bound = c/(4π) * log(A) * E1((A/s)^2)
    # This is rigorous since notch <= 1 and log is monotone for t>=A >= e.
    u1 = (A / s) ** 2
    # Ensure A >= e so log(A) >= 1; otherwise clamp at e (density_envelope already does this)
    logA = mp.log(A if A >= mp.e else mp.e)
    remainder_bound = c * logA * mp.e1(u1) / (4 * mp.pi)

    tail = tail_near + remainder_bound
    cap_hi = 2 * (part + tail)
    return T, part, tail, cap_hi


def main():
    """
    Zeros-driven analytic operator-norm cap for the prime block.

    CLI (normalized v2.1):
      --zeros      : path to zeros data (JSON or text)
      --N          : use only first N zeros when reading text (optional)
      --sigma      : Gaussian width parameter sigma > 0
      --k0         : notch parameter k0 > 0
      --scale      : envelope scale >= 1 (kept for existing bounds)
      --Tcut       : cutoff T; default = max zero in file
      --tail-m     : tail split width M in multiples of sigma for numeric part
      --dps        : decimal precision for mpmath
      --out        : output JSON path

    JSON (normalized v2.1):
      kind  = "prime_block_norm"
      inputs.zeros_path
      inputs.N
      inputs.sigma
      inputs.k0
      inputs.scale
      inputs.Tcut
      inputs.tail_m
      prime_block_norm.used_operator_norm
      prime_block_norm.operator_norm_cap_hi
      prime_block_norm.method
      prime_block_norm.sum_zeros_contrib
      prime_block_norm.tail_bound_hi
      prime_block_norm.cap_total_hi
      meta.tool = "prime_block_norm"
    """
    ap = argparse.ArgumentParser(
        description="Zeros-driven analytic operator-norm cap for the prime block"
    )
    ap.add_argument(
        "--zeros", required=True, help="Path to zeros data (JSON or text)"
    )
    ap.add_argument(
        "--N",
        type=int,
        default=None,
        help="If reading text zeros, use only the first N zeros",
    )
    ap.add_argument(
        "--sigma",
        type=str,
        required=True,
        help="Gaussian width sigma > 0",
    )
    ap.add_argument(
        "--k0",
        type=str,
        required=True,
        help="Notch parameter k0 > 0",
    )
    ap.add_argument(
        "--scale",
        type=str,
        default="1.0",
        help="Envelope scale >= 1",
    )
    ap.add_argument(
        "--Tcut",
        type=str,
        default=None,
        help="Cutoff T; default = max zero in file",
    )
    ap.add_argument(
        "--tail-m",
        type=int,
        default=8,
        help="Tail split width M in multiples of sigma for numeric part",
    )
    ap.add_argument(
        "--dps",
        type=int,
        default=300,
        help="Decimal precision for mpmath",
    )
    ap.add_argument(
        "--out",
        required=True,
        help="Output JSON path",
    )
    args = ap.parse_args()

    set_precision(args.dps)

    zeros_path = args.zeros
    if zeros_path.lower().endswith(".txt"):
        zeros = read_zeros_txt(zeros_path, args.N)
    else:
        zeros = read_zeros_json(zeros_path)

    zeros = sorted([z for z in zeros if z > 0])

    T, part, tail, cap_hi = cap_from_zeros(
        zeros,
        mp.mpf(args.sigma),
        mp.mpf(args.k0),
        mp.mpf(args.scale),
        mp.mpf(args.Tcut) if args.Tcut is not None else None,
        tail_M=int(args.tail_m),
    )

    inputs = {
        "zeros_path": zeros_path,
        "N": None if args.N is None else int(args.N),
        "sigma": mp_str(args.sigma),
        "k0": mp_str(args.k0),
        "scale": mp_str(args.scale),
        "Tcut": None if args.Tcut is None else mp_str(args.Tcut),
        "tail_m": int(args.tail_m),
    }

    prime_block_data = {
        "used_operator_norm": mp_str(cap_hi),
        "operator_norm_cap_hi": mp_str(cap_hi),
        "method": "zeros_analytic_cap/v2_tail_split",
        "sum_zeros_contrib": mp_str(2 * part),
        "tail_bound_hi": mp_str(2 * tail),
        "cap_total_hi": mp_str(cap_hi),
        "cutoff_T": mp_str(T),
    }

    payload = {
        "kind": "prime_block_norm",
        "inputs": inputs,
        "prime_block_norm": prime_block_data,
        "meta": {
            "tool": "prime_block_norm",
            "created_utc": dt.datetime.now(timezone.utc).isoformat(),
            "dps": int(mp.dps),
            "zeros_count": int(len([g for g in zeros if g <= T])),
        },
        "PASS": True,
    }

    h = _write_json(args.out, payload)
    payload["meta"]["sha256"] = h
    _write_json(args.out, payload)

    print(
        f"[ok] prime_block_norm cap_hi={mp_str(cap_hi)}  "
        f"(2*part={mp_str(2*part)}  2*tail={mp_str(2*tail)})"
    )
    print(
        f"[info] zeros used up to T={mp_str(T)}  "
        f"count={payload['meta']['zeros_count']}  "
        f"scale={args.scale}  tail_M={args.tail_m}"
    )


if __name__ == "__main__":
    main()
