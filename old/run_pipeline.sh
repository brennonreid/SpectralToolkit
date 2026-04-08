#!/usr/bin/env bash
set -euo pipefail

# run_pipeline.sh
#
# Pipeline:
#   1) heatsonar.exe --mode range ... -> writes heatmap + seeds (defaults)
#   2) for each seed T: nsniperv6 --t T --depth ... --max_step ... -> capture stdout
#
# Policy:
# - Strict: missing executables or files hard-fail.
# - No fallbacks.
# - Sniper output: stdout only (captured to logs). No --out. No extra derived parsing.

# -----------------------------
# Defaults (override via CLI)
# -----------------------------
HEATSONAR="${HEATSONAR:-./heatsonar.exe}"
NSNIPER="${NSNIPER:-./nsniperv6.exe}"

OUTDIR="${OUTDIR:-sniper_runs}"

# Heatsonar defaults (range)
MODE="range"
INT_START="${INT_START:-10}"
INT_END="${INT_END:-100}"
T_STEP="${T_STEP:-0.001}"
DPS="${DPS:-80}"
PRETTY_DIGITS="${PRETTY_DIGITS:-25}"
MAX_STEP_HEAT="${MAX_STEP_HEAT:-200}"
TARGET_ENERGY="${TARGET_ENERGY:-}"     # empty disables
ORDER="${ORDER:-T}"
CSV_HITS="${CSV_HITS:-}"               # empty disables

# Heatsonar outputs (defaults match your heatsonar.cpp)
HEATMAP_PATH="${HEATMAP_PATH:-heatmap.csv}"
SEEDS_PATH="${SEEDS_PATH:-seeds.txt}"

# nsniperv6 defaults (as you stated)
SNIPER_DEPTH="${SNIPER_DEPTH:-25}"
SNIPER_MAX_STEP="${SNIPER_MAX_STEP:-100}"
SNIPER_EXTRA="${SNIPER_EXTRA:-}"       # appended verbatim (use carefully)

SKIP_HEATSONAR=0
QUIET_PROGRESS=0

# -----------------------------
# Help
# -----------------------------
usage() {
  cat << 'USAGE'
Usage: ./run_pipeline.sh [options]

Core options:
  --skip_heatsonar            Do not run heatsonar; just batch existing --seeds file
  --outdir <dir>              Output directory for sniper logs (default: sniper_runs)
  --quiet                     Reduce progress output to minimal

Heatsonar (range mode) options:
  --heatsonar <path>          heatsonar executable (default: ./heatsonar.exe)
  --int_start <int>           (default: 10)
  --int_end <int>             (default: 100)
  --T_step <num>              (default: 0.001)
  --dps <int>                 (default: 80)
  --pretty_digits <int>       (default: 25)
  --max_step_heat <num>       (default: 200)
  --target_energy <dbl>       If set, passed to heatsonar
  --order <T|energy>          (default: T)
  --csv <path>                If set, heatsonar writes hits/report CSV there
  --heatmap <path>            Full heatmap CSV (default: heatmap.csv)
  --seeds <path>              Seeds output file (default: seeds.txt)

nsniperv6 options:
  --nsniper <path>            nsniperv6 executable (default: ./nsniperv6.exe)
  --sniper_depth <int>        (default: 25)
  --sniper_max_step <num>     (default: 100)
  --sniper_extra "<flags>"    Extra flags appended verbatim to nsniperv6 call

Notes:
- This script captures nsniperv6 stdout only. No --out is used.
- Per-seed log:  <outdir>/run_<token>.stdout.txt
- Combined log:  <outdir>/all_runs.stdout.txt
USAGE
}

# -----------------------------
# Strict arg parse
# -----------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --help|-h) usage; exit 0;;

    --skip_heatsonar) SKIP_HEATSONAR=1; shift;;
    --outdir) OUTDIR="$2"; shift 2;;
    --quiet) QUIET_PROGRESS=1; shift;;

    --heatsonar) HEATSONAR="$2"; shift 2;;
    --int_start) INT_START="$2"; shift 2;;
    --int_end) INT_END="$2"; shift 2;;
    --T_step) T_STEP="$2"; shift 2;;
    --dps) DPS="$2"; shift 2;;
    --pretty_digits) PRETTY_DIGITS="$2"; shift 2;;
    --max_step_heat) MAX_STEP_HEAT="$2"; shift 2;;
    --target_energy) TARGET_ENERGY="$2"; shift 2;;
    --order) ORDER="$2"; shift 2;;
    --csv) CSV_HITS="$2"; shift 2;;
    --heatmap) HEATMAP_PATH="$2"; shift 2;;
    --seeds) SEEDS_PATH="$2"; shift 2;;

    --nsniper) NSNIPER="$2"; shift 2;;
    --sniper_depth) SNIPER_DEPTH="$2"; shift 2;;
    --sniper_max_step) SNIPER_MAX_STEP="$2"; shift 2;;
    --sniper_extra) SNIPER_EXTRA="$2"; shift 2;;

    *)
      echo "ERROR: unknown arg: $1" >&2
      echo "Use --help for usage." >&2
      exit 1
      ;;
  esac
done

# -----------------------------
# Preconditions
# -----------------------------
mkdir -p "$OUTDIR"

if [[ ! -x "$NSNIPER" ]]; then
  echo "ERROR: nsniper executable not found/executable at: $NSNIPER" >&2
  exit 1
fi

if [[ "$SKIP_HEATSONAR" -eq 0 ]]; then
  if [[ ! -x "$HEATSONAR" ]]; then
    echo "ERROR: heatsonar executable not found/executable at: $HEATSONAR" >&2
    exit 1
  fi
fi

# -----------------------------
# Run heatsonar (optional)
# -----------------------------
if [[ "$SKIP_HEATSONAR" -eq 0 ]]; then
  hs_args=( "--mode" "range"
            "--int_start" "$INT_START"
            "--int_end" "$INT_END"
            "--T_step" "$T_STEP"
            "--dps" "$DPS"
            "--pretty_digits" "$PRETTY_DIGITS"
            "--max_step" "$MAX_STEP_HEAT"
            "--order" "$ORDER"
            "--heatmap" "$HEATMAP_PATH"
            "--seeds" "$SEEDS_PATH" )

  if [[ -n "$TARGET_ENERGY" ]]; then
    hs_args+=( "--target_energy" "$TARGET_ENERGY" )
  fi

  if [[ -n "$CSV_HITS" ]]; then
    hs_args+=( "--csv" "$CSV_HITS" )
  fi

  if [[ "$QUIET_PROGRESS" -eq 0 ]]; then
    echo "[PIPE] heatsonar: $HEATSONAR ${hs_args[*]}"
  fi

  "$HEATSONAR" "${hs_args[@]}"
fi

# -----------------------------
# Seeds must exist now
# -----------------------------
if [[ ! -f "$SEEDS_PATH" ]]; then
  echo "ERROR: seeds file not found: $SEEDS_PATH" >&2
  exit 1
fi

# -----------------------------
# Batch seeds through nsniperv6
# -----------------------------
ALL_LOG="$OUTDIR/all_runs.stdout.txt"
: > "$ALL_LOG"

if [[ "$QUIET_PROGRESS" -eq 0 ]]; then
  echo "[PIPE] seeds:  $SEEDS_PATH"
  echo "[PIPE] nsniper: $NSNIPER --depth $SNIPER_DEPTH --max_step $SNIPER_MAX_STEP $SNIPER_EXTRA"
  echo "[PIPE] outdir: $OUTDIR"
fi

seed_count=0

while IFS= read -r T0; do
  [[ -z "$T0" ]] && continue
  [[ "$T0" =~ ^[[:space:]]*# ]] && continue

  seed_count=$((seed_count + 1))

  # File-safe token
  token="$(echo "$T0" | tr ':/\ ' '____')"
  run_log="$OUTDIR/run_${token}.stdout.txt"

  if [[ "$QUIET_PROGRESS" -eq 0 ]]; then
    echo "[SEED $seed_count] T0=$T0 -> $run_log"
  fi

  {
    echo "============================================================"
    echo "SEED_INDEX: $seed_count"
    echo "T0: $T0"
    echo "CMD: $NSNIPER --t \"$T0\" --depth $SNIPER_DEPTH --max_step $SNIPER_MAX_STEP $SNIPER_EXTRA"
    echo "------------------------------------------------------------"
    "$NSNIPER" --t "$T0" --depth "$SNIPER_DEPTH" --max_step "$SNIPER_MAX_STEP" $SNIPER_EXTRA
    echo
  } > "$run_log"

  # Append to combined log
  cat "$run_log" >> "$ALL_LOG"

done < "$SEEDS_PATH"

if [[ "$QUIET_PROGRESS" -eq 0 ]]; then
  echo "[DONE] seeds_processed=$seed_count"
  echo "[DONE] combined_log=$ALL_LOG"
fi
