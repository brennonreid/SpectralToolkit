[CmdletBinding()]
Param(
  # Paths
  [string]$Root     = $(Get-Location).Path,
  [string]$ToolsDir = "",
  [string]$CertsDir = "",

  # Numerics / grid
  [int]$Grid   = 8000,
  [int]$Digits = 220,    # decimal digits for printed output (if used by tool)
  [int]$Dps    = 950,    # working precision

  # Cone center
  [string]$SigmaMid = "6.0",
  [string]$K0Mid    = "0.25",

  # Cone extents
  [string]$SigmaSpan  = "0.5",  # total span around SigmaMid (so [mid-span/2, mid+span/2])
  [string]$K0Width    = "0.1",  # total width around K0Mid (so [mid-width/2, mid+width/2])

  # Discrete grid in sigma / k0
  [int]$SigmaSteps = 15,
  [int]$K0Steps    = 15,

  # Cone label
  [string]$ConeName = "cone_main",

  # Critical / inner intervals (for tool CLI)
  [string]$CriticalLeft  = "0.3",
  [string]$CriticalRight = "2.8",
  [string]$InnerLeft     = "0.3",
  [string]$InnerRight    = "2.8",

  # CSV subdir name (under CertsDir)
  [string]$CsvSubDirName = "param_cone_uniform_csv",

  # Python / venv
  [switch]$UseVenv = $true,
  [string]$VenvDir = "",

  # Explicit python override (optional)
  [string]$PyExe = "",

  # Parallelism hint (passed through to tool)
  [int]$Jobs = 0,
  [string]$Executor = "thread"
)

# ---------------- helpers ----------------

$Tag = "11_cone_uniform"

function Msg {
  param([string]$Message)
  Write-Host "[$Tag] $Message"
}

function Warn {
  param([string]$Message)
  Write-Warning "[$Tag] $Message"
}

function Die {
  param([string]$Message)
  Write-Error "[$Tag] $Message"
  exit 1
}

function Assert-File {
  param(
    [string]$Path,
    [string]$Err
  )
  if (-not (Test-Path $Path)) {
    Die $Err
  }
}

function New-Dir {
  param([string]$Path)
  if (-not (Test-Path $Path)) {
    New-Item -ItemType Directory -Path $Path | Out-Null
  }
}

function Resolve-Python {
  param(
    [string]$Root,
    [switch]$UseVenv,
    [string]$VenvDir,
    [string]$PyExe
  )

  if ($PyExe) {
    return $PyExe
  }

  if ($UseVenv) {
    if (-not $VenvDir) {
      $VenvDir = Join-Path $Root "rh-env"
    }
    $venvPy = Join-Path $VenvDir "Scripts\python.exe"
    if (Test-Path $venvPy) {
      return $venvPy
    } else {
      Warn "venv requested but $venvPy not found, using python on PATH"
    }
  }

  return "python"
}

# ---------------- defaults / paths ----------------

if (-not $ToolsDir) {
  $ToolsDir = Join-Path $Root "tools"
}
if (-not $CertsDir) {
  $CertsDir = Join-Path $Root "PROOF_PACKET"
}

Msg "Cone uniform cert stage"
Msg "Root=$Root, ToolsDir=$ToolsDir, CertsDir=$CertsDir, Dps=$Dps"

# Basic sigma / k0 range diagnostics
$SigmaMidNum  = [double]$SigmaMid
$SigmaSpanNum = [double]$SigmaSpan
$SigmaLo      = $SigmaMidNum - $SigmaSpanNum / 2.0
$SigmaHi      = $SigmaMidNum + $SigmaSpanNum / 2.0

$K0MidNum   = [double]$K0Mid
$K0WidthNum = [double]$K0Width
$K0Lo       = $K0MidNum - $K0WidthNum / 2.0
$K0Hi       = $K0MidNum + $K0WidthNum / 2.0

Msg ("Sigma in [{0}, {1}], K0 in [{2}, {3}], ConeName={4}" -f $SigmaLo, $SigmaHi, $K0Lo, $K0Hi, $ConeName)

# Upstream certs
$UniformCert = Join-Path $CertsDir "uniform_certificate.json"
$Density     = Join-Path $CertsDir "density_prover.json"

# ---------------- LhsTotal from uniform_certificate.json ----------------

# Default
$LhsTotal = "1"

if (Test-Path $UniformCert) {
  Msg "uniform_certificate.json found; extracting lhs_total"
  try {
    $UniformObj = Get-Content $UniformCert -Raw | ConvertFrom-Json

    $lhs = $UniformObj.uniform_certificate.lhs_total
    if ($lhs -and ($lhs -ne "")) {
      $LhsTotal = "$lhs"
      Msg "Using lhs_total=$LhsTotal from uniform_certificate.uniform_certificate.lhs_total"
    } else {
      Warn "uniform_certificate.lhs_total missing or empty; falling back to LhsTotal=1"
      $LhsTotal = "1"
    }

    if ($null -ne $UniformObj.PASS -and (-not [bool]$UniformObj.PASS)) {
      Warn "uniform_certificate.PASS is false; cone uniform cert will still run using lhs_total, but upstream uniform certificate did not pass."
    }
  } catch {
    Warn ("Failed to parse lhs_total from uniform_certificate.json; falling back to LhsTotal=1 ({0})" -f $_.Exception.Message)
    $LhsTotal = "1"
  }
} else {
  Warn "uniform_certificate.json not found (run 08_core_integral.ps1 if you want a matched lhs_total); proceeding with LhsTotal=1"
}

# Density is optional; we just warn if missing (kept for compatibility / diagnostics)
if (-not (Test-Path $Density)) {
  Warn "density_prover.json not found; continuing without density reference"
}

# ---------------- Python tool setup ----------------

$Tool = Join-Path $ToolsDir "param_cone_uniform_cert_v4c.py"
Assert-File $Tool "Missing param_cone_uniform_cert_v4c.py in $ToolsDir"

$Py = Resolve-Python -Root $Root -UseVenv:$UseVenv -VenvDir $VenvDir -PyExe $PyExe

# CSV dir and output paths
$CsvDir  = Join-Path $CertsDir ("{0}_{1}" -f $CsvSubDirName, $ConeName)
New-Dir $CsvDir

$OutJson = Join-Path $CertsDir ("param_cone_uniform_cert_{0}.json" -f $ConeName)

# ---------------- Build CLI args ----------------

$args = @(
  "--sigma-mid",  "$SigmaMid",
  "--sigma-span", "$SigmaSpan",
  "--sigma-steps", $SigmaSteps,
  "--k0-mid",     "$K0Mid",
  "--k0-width",   "$K0Width",
  "--k0-steps",   $K0Steps,
  "--lhs-total",  "$LhsTotal",
  "--grid",       $Grid,
  "--digits",     $Digits,
  "--critical-left",  "$CriticalLeft",
  "--critical-right", "$CriticalRight",
  "--inner-left",     "$InnerLeft",
  "--inner-right",    "$InnerRight",
  "--csv-dir",    $CsvDir,
  "--dps",        $Dps,
  "--out",        $OutJson
)

if ($Jobs -gt 0) {
  $args += @("--jobs", ("{0}" -f $Jobs))
}
if ($Executor) {
  $args += @("--executor", $Executor)
}

Msg ("python {0} {1}" -f $Tool, ($args -join " "))
& $Py $Tool $args

Assert-File $OutJson "param_cone_uniform_cert_v4c did not produce $OutJson"

Msg "Cone uniform cert stage complete."
