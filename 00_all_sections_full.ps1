Param(
    [string]$Root      = (Get-Location).Path,
    [string]$ToolsDir  = "",
    [string]$InputsDir = "",
    [string]$CertsDir  = "",
    [int]$Dps          = 950,

    [int]$Grid         = 12000,
    [double]$CritLeft  = 0.30,
    [double]$CritRight = 2.80,
    [double]$Sigma     = 6.0,
    [double]$K0        = 0.25,

    [string]$T0        = "10000000000000000",
    [string]$T1        = "100000000000000000",
    [string]$TMin      = "10000000000000000",
    [string]$TMax      = "100000000000000000",
    [string]$TStep     = "1000000000000000",

    [double]$X0        = 40.0,
    [double]$Ap        = 1.0,
    [double]$Ag        = 1.0,

    [int]$Digits       = 950,

    [string]$SigmaMid  = "6.0",
    [string]$K0Mid     = "0.25",
    [double]$SigmaSpan = 0.5,
    [double]$K0Width   = 0.10,
    [int]$SigmaSteps   = 25,
    [int]$K0Steps      = 25,
    [double]$LhsTotal  = 1.0,
    [double]$InnerLeft = 0.3,
    [double]$InnerRight = 2.8,
    [string]$ConeName  = "cone_main",
    [int]$Jobs         = 0,
    [string]$Executor  = "thread",

    [string]$DeltaTarget = "1e-14",
    [int]$MeshInitial    = 512,
    [int]$MeshMax        = 262144,

    [switch]$UseVenv = $true,
    [string]$VenvDir = ""
)

$ErrorActionPreference = "Stop"

function Msg($m) {
    Write-Host "[00_all_sections_full] $m" -ForegroundColor Cyan
}

function Die($m) {
    Write-Host "[00_all_sections_full] FATAL: $m" -ForegroundColor Red
    exit 1
}

# Resolve base paths
if (-not $ToolsDir) {
    $ToolsDir = Join-Path $Root "tools"
}
if (-not $CertsDir) {
    $CertsDir = Join-Path $Root "PROOF_PACKET"
}
if (-not $InputsDir) {
    $InputsDir = Join-Path $Root ""
}
if (-not $VenvDir) {
    $VenvDir = Join-Path $Root "rh-env"
}

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path

Msg "Root=$Root"
Msg "ToolsDir=$ToolsDir"
Msg "InputsDir=$InputsDir"
Msg "CertsDir=$CertsDir"
Msg "Dps=$Dps"

# Helper to run each stage and abort on failure
function Invoke-Stage {
    Param(
        [string]$Name,
        [scriptblock]$Block
    )
    Msg "==== Running $Name ===="
    try {
        & $Block
        Msg "==== $Name completed ===="
    }
    catch {
        Die "$Name failed: $_"
    }
}

# 01_inputs.ps1
Invoke-Stage "01_inputs.ps1" {
    & (Join-Path $ScriptDir "01_inputs.ps1") `
        -Root $Root `
        -ToolsDir $ToolsDir `
        -InputsDir $InputsDir `
        -CertsDir $CertsDir `
        -Dps $Dps `
        -Grid $Grid `
        -CritLeft $CritLeft `
        -CritRight $CritRight `
        -Sigma $Sigma `
        -K0 $K0 `
        -UseVenv:$UseVenv `
        -VenvDir $VenvDir
}

# 02_tails.ps1
Invoke-Stage "02_tails.ps1" {
    & (Join-Path $ScriptDir "02_tails.ps1") `
        -Root $Root `
        -ToolsDir $ToolsDir `
        -InputsDir $InputsDir `
        -CertsDir $CertsDir `
        -Dps $Dps `
        -Sigma $Sigma `
        -K0 $K0 `
        -T0 $T0 `
        -UseVenv:$UseVenv `
        -VenvDir $VenvDir
}

# 03_prime_block_and_grid.ps1
Invoke-Stage "03_prime_block_and_grid.ps1" {
    & (Join-Path $ScriptDir "03_prime_block_and_grid.ps1") `
        -Root $Root `
        -ToolsDir $ToolsDir `
        -InputsDir $InputsDir `
        -CertsDir $CertsDir `
        -Dps $Dps `
        -Sigma $Sigma `
        -K0 $K0 `
        -Grid $Grid `
        -UseVenv:$UseVenv `
        -VenvDir $VenvDir
}

# 04_weil_psd.ps1
Invoke-Stage "04_weil_psd.ps1" {
    & (Join-Path $ScriptDir "04_weil_psd.ps1") `
        -Root $Root `
        -ToolsDir $ToolsDir `
        -InputsDir $InputsDir `
        -CertsDir $CertsDir `
        -Dps $Dps `
        -UseVenv:$UseVenv `
        -VenvDir $VenvDir
}

# 05_continuum_rollup.ps1
Invoke-Stage "05_continuum_rollup.ps1" {
    & (Join-Path $ScriptDir "05_continuum_rollup.ps1") `
        -Root $Root `
        -ToolsDir $ToolsDir `
        -CertsDir $CertsDir `
        -Dps $Dps `
        -UseVenv:$UseVenv `
        -VenvDir $VenvDir
}

# 06_lipschitz_and_density.ps1
Invoke-Stage "06_lipschitz_and_density.ps1" {
    & (Join-Path $ScriptDir "06_lipschitz_and_density.ps1") `
        -Root $Root `
        -ToolsDir $ToolsDir `
        -CertsDir $CertsDir `
        -Dps $Dps `
        -UseVenv:$UseVenv `
        -VenvDir $VenvDir
}

# 07_tails_bundle_and_explicit.ps1
Invoke-Stage "07_tails_bundle_and_explicit.ps1" {
    & (Join-Path $ScriptDir "07_tails_bundle_and_explicit.ps1") `
        -Root $Root `
        -ToolsDir $ToolsDir `
        -InputsDir $InputsDir `
        -CertsDir $CertsDir `
        -Dps $Dps `
        -Sigma $Sigma `
        -K0 $K0 `
        -T0 $T0 `
        -X0 $X0 `
        -Ap $Ap `
        -Ag $Ag `
        -UseVenv:$UseVenv `
        -VenvDir $VenvDir
}

# 08_core_integral.ps1
Invoke-Stage "08_core_integral.ps1" {
    & (Join-Path $ScriptDir "08_core_integral.ps1") `
        -Root $Root `
        -ToolsDir $ToolsDir `
        -InputsDir $InputsDir `
        -CertsDir $CertsDir `
        -Dps $Dps `
        -T0 $T0 `
        -UseVenv:$UseVenv `
        -VenvDir $VenvDir
}

# 09_uniform_rollup.ps1 (currently mostly STP test)
Invoke-Stage "09_uniform_rollup.ps1" {
    & (Join-Path $ScriptDir "09_uniform_rollup.ps1") `
        -Root $Root `
        -ToolsDir $ToolsDir `
        -InputsDir $InputsDir `
        -CertsDir $CertsDir `
        -Dps $Dps `
        -T0 $T0 `
        -TMin $TMin `
        -TMax $TMax `
        -TStep $TStep `
        -UseVenv:$UseVenv `
        -VenvDir $VenvDir
}

# 10_infinite_analysis.ps1
Invoke-Stage "10_infinite_analysis.ps1" {
    & (Join-Path $ScriptDir "10_infinite_analysis.ps1") `
        -Root $Root `
        -ToolsDir $ToolsDir `
        -CertsDir $CertsDir `
        -Dps $Dps `
        -T0 $T0 `
        -X0 $X0 `
        -UseVenv:$UseVenv `
        -VenvDir $VenvDir
}

# 11_cone_uniform.ps1
Invoke-Stage "11_cone_uniform.ps1" {
    & (Join-Path $ScriptDir "11_cone_uniform.ps1") `
        -Root $Root `
        -ToolsDir $ToolsDir `
        -CertsDir $CertsDir `
        -Grid $Grid `
        -Digits $Digits `
        -Dps $Dps `
        -SigmaMid ("{0}" -f $Sigma) `
        -K0Mid ("{0}" -f $K0) `
        -SigmaSpan ("{0}" -f $SigmaSpan) `
        -K0Width ("{0}" -f $K0Width) `
        -SigmaSteps $SigmaSteps `
        -K0Steps $K0Steps `
        -CriticalLeft ("{0}" -f $CritLeft) `
        -CriticalRight ("{0}" -f $CritRight) `
        -InnerLeft ("{0}" -f $InnerLeft) `
        -InnerRight ("{0}" -f $InnerRight) `
        -ConeName $ConeName `
        -Jobs $Jobs `
        -UseVenv:$UseVenv `
        -VenvDir $VenvDir `
        -PyExe $PyExe
}

# 12_rolling_t.ps1
Invoke-Stage "12_rolling_t.ps1" {
    & (Join-Path $ScriptDir "12_rolling_t.ps1") `
        -Root $Root `
        -ToolsDir $ToolsDir `
        -CertsDir $CertsDir `
        -T0 $T0 `
        -T1 $T1 `
        -DeltaTarget $DeltaTarget `
        -MeshInitial $MeshInitial `
        -MeshMax $MeshMax `
        -Dps $Dps `
        -UseVenv:$UseVenv `
        -VenvDir $VenvDir
}

# 13_report_wrap.ps1
Invoke-Stage "13_report_wrap.ps1" {
    & (Join-Path $ScriptDir "13_report_wrap.ps1") `
        -Root $Root `
        -ToolsDir $ToolsDir `
        -CertsDir $CertsDir `
        -InputsDir $InputsDir `
        -Dps $Dps `
        -UseVenv:$UseVenv `
        -VenvDir $VenvDir
}

Msg "All stages complete."
