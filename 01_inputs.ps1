Param(
  [string]$Root = $(Get-Location).Path,
  [string]$ToolsDir = "",
  [string]$InputsDir = "",
  [string]$CertsDir = "",
  [int]$Dps = 220,
  [int]$Grid = 6000,
  [double]$CritLeft = 0.30,
  [double]$CritRight = 2.80,
  [double]$Sigma = 6.0,
  [double]$K0 = 0.025,
  [switch]$UseVenv = $true,
  [string]$VenvDir = ""
)

function Msg($m) { Write-Host "[01_inputs] $m" -ForegroundColor Cyan }
function Die($m) { throw $m }
function Assert-File($p, $m) { if (-not (Test-Path $p)) { Die $m } }

if (-not $ToolsDir)  { $ToolsDir  = Join-Path $Root "tools" }
if (-not $InputsDir) { $InputsDir = Join-Path $Root "packs\rh\inputs" }
if (-not $CertsDir)  { $CertsDir  = Join-Path $Root "PROOF_PACKET" }
if (-not $VenvDir)   { $VenvDir   = Join-Path $Root "rh-env" }

$Py = "python"
if ($UseVenv) {
  $PyCandidate = Join-Path $VenvDir "Scripts\python.exe"
  if (Test-Path $PyCandidate) {
    $Py = $PyCandidate
  } else {
    Msg "Warning: venv requested but $PyCandidate not found, falling back to python on PATH."
  }
}

Msg "Inputs stage (window + bands + band_cert)"
Msg "Root=$Root, ToolsDir=$ToolsDir, CertsDir=$CertsDir, Dps=$Dps, Grid=$Grid"

# 1) window_gen -> window.json
$WindowJson = Join-Path $CertsDir "window.json"
$wg = Join-Path $ToolsDir "window_gen.py"
if (Test-Path $wg) {
  Msg "python window_gen.py --mode gauss --sigma $Sigma --k0 $K0 --out $WindowJson --dps $Dps"
  & $Py $wg --mode "gauss" --sigma "$Sigma" --k0 "$K0" --out $WindowJson --dps "$Dps"
}
Assert-File $WindowJson "window_gen must produce window.json or provide it manually"

# 2) bands_make -> bands_auto.json
$AutoBands = Join-Path $CertsDir "bands_auto.json"
$bm = Join-Path $ToolsDir "bands_make.py"
if (Test-Path $bm) {
  Msg "python bands_make.py --window-config $WindowJson --critical-left $CritLeft --critical-right $CritRight --grid $Grid --out $AutoBands --dps $Dps"
  & $Py $bm `
    --window-config $WindowJson `
    --critical-left "$CritLeft" `
    --critical-right "$CritRight" `
    --grid "$Grid" `
    --out $AutoBands `
    --dps "$Dps"
}
Assert-File $AutoBands "bands_make must produce bands_auto.json"

# 3) band_cert -> band_cert.json
$BandCert = Join-Path $CertsDir "band_cert.json"
$bc = Join-Path $ToolsDir "band_cert.py"
if (Test-Path $bc) {
  Msg "python band_cert.py --window-config $WindowJson --bands $AutoBands --out $BandCert --dps $Dps --tqdm"
  & $Py $bc `
    --window-config $WindowJson `
    --bands $AutoBands `
    --out $BandCert `
    --dps "$Dps" `
    --tqdm
}
Assert-File $BandCert "band_cert must produce band_cert.json"

Msg "Inputs stage complete."
