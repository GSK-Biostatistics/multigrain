<#
.SYNOPSIS
  Build every diagram spec in specs/ to out/<name>.excalidraw and out/<name>.svg.

.DESCRIPTION
  Stage 1 (build_diagrams.R)  spec YAML  -> Excalidraw scene, with validation
  Stage 2 (export.mjs)        scene      -> SVG, via Playwright + Chromium

  Both stages are offline: stage 2 serves @excalidraw/excalidraw and its fonts
  from node_modules over a throwaway localhost port and fails if the page tries
  to reach anything else.

.PARAMETER Spec
  Build only this spec (name, with or without .yaml). Default: all of specs/.

.PARAMETER Png
  Also rasterise each SVG to out/<name>.png, for a quick visual check.

.EXAMPLE
  .\build.ps1
  .\build.ps1 -Spec gsd-architecture -Png
#>
[CmdletBinding()]
param(
  [string] $Spec,
  [switch] $Png
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$here = $PSScriptRoot
$specDir = Join-Path $here "specs"
$outDir = Join-Path $here "out"

function Assert-Tool($name, $probe) {
  if (-not (Get-Command $name -ErrorAction SilentlyContinue)) {
    throw "$name is not on PATH; $probe"
  }
}
Assert-Tool "Rscript" "install R, or add R's bin\x64 directory to PATH"
Assert-Tool "node" "install Node.js"

if (-not (Test-Path (Join-Path $here "node_modules\@excalidraw\excalidraw"))) {
  throw "node_modules is incomplete; run 'npm install' in $here"
}

# @(...) around the whole assignment: the output of an if-block with a single
# item is otherwise unrolled to a scalar, which has no .Count under StrictMode.
$specs = @(if ($Spec) {
  $name = if ($Spec.EndsWith(".yaml") -or $Spec.EndsWith(".yml")) { $Spec } else { "$Spec.yaml" }
  $p = Join-Path $specDir $name
  if (-not (Test-Path $p)) { throw "no such spec: $p" }
  Get-Item $p
} else {
  Get-ChildItem -Path $specDir -Filter "*.yaml" -File
})
if ($specs.Count -eq 0) { throw "no specs found in $specDir" }

New-Item -ItemType Directory -Force -Path $outDir | Out-Null

# Stage 1: every spec in one R process, so the R startup cost is paid once.
& Rscript (Join-Path $here "build_diagrams.R") @($specs | ForEach-Object { $_.FullName })
if ($LASTEXITCODE -ne 0) { throw "stage 1 failed (exit $LASTEXITCODE)" }

# Stage 2: one browser launch per scene.
Write-Host "stage 2: excalidraw -> svg"
foreach ($s in $specs) {
  $name = [IO.Path]::GetFileNameWithoutExtension($s.Name)
  $scene = Join-Path $outDir "$name.excalidraw"
  $svg = Join-Path $outDir "$name.svg"
  if (-not (Test-Path $scene)) { throw "stage 1 produced no scene for $name" }
  $argv = @((Join-Path $here "export.mjs"), $scene, $svg)
  if ($Png) { $argv += @("--png", (Join-Path $outDir "$name.png")) }
  & node $argv
  if ($LASTEXITCODE -ne 0) { throw "stage 2 failed for $name (exit $LASTEXITCODE)" }
}

Write-Host "done -> $outDir"
