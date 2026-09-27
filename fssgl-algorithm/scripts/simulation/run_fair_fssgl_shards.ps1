param([int]$ShardCount = 5)

$ErrorActionPreference = 'Stop'
if ($ShardCount -lt 2) { throw 'ShardCount must be at least two.' }
$projectRoot = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
$rscriptCommand = Get-Command Rscript -ErrorAction SilentlyContinue
if ($rscriptCommand) {
  $rscriptExe = $rscriptCommand.Source
} else {
  $candidates = Get-ChildItem 'C:\Program Files\R\R-*\bin\Rscript.exe' -ErrorAction SilentlyContinue |
    Sort-Object FullName -Descending
  if (-not $candidates) { throw 'Rscript was not found.' }
  $rscriptExe = $candidates[0].FullName
}
$logDir = Join-Path $projectRoot 'tmp\fssgl_shards'
New-Item -ItemType Directory -Force -Path $logDir | Out-Null
foreach ($variableName in @('LANG', 'LC_ALL', 'LC_COLLATE', 'LC_CTYPE', 'LC_MONETARY', 'LC_TIME')) {
  if ([Environment]::GetEnvironmentVariable($variableName) -eq 'C.UTF-8') {
    Remove-Item "Env:$variableName"
  }
}
$env:FSSGL_FAIR_SHARD_COUNT = [string]$ShardCount
$processes = @()
for ($id = 1; $id -le $ShardCount; $id++) {
  $env:FSSGL_FAIR_SHARD_ID = [string]$id
  $stdout = Join-Path $logDir "shard$id.stdout.log"
  $stderr = Join-Path $logDir "shard$id.stderr.log"
  $processes += Start-Process -FilePath $rscriptExe `
    -ArgumentList 'scripts/simulation/11_fair_fssgl_comparison.R' `
    -WorkingDirectory $projectRoot -WindowStyle Hidden -PassThru `
    -RedirectStandardOutput $stdout -RedirectStandardError $stderr
  Write-Output "Started shard $id/$ShardCount (PID $($processes[-1].Id))."
}
Remove-Item Env:FSSGL_FAIR_SHARD_ID
Remove-Item Env:FSSGL_FAIR_SHARD_COUNT
foreach ($process in $processes) { $process.WaitForExit() }
$failed = @()
for ($id = 1; $id -le $ShardCount; $id++) {
  $processes[$id - 1].Refresh()
  $code = $processes[$id - 1].ExitCode
  $suffix = "_shard$($id)of$ShardCount"
  $manifest = Join-Path $projectRoot "data\processed\simulation\v3_fair\fssgl_v3_fair_manifest$suffix.rds"
  $checkpoint = Join-Path $projectRoot "data\processed\simulation\v3_fair\fssgl_v3_fair_running$suffix.rds"
  Write-Output "Shard $id/$ShardCount completed (exit code $code)."
  if (($null -ne $code -and $code -ne 0) -or
      -not (Test-Path -LiteralPath $manifest) -or
      (Test-Path -LiteralPath $checkpoint)) { $failed += $id }
}
if ($failed.Count) {
  throw "FSSGL shards failed: $($failed -join ', '). See $logDir."
}
