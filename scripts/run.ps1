param(
  [Parameter(Position = 0)]
  [string]$Task = "help"
)

$ErrorActionPreference = "Stop"
$ProjectRoot = Split-Path -Parent $PSScriptRoot

# R for Windows does not support the Unix-style C.UTF-8 locale name. Clear it
# only for this child process; this does not change the user's global settings.
foreach ($VariableName in @("LANG", "LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME")) {
  if ([Environment]::GetEnvironmentVariable($VariableName) -eq "C.UTF-8") {
    Remove-Item "Env:$VariableName"
  }
}

$RscriptCommand = Get-Command Rscript -ErrorAction SilentlyContinue
if ($RscriptCommand) {
  $RscriptExe = $RscriptCommand.Source
} else {
  $Candidates = Get-ChildItem "C:\Program Files\R\R-*\bin\Rscript.exe" -ErrorAction SilentlyContinue |
    Sort-Object FullName -Descending
  if (-not $Candidates) {
    throw "Rscript was not found on PATH or under C:\Program Files\R."
  }
  $RscriptExe = $Candidates[0].FullName
}

& $RscriptExe (Join-Path $ProjectRoot "scripts\run.R") $Task
exit $LASTEXITCODE
