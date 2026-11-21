param(
    [Parameter(Mandatory = $true)]
    [string] $OutputDir
)

$ErrorActionPreference = 'Stop'
$target = 'x86_64-pc-windows-msvc'
$targetDir = Join-Path $PSScriptRoot '../../target-stata'
$output = (New-Item -ItemType Directory -Force -Path $OutputDir).FullName

if ($env:VSCMD_ARG_TGT_ARCH -ne 'x64') {
    $vswhere = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio/Installer/vswhere.exe'
    $install = & $vswhere -latest -products '*' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
    if (-not $install) {
        throw 'Install the Visual Studio C++ build tools for x64 before building the Stata plugins.'
    }
    Import-Module (Join-Path $install 'Common7/Tools/Microsoft.VisualStudio.DevShell.dll')
    Enter-VsDevShell -VsInstallPath $install -SkipAutomaticLocation -DevCmdArguments '-arch=x64 -host_arch=x64'
}

Push-Location $PSScriptRoot
try {
    # Ask rustc for the system libraries required by this toolchain's static archive.
    $build = & cargo rustc --color never --release --locked --offline --target $target --target-dir $targetDir --lib -- --print native-static-libs 2>&1
    $build | ForEach-Object { Write-Host $_ }
    if ($LASTEXITCODE -ne 0) {
        throw 'Building the Stata Rust archive failed.'
    }
} finally {
    Pop-Location
}

$libraryPrefix = 'note: native-static-libs: '
$libraryLine = $build | ForEach-Object { "$_" } | Where-Object { $_.StartsWith($libraryPrefix) } | Select-Object -Last 1
if (-not $libraryLine) {
    throw 'Rust did not report the libraries required to link the Stata plugin.'
}
$libraries = $libraryLine.Substring($libraryPrefix.Length).Trim() -split '\s+'
$rustArchive = Join-Path $targetDir "$target/release/spatial_stata.lib"
$source = Join-Path $PSScriptRoot 'c'

Push-Location $output
try {
    & cl.exe /nologo /LD /MD /O2 /DSYSTEM=STWIN32 "/I$source" "$source/stplugin.c" "$source/stata_bridge.c" $rustArchive @libraries /Fe:scpc_rust.plugin
    if ($LASTEXITCODE -ne 0) {
        throw 'Linking the Stata plugin failed.'
    }
    $exports = & dumpbin.exe /nologo /exports scpc_rust.plugin
    $exportedNames = $exports | ForEach-Object {
        if ($_ -match '^\s+\d+\s+[0-9a-fA-F]+\s+[0-9a-fA-F]+\s+(\w+)\s*$') { $Matches[1] }
    }
    if ($LASTEXITCODE -ne 0 -or ($exportedNames -join ' ') -ne 'pginit stata_call') {
        throw 'The Stata plugin must export pginit and stata_call.'
    }
    Copy-Item scpc_rust.plugin spur_rust.plugin -Force
} finally {
    Pop-Location
}
