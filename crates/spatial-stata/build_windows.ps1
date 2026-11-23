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

$cargoArgs = @('rustc', '--color', 'never', '--release', '--locked', '--target', $target, '--target-dir', $targetDir, '--lib')
if ($env:CARGO_NET_OFFLINE -eq 'true') {
    $cargoArgs += '--offline'
}
$cargoArgs += @('--', '--print', 'native-static-libs')

$rustArchive = Join-Path $targetDir "$target/release/spatial_stata.lib"
$libraryPrefix = 'note: native-static-libs: '
$libraryLine = $null

Push-Location $PSScriptRoot
try {
    # Ask rustc for the system libraries required by this toolchain's static archive.
    for ($attempt = 0; $attempt -lt 2; $attempt++) {
        $prevPref = $ErrorActionPreference
        $ErrorActionPreference = 'Continue'
        $build = & cargo @cargoArgs 2>&1
        $exitCode = $LASTEXITCODE
        $ErrorActionPreference = $prevPref
        $build | ForEach-Object { Write-Host $_ }
        if ($exitCode -ne 0) {
            throw 'Building the Stata Rust archive failed.'
        }
        $libraryLine = $build | ForEach-Object { "$_" } | Where-Object { $_.StartsWith($libraryPrefix) } | Select-Object -Last 1
        if ($libraryLine) {
            break
        }
        # If the archive was already fresh, remove it so the second invocation prints native-static-libs.
        if ($attempt -eq 0 -and (Test-Path $rustArchive)) {
            Remove-Item $rustArchive -Force
        }
    }
} finally {
    Pop-Location
}

if (-not $libraryLine) {
    throw 'Rust did not report the libraries required to link the Stata plugin.'
}
$libraries = $libraryLine.Substring($libraryPrefix.Length).Trim() -split '\s+'
$source = Join-Path $PSScriptRoot 'c'
$tmpLinkDir = (New-Item -ItemType Directory -Force -Path (Join-Path $targetDir 'win-link-tmp')).FullName

Push-Location $tmpLinkDir
try {
    & cl.exe /nologo /LD /MD /O2 /DSYSTEM=STWIN32 "/I$source" "$source/stplugin.c" "$source/stata_bridge.c" $rustArchive @libraries /Fe:spatial_spur.plugin
    if ($LASTEXITCODE -ne 0) {
        throw 'Linking the Stata plugin failed.'
    }
    $exports = & dumpbin.exe /nologo /exports spatial_spur.plugin
    $exportedNames = $exports | ForEach-Object {
        if ($_ -match '^\s+\d+\s+[0-9a-fA-F]+\s+[0-9a-fA-F]+\s+(\w+)\s*$') { $Matches[1] }
    }
    if ($LASTEXITCODE -ne 0 -or ($exportedNames -join ' ') -ne 'pginit stata_call') {
        throw 'The Stata plugin must export pginit and stata_call.'
    }
    Copy-Item spatial_spur.plugin (Join-Path $output 'spatial_spur.plugin') -Force
} finally {
    Pop-Location
}
