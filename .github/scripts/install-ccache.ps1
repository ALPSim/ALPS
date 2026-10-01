# Download the native Windows compiler cache; never build the tool in CI.
param(
    [Parameter(Mandatory)][ValidateSet('x64', 'arm64')][string]$Architecture,
    [Parameter(Mandatory)][string]$Destination
)
$ErrorActionPreference = 'Stop'
$packages = @{
    x64 = @{ Suffix = 'x86_64'; Sha256 = '2568347a697e103ca1b073981c704ad76fb2507d066c38dba038dd73399d968f' }
    arm64 = @{ Suffix = 'aarch64'; Sha256 = 'baa07085f68f75ea033d6f8d8fe6b53b7f353050bd191e299ff27e131cc3d787' }
}
$package = $packages[$Architecture]
$name = "ccache-4.14-windows-$($package.Suffix)"
New-Item -ItemType Directory -Force -Path $Destination | Out-Null
$archive = Join-Path $Destination "$name.zip"
Invoke-WebRequest "https://github.com/ccache/ccache/releases/download/v4.14/$name.zip" -OutFile $archive
if ((Get-FileHash -LiteralPath $archive -Algorithm SHA256).Hash -ne $package.Sha256) {
    throw 'ccache archive checksum mismatch'
}
Expand-Archive -LiteralPath $archive -DestinationPath $Destination -Force
$binaryDirectory = Join-Path (Resolve-Path -LiteralPath $Destination).Path $name
& "$binaryDirectory/ccache.exe" --version
if ($LASTEXITCODE) { throw 'ccache installation failed' }
if ($env:GITHUB_PATH) { $binaryDirectory >> $env:GITHUB_PATH }
Write-Output $binaryDirectory
