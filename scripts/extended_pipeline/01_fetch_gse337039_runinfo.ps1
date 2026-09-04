param(
    [string]$Study = "SRP712686",
    [string]$Output = (Join-Path $PSScriptRoot "..\..\..\05_数据与元数据\validation_metadata\metadata\GSE337039_SraRunInfo.csv")
)

$ErrorActionPreference = "Stop"
[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12

$searchUrl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=$Study&retmax=1000"
[xml]$search = (Invoke-WebRequest -Uri $searchUrl -UseBasicParsing -TimeoutSec 120).Content
$ids = @($search.eSearchResult.IdList.Id)
if ($ids.Count -ne 60) {
    throw "Expected 60 SRA records for $Study, found $($ids.Count)."
}

$idText = $ids -join ","
$fetchUrl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=$idText&rettype=runinfo&retmode=text"
$content = (Invoke-WebRequest -Uri $fetchUrl -UseBasicParsing -TimeoutSec 180).Content

$parent = Split-Path -Parent $Output
New-Item -ItemType Directory -Force -Path $parent | Out-Null
[IO.File]::WriteAllText($Output, $content, [Text.UTF8Encoding]::new($false))

$rows = Import-Csv -LiteralPath $Output
if ($rows.Count -ne 60) {
    throw "RunInfo validation failed: expected 60 rows, found $($rows.Count)."
}
if (($rows.LibraryLayout | Sort-Object -Unique) -ne "SINGLE") {
    throw "RunInfo validation failed: not all libraries are SINGLE."
}

$sizeMb = ($rows | Measure-Object -Property size_MB -Sum).Sum
$spots = ($rows | Measure-Object -Property spots -Sum).Sum
$bases = ($rows | Measure-Object -Property bases -Sum).Sum
$sha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $Output).Hash.ToLowerInvariant()

Write-Output "RUNINFO_OK"
Write-Output "runs=$($rows.Count)"
Write-Output "normalized_sra_size_mb=$sizeMb"
Write-Output "spots=$spots"
Write-Output "bases=$bases"
Write-Output "sha256=$sha256"
Write-Output "output=$Output"
