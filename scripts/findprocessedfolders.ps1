$searchTerm = "My_V4_Miniscope"
$rootDir = "C:\Users\HomePC\OneDrive - Northwestern University\Desktop\videos\eyeblink\031423"
$folders = Get-ChildItem -Recurse -Directory -Path $rootDir | Where-Object { $_.FullName -match $searchTerm } | Sort-Object CreationTime

foreach ($folder in $folders) {
    $files = Get-ChildItem -Path $folder.FullName -File | Where-Object { $_.Name -match "sorted" }
    if ($files.Count -gt 0) {
        Write-Output "$($folder.FullName)"
    }
}
