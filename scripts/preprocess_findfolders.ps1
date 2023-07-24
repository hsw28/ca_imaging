$searchTerm = "My_V4_Miniscope"
$rootDir = "C:\Users\HomePC\OneDrive - Northwestern University\Desktop\videos\eyeblink\031423"
$folders = Get-ChildItem -Recurse -Directory -Path $rootDir | Where-Object { $_.FullName -match $searchTerm } | Sort-Object CreationTime

foreach ($folder in $folders) {
    $folderPath = $folder.FullName
    $matFiles = Get-ChildItem -Path $folderPath -Filter "*.mat" -File -ErrorAction SilentlyContinue

    if (-not $matFiles) {
        $itemCount = (Get-ChildItem -Path $folderPath -File | Measure-Object).Count
        if ($itemCount -ge 5) {
            $aviFiles = Get-ChildItem -Path $folderPath -Filter "*.avi" -File
            foreach ($aviFile in $aviFiles) {
                if ($aviFile.Name -notlike "raw*.avi") {
                    $newFileName = "raw" + $aviFile.Name
                    $newFilePath = Join-Path -Path $folderPath -ChildPath $newFileName
                    $aviFile | Rename-Item -NewName $newFileName -Force
                    Write-Host "Renamed file: $aviFile to $newFileName"
                }
            }
            Write-Host "Folder found: $folderPath"
        }
    }
}
