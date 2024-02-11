# Set the root folder path where the search will start, be sure to replace online folder with hswlab
$rootFolder = "C:\Users\HomePC\hswlab\Desktop\videos\eyeblink\022223"
$webcamFolders = Get-ChildItem -Path $rootFolder -Filter "My_WebCam" -Directory -Recurse

foreach ($folder in $webcamFolders) {
    $aviFiles = Get-ChildItem -Path $folder.FullName -Filter "*.avi" | Sort-Object Name

    if ($aviFiles.Count -ge 5) {
        $videoFiles = $aviFiles.FullName

        Write-Host "Input Files for Concatenation in $($folder.FullName):"
        foreach ($file in $aviFiles) {
            Write-Host $file.Name
        }

        $creationDate = $folder.CreationTime.ToString("MMddyyyy")
        $outputFile = Join-Path -Path $folder.FullName -ChildPath "concat_webcam_$creationDate.avi"

        $concatInput = ($videoFiles | ForEach-Object { "file '$($_)'" }) -join [Environment]::NewLine
        $concatInputFile = [System.IO.Path]::GetTempFileName()
        $concatInput | Set-Content -Path $concatInputFile

        $concatArguments = "-f", "concat", "-safe", "0", "-i", "`"concat:$concatInputFile`"", "-c", "copy", $outputFile

        try {
            $process = Start-Process -FilePath "ffmpeg" -ArgumentList $concatArguments -NoNewWindow -PassThru -Wait -ErrorAction Stop
            Write-Host "Concatenation complete for folder: $($folder.FullName). Output file: $outputFile"

            # Get the frame count using ffprobe
            $ffprobeArguments = "-v", "error", "-count_frames", "-select_streams", "v:0", "-show_entries", "stream=nb_read_frames", "-of", "default=nokey=1:noprint_wrappers=1", $outputFile
            $frameCount = & "ffprobe" $ffprobeArguments

            Write-Host "Number of Frames: $frameCount"
        } catch {
            Write-Host "Error occurred while executing FFmpeg command for folder: $($folder.FullName)"
            Write-Host "Error message: $($_.Exception.Message)"
        }

        # Remove the temporary concat input file
        Remove-Item -Path $concatInputFile -Force
    }
}
