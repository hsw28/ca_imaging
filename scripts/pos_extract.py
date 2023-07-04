import os
import glob
import deeplabcut
from deeplabcut import analyze_videos

root_folder = r"C:\Users\HomePC\hswlab\Desktop\videos\DLC_cage_DEP"
config_file = os.path.join(root_folder, "config.yaml")

# Find all files starting with "concat" in the "My_WebCam" folders
search_pattern = os.path.join(root_folder, "My_WebCam", "*", "concat*")
file_paths = glob.glob(search_pattern)

# Loop through each file and run deeplabcut.analyze_videos()
for file_path in file_paths:
    # Replace 'FILEPATH' with the actual file path in the command
    command = f"deeplabcut.analyze_videos('{config_file}', ['{file_path}'], save_as_csv=True)"
    print(f"Running command: {command}")
    analyze_videos(config_file, [file_path], save_as_csv=True)
