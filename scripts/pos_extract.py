import os
import glob
from deeplabcut import analyze_videos

root_folder = r"C:\Users\HomePC\hswlab\Desktop\videos\eyeblink\031423\2023_05_25"
config_file = r"C:\Users\HomePC\hswlab\Desktop\videos\DLC_oval\config.yaml"
#for square: C:\Users\HomePC\hswlab\Desktop\videos\DLC_cage\config.yaml"
#for oval: "C:\Users\HomePC\hswlab\Desktop\videos\DLC_oval\config.yaml"

# Find all folders named "My_WebCam" recursively within the root folder
search_pattern = os.path.join(root_folder, "**", "My_WebCam")
my_webcam_folders = glob.glob(search_pattern, recursive=True)

# Print the available My_WebCam folders
print("Available My_WebCam folders:")
for i, folder in enumerate(my_webcam_folders):
    print(f"{i + 1}. {folder}")

# Prompt user for selection
while True:
    try:
        selection = input("Select the My_WebCam folders you want to search (comma-separated numbers): ")
        selected_indices = [int(index) - 1 for index in selection.split(",")]
        selected_folders = [my_webcam_folders[index] for index in selected_indices]
        break
    except (ValueError, IndexError):
        print("Invalid selection. Please try again.")

# List to store selected file paths
selected_files = []

# Loop through each selected folder
for folder in selected_folders:
    # Find all files starting with "concat" within the current My_WebCam folder
    search_pattern = os.path.join(folder, "concat*")
    file_paths = glob.glob(search_pattern)

    # Print the files found in the current My_WebCam folder
    print(f"Files found in {folder}:")
    for i, file_path in enumerate(file_paths):
        print(f"{i + 1}. {file_path}")

    # Add the files to the selected_files list
    selected_files.extend(file_paths)

# Loop through each selected file and run deeplabcut.analyze_videos()
for file_path in selected_files:
    # Replace 'FILEPATH' with the actual file path in the command
    command = f"deeplabcut.analyze_videos('{config_file}', ['{file_path}'], save_as_csv=True)"
    print(f"Running command: {command}")
    analyze_videos(config_file, [file_path], save_as_csv=True)
