import os
import shutil
from datetime import datetime

# Python equivalent of tid_psam_check_folder.m

def tid_psam_check_folder_TD(*folder_paths):
    """
    tid_psam_check_folder_TD - Checks if the specified folder paths exist. If a path does not exist, the function creates it.

    Usage:
        tid_psam_check_folder_TD('Path1', 'Path2', ...)

    Inputs:
        *folder_paths - any number of folder paths as strings

    Tim Dressler, 29.05.2024
    """
    for folder_path in folder_paths:
        if not os.path.exists(folder_path):
            print(f'The folder "{folder_path}" does not exist.')
            try:
                os.makedirs(folder_path)
                print(f'Folder "{folder_path}" was successfully created.')
            except Exception as e:
                print(f'Error: Could not create folder "{folder_path}". Exception: {e}')
        else:
            print(f'Folder "{folder_path}" exists.')

# Python equivalent of tid_psam_clean_up_folder.m

def tid_psam_clean_up_folder_TD(folder_path):
    """
    tid_psam_clean_up_folder_TD - Moves files and subfolders to a timestamped archive subfolder within an "archive" folder

    Usage:
        tid_psam_clean_up_folder_TD('Path')

    Input:
        folder_path - The path of the folder to check and clean up.

    Tim Dressler, 29.05.2024 (Updated)
    """
    # Check if the folder exists
    if not os.path.isdir(folder_path):
        raise FileNotFoundError("Specified folder does not exist")

    # List everything in the folder except 'archive' and hidden files
    items = [f for f in os.listdir(folder_path)
             if f not in {'.', '..', 'archive'} and not f.startswith('.')]
    
    # If there's anything to archive
    if items:
        archive_folder_path = os.path.join(folder_path, 'archive')
        os.makedirs(archive_folder_path, exist_ok=True)

        timestamp = datetime.now().strftime("%d_%m_%Y_%H-%M-%S")
        timestamped_archive_folder = os.path.join(archive_folder_path, f"archive_{timestamp}")
        os.makedirs(timestamped_archive_folder)

        # Move each item (file or folder) to the archive
        for item in items:
            src = os.path.join(folder_path, item)
            dst = os.path.join(timestamped_archive_folder, item)
            shutil.move(src, dst)

        print(f'All contents moved to: {timestamped_archive_folder}')
    else:
        print('The folder is empty. No action taken.')
