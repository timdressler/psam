# tid_psam_svm_analysis.py
#
# Performs SVM analysis, creates plots and exports data for further analysis.
#
# Tim Dressler, 17.04.2025

import os
import re
import time
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from tqdm import tqdm

# Set up paths
SCRIPTPATH = os.path.dirname(os.path.abspath(__file__))
expected_subpath = os.path.join('psam', 'analysis_script', 'python')
if re.search(re.escape(expected_subpath) + r'$', SCRIPTPATH):
    print('Path OK')
else:
    raise EnvironmentError('Path not OK')

MAINPATH = os.path.abspath(os.path.join(SCRIPTPATH, '..', '..'))
INPATH = os.path.join(MAINPATH, 'data', 'processed_data', 'svm_prepared_clean')
#INPATH = Path(INPATH)
OUTPATH = os.path.join(MAINPATH, 'data', 'analysis_data', 'svm_analysis')
# FUNPATH = os.path.join(MAINPATH, 'functions')  # Placeholder for addpath

os.makedirs(OUTPATH, exist_ok=True)
# Placeholder: tid_psam_check_folder_TD(MAINPATH, INPATH, OUTPATH)
# Placeholder: tid_psam_clean_up_folder_TD(OUTPATH)

# Variables to edit
INDIVIDUAL_PLOTS = True  # Whether or not to create individual ERP plots and Topoplots

# Set colors
def hex_to_rgb(color):
    color = color.lstrip('#')
    return tuple(int(color[i:i+2], 16) / 255 for i in (0, 2, 4))

main_blue = hex_to_rgb('#004F9F')
main_red = hex_to_rgb('#D53D0E')
main_green = hex_to_rgb('#00786B')
light_blue = hex_to_rgb('#5BC5F2')
main_yellow = hex_to_rgb('#FDC300')

# Get directory content
dircont_subj = [f for f in Path(INPATH).glob("sub*.csv")]

# Initialize sanity check variables
marked_subj = []
protocol = []

# Setup progress bar
print("Starting subject loop...")

# clear subj_idx
subj_idx = None
counter = 1
cor_counter = 1
con_counter = 1

for subj_idx, file in enumerate(tqdm(dircont_subj, desc="SVM Analysis")):

    # Get current ID
    subj = file.name
    print(file)
    print(subj)
    subj = re.search(r'sub-\d+', subj).group(0)

    # Update progress bar
    print(f"Processing {subj}...")

    tic = time.time()

    # Placeholder for main ERP logic

    # Update Protocol
    subj_time = time.time() - tic
    protocol.append([subj, subj_time, 'MARKED' if subj in [x[0] for x in marked_subj] else 'OK'])


# Get grandaverage ERPs (Placeholder)

# End of processing

protocol = pd.DataFrame(protocol, columns=['subj', 'time', 'status'])
protocol_path = os.path.join(OUTPATH, 'svm_analysis_protocol.xlsx')
protocol.to_excel(protocol_path, index=False)
quit()
print(f"Protocol saved to: {protocol_path}")

if marked_subj:
    marked_subj_df = pd.DataFrame(marked_subj, columns=['subj', 'issue'])
    marked_path = os.path.join(OUTPATH, 'tid_psam_svm_analysis_marked_subj.xlsx')
    marked_subj_df.to_excel(marked_path, index=False)
    print(f"Marked subjects saved to: {marked_path}")

check_done = 'tid_psam_svm_analysis_DONE'
print(check_done)

# Plots
plt.close('all')
