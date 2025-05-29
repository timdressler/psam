import os
import sys
import re
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FuncFormatter
from pathlib import Path
from tqdm import tqdm

from sklearn.model_selection import StratifiedKFold
from sklearn.svm import SVC
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import accuracy_score

# Set up paths
SCRIPTPATH = os.path.dirname(os.path.abspath(__file__))
expected_subpath = os.path.join('psam', 'analysis_script', 'python')
if not re.search(re.escape(expected_subpath) + r'$', SCRIPTPATH):
    raise EnvironmentError('Path not OK')

MAINPATH = os.path.abspath(os.path.join(SCRIPTPATH, '..', '..'))
INPATH = os.path.join(MAINPATH, 'data', 'processed_data', 'svm_prepared_clean')
OUTPATH = os.path.join(MAINPATH, 'data', 'analysis_data', 'svm_analysis_TEST')

FUNPATH = os.path.join(MAINPATH, 'functions')
sys.path.append(FUNPATH)

# Load costum functions
from tid_psam_check_folder_clean_up_folder_TD import tid_psam_check_folder_TD, tid_psam_clean_up_folder_TD

tid_psam_check_folder_TD(MAINPATH, INPATH, OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH)
#os.makedirs(OUTPATH, exist_ok=True)

# Variables to edit
C_range_lower = -2
C_range_upper = 10
gamma_range_lower = -14
gamma_range_upper = -1
fixed_gamma_val = 1e-10
fixed_C_val = 1e4
grid_size = 3
sig_tresh = 0.55 # Threshold used for determined above-chance-level classification
n_splits = 5 # Number splits in k-fold CV

# Get directory content
dircont_subj_early = sorted([f for f in Path(INPATH).glob("sub-*early.csv")])
dircont_subj_late = sorted([f for f in Path(INPATH).glob("sub-*late.csv")])

# Create C-Range and Gamma-Range 
C_range = np.logspace(C_range_lower, C_range_upper, grid_size) 
gamma_range = np.logspace(gamma_range_lower, gamma_range_upper, grid_size) 

# Look for nearest value in the selected range to the selected fixed value (only relevant for plotting)
fixed_gamma_idx = np.argmin(np.abs(gamma_range - fixed_gamma_val))
fixed_C_idx = np.argmin(np.abs(C_range - fixed_C_val))
fixed_gamma_val = gamma_range[fixed_gamma_idx] # get gamma value (as selected one might be not exactly present in the gamma range)
fixed_C_val = C_range[fixed_C_idx] # get C value (as selected one might be not exactly present in the gamma range)

# Initialize variables
protocol = []
acc_collection_early, acc_collection_late = [], []
best_gamma_early_list, best_gamma_late_list = [], []
best_C_early_list, best_C_late_list = [], []

# Function: Performs stratified k-fold CV, PCA and fits SVMs for each subject and condition
def process_subject(file, condition):
    # Get current ID
    subj = re.search(r'sub-\d+', file.name).group(0)
    subj_outpath = os.path.join(OUTPATH, subj)
    # Get subject outpath
    os.makedirs(subj_outpath, exist_ok=True)
    # Load data
    df = pd.read_csv(file)
    # Extract features
    X = df.iloc[:, :-1].values
    # Extract outcome
    y = df.iloc[:, -1].values

    # Scale features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Fit 2D PCA (only for plotting)
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)

    # Plot: Individual 2D PCA projections including class labels
    y_encoded = LabelEncoder().fit_transform(y)
    plt.figure(figsize=(6, 5))
    plt.scatter(X_pca[:, 0], X_pca[:, 1], c=y_encoded, cmap='coolwarm', edgecolors='k', alpha=0.7)
    plt.title(f"{subj} - {condition} - PCA Projection")
    plt.xlabel("PC 1")
    plt.ylabel("PC 2")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(subj_outpath, f"{subj}_{condition}_pca_projection.png"), dpi=300)
    plt.close()

    # Initialize accuracy matrix
    acc_matrix = np.zeros((len(gamma_range), len(C_range)))
    # Prepare stratified k-fold CV
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=123)

    # Start grid search
    for i, gamma_val in enumerate(gamma_range):
        for j, C_val in enumerate(C_range):
            fold_accuracies = []
            for train_idx, test_idx in skf.split(X, y): # Get train and test indices
                # Scale features
                scaler = StandardScaler()
                X_train_scaled = scaler.fit_transform(X[train_idx])
                X_test_scaled = scaler.transform(X[test_idx])
                # Fit PCA (only fitted using the train split)
                pca = PCA(n_components=100)
                X_train_pca = pca.fit_transform(X_train_scaled)
                # Apply PCA weights also to the test split
                X_test_pca = pca.transform(X_test_scaled)
                # Fit SVM
                clf = SVC(C=C_val, kernel='rbf', gamma=gamma_val)
                clf.fit(X_train_pca, y[train_idx])
                # Get accuracy based on the test split 
                acc = accuracy_score(y[test_idx], clf.predict(X_test_pca))
                fold_accuracies.append(acc)
            acc_matrix[i, j] = np.mean(fold_accuracies)
            print(str(np.mean(fold_accuracies)))

    # Save accuracy grid
    acc_df = pd.DataFrame(acc_matrix,
                          index=[f"{g:.0e}" for g in gamma_range],
                          columns=[f"{c:.5f}" for c in C_range])
    acc_df.to_csv(os.path.join(subj_outpath, f"{subj}_{condition}_accuracy_grid.csv"))

    # Get individual best gamma and best C value (not further used)
    best_idx = np.unravel_index(np.argmax(acc_matrix), acc_matrix.shape)
    best_gamma = gamma_range[best_idx[0]]
    best_C = C_range[best_idx[1]]

    return acc_matrix, best_idx, subj_outpath, best_gamma, best_C

# Start processing
for file_early, file_late in tqdm(zip(dircont_subj_early, dircont_subj_late), total=len(dircont_subj_early), desc="Processing subjects"):
    # Get current ID
    subj = re.search(r'sub-\d+', file_early.name).group(0)
    tic = time.time()

    # Run process_subject function (see above)
    acc_early, best_idx_early, outpath, best_gamma_early, best_C_early = process_subject(file_early, 'early')
    acc_late, best_idx_late, _, best_gamma_late, best_C_late = process_subject(file_late, 'late')

    # Store values
    acc_collection_early.append(acc_early)
    acc_collection_late.append(acc_late)
    best_gamma_early_list.append(best_gamma_early)
    best_gamma_late_list.append(best_gamma_late)
    best_C_early_list.append(best_C_early)
    best_C_late_list.append(best_C_late)
    protocol.append([subj, time.time() - tic, 'OK'])

    # Plot: Individual accuracy heatmaps
    fig, axs = plt.subplots(1, 2, figsize=(16, 6), constrained_layout=True)
    for i, (acc_matrix, label) in enumerate(zip([acc_early, acc_late], ['early', 'late'])):
        acc_df = pd.DataFrame(acc_matrix, index=[f"{g:.0e}" for g in gamma_range], columns=[f"{c:.5f}" for c in C_range])
        sns.heatmap(acc_df, ax=axs[i], cmap='magma', vmin=0.5, vmax=0.95,
                    cbar=(i == 1), cbar_kws={"label": "Validation Accuracy"}, xticklabels=True, yticklabels=True)
        axs[i].set_title(f"{subj} - {label} - Accuracy Grid")
        axs[i].set_xlabel("C")
        axs[i].set_ylabel("Gamma")
        axs[i].tick_params(axis='x', rotation=90, labelsize=8)
        axs[i].tick_params(axis='y', labelsize=8)
    plt.savefig(os.path.join(outpath, f"{subj}_accuracy_heatmaps.png"), dpi=300)
    plt.close()

    # CHATGPT: Plot: Accuracy vs C (with fixed Gamma) for early and late
    plt.figure(figsize=(10, 5))
    plt.plot(C_range, acc_early[fixed_gamma_idx, :], label='Early', marker='o')
    plt.plot(C_range, acc_late[fixed_gamma_idx, :], label='Late', marker='s')
    plt.xscale('log')
    plt.xlabel("C")
    plt.ylabel("Accuracy")
    plt.title(f"{subj} - Accuracy vs C (Fixed Gamma: {gamma_range[fixed_gamma_idx]:.5e})")
    plt.xticks(C_range, [f"{c:.5f}" for c in C_range], rotation=90)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(outpath, f"{subj}_c_accuracy.png"), dpi=300)
    plt.close()

    # CHATGPT: Plot: Accuracy vs Gamma (with fixed C) for early and late
    plt.figure(figsize=(10, 5))
    plt.plot(gamma_range, acc_early[:, fixed_C_idx], label='Early', marker='o')
    plt.plot(gamma_range, acc_late[:, fixed_C_idx], label='Late', marker='s')
    plt.xscale('log')
    plt.xlabel("Gamma")
    plt.ylabel("Accuracy")
    plt.title(f"{subj} - Accuracy vs Gamma (Fixed C: {C_range[fixed_C_idx]:.2f})")
    plt.xticks(gamma_range, [f"{g:.0e}" for g in gamma_range], rotation=90)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(outpath, f"{subj}_gamma_accuracy.png"), dpi=300)
    plt.close()

# Get grandaverages of the accuracy matrices
ga_early = np.mean(acc_collection_early, axis=0)
ga_late = np.mean(acc_collection_late, axis=0)

# Plot: Grandaverage accuracy heatmaps
fig, axs = plt.subplots(1, 2, figsize=(16, 6), constrained_layout=True)
for i, (ga_matrix, label) in enumerate(zip([ga_early, ga_late], ['early', 'late'])):
    ga_df = pd.DataFrame(ga_matrix, index=[f"{g:.0e}" for g in gamma_range], columns=[f"{c:.5f}" for c in C_range])
    sns.heatmap(ga_df, ax=axs[i], cmap='magma', vmin=0.5, vmax=0.95,
                cbar=(i == 1), cbar_kws={"label": "Mean Accuracy"}, xticklabels=True, yticklabels=True)
    axs[i].set_title(f"Grand Average - {label.capitalize()}")
    axs[i].set_xlabel("C")
    axs[i].set_ylabel("Gamma")
    axs[i].tick_params(axis='x', rotation=90, labelsize=8)
    axs[i].tick_params(axis='y', labelsize=8)
plt.savefig(os.path.join(OUTPATH, "grandaverage_accuracy_heatmaps.png"), dpi=300)
plt.close()

# Plot: Grandaverage accuracy as a function of C (with fixed Gamma)
plt.figure(figsize=(10, 5))
plt.plot(C_range, ga_early[fixed_gamma_idx, :], label='Early', marker='o')
plt.plot(C_range, ga_late[fixed_gamma_idx, :], label='Late', marker='s')
plt.xscale('log')
plt.xlabel("C")
plt.ylabel("Accuracy")
plt.title(f"GA - Accuracy vs C (Fixed Gamma: {gamma_range[fixed_gamma_idx]:.5e})")
plt.xticks(C_range, [f"{c:.5f}" for c in C_range], rotation=90)
plt.tick_params(axis='both', labelsize=8)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(OUTPATH, "grandaverage_c_accuracy.png"), dpi=300)
plt.close()

# Plot: Grandaverage accuracy as a function of Gamma (with fixed C)
plt.figure(figsize=(10, 5))
plt.plot(gamma_range, ga_early[:, fixed_C_idx], label='Early', marker='o')
plt.plot(gamma_range, ga_late[:, fixed_C_idx], label='Late', marker='s')
plt.xscale('log')
plt.xlabel("Gamma")
plt.ylabel("Accuracy")
plt.title(f"GA - Accuracy vs Gamma (Fixed C: {C_range[fixed_C_idx]:.2f})")
plt.xticks(gamma_range, [f"{g:.0e}" for g in gamma_range], rotation=90)
plt.tick_params(axis='both', labelsize=8)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(OUTPATH, "grandaverage_gamma_accuracy.png"), dpi=300)
plt.close()

# Save accuracy data
pd.DataFrame(ga_early, index=[f"{g:.0e}" for g in gamma_range], columns=[f"{c:.5f}" for c in C_range]) \
    .to_excel(os.path.join(OUTPATH, 'grandaverage_accuracy_grid_early.xlsx'))
pd.DataFrame(ga_late, index=[f"{g:.0e}" for g in gamma_range], columns=[f"{c:.5f}" for c in C_range]) \
    .to_excel(os.path.join(OUTPATH, 'grandaverage_accuracy_grid_late.xlsx'))

# Get proportion of the possible hyperparameter pairs that resulted in a above-chance-level classification for both early and late conditions
sig_props = []
for subj_idx, subj in enumerate([re.search(r'sub-\d+', f.name).group(0) for f in dircont_subj_early]):
    acc_early = acc_collection_early[subj_idx]
    acc_late = acc_collection_late[subj_idx]
    prop_early = np.mean(acc_early > sig_tresh)
    prop_late = np.mean(acc_late > sig_tresh)
    sig_props.append([subj, prop_early, prop_late])

# Save proportion data
sig_prop_df = pd.DataFrame(sig_props, columns=['subj', 'prop_sig_early', 'prop_sig_late'])
sig_prop_df.to_excel(os.path.join(OUTPATH, 'all_subj_accuracy_proportions.xlsx'), index=False)

# End of processing
protocol_df = pd.DataFrame(protocol, columns=['subj', 'time', 'status'])
protocol_df.to_excel(os.path.join(OUTPATH, 'svm_analysis_protocol.xlsx'), index=False)

check_done = "tid_psam_svm_analysis_DONE"

print(check_done)
