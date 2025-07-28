import os
import sys
import re
import time
import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FuncFormatter
from pathlib import Path
from tqdm import tqdm
from scipy.stats import binom
from sklearn.model_selection import StratifiedKFold
from sklearn.svm import SVC
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import accuracy_score
from itertools import product

"""
tid_psam_svm_analysis.py

Performs subject-level classification using
SVMs with grid search and stratified k-fold cross-validation.

Processing includes:
- Loading early and late time-window features for each subject
- Feature scaling and PCA
- Grid search across C and gamma hyperparameters
- Accuracy evaluation and statistical thresholding
- Generation of per-subject and grand-average accuracy plots
- Export of accuracy grids, plots, and summary metrics

Tim Dressler, 07.06.2025
"""

random.seed(123)

# Set up paths
SCRIPTPATH = os.path.dirname(os.path.abspath(__file__))
expected_subpath = os.path.join('psam', 'analysis_script', 'python')

if not SCRIPTPATH.endswith(expected_subpath):
    raise EnvironmentError('Path not OK')

MAINPATH = os.path.abspath(os.path.join(SCRIPTPATH, '..', '..'))
INPATH = os.path.join(MAINPATH, 'data', 'processed_data', 'svm_prepared_clean')
OUTPATH = os.path.join(MAINPATH, 'data', 'analysis_data', 'svm_analysis')
FUNPATH = os.path.join(MAINPATH, 'functions')

sys.path.append(FUNPATH)

# Load costum functions
from tid_psam_check_folder_clean_up_folder_TD import tid_psam_check_folder_TD, tid_psam_clean_up_folder_TD # type: ignore

tid_psam_check_folder_TD(MAINPATH, INPATH, OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH)

# Variables to edit
C_RANGE_LOWER = -2 # Set to -2
C_RANGE_UPPER = 10 # Set to 10
GAMMA_RANGE_LOWER = -14 # Set to -14
GAMMA_RANGE_UPPER = -1 # Set to -1
FIXED_GAMMA_VAL = 1e-10 # Used to plot validation accuracy as a function of C (with a fixed value for gamma) # Set to 1e-10
FIXED_C_VAL = 1e4 # Used to plot validation accuracy as a function of gamma (with a fixed value for C) # Set to 1e4
GRID_SIZE = 30 # Number of values used for gamma and C, respectively # Set to 30
#sig_thresh = 0.55 # Threshold used for determined above-chance-level classification (not further used, individual tresholds used)
N_SPLITS = 5 # Number splits in k-fold CV # Set to 5
N_COMPS = 100 # Number of components extracted by the PCA # Set to 100

# Get directory content
dircont_subj_early = sorted([f for f in Path(INPATH).glob("sub-*early.csv")])
dircont_subj_late = sorted([f for f in Path(INPATH).glob("sub-*late.csv")])

# Create C-Range and Gamma-Range 
C_range = np.logspace(C_RANGE_LOWER, C_RANGE_UPPER, GRID_SIZE) 
gamma_range = np.logspace(GAMMA_RANGE_LOWER, GAMMA_RANGE_UPPER, GRID_SIZE) 

# Look for nearest value in the selected range to the selected fixed value (only relevant for plotting)
fixed_gamma_idx = np.argmin(np.abs(gamma_range - FIXED_GAMMA_VAL))
fixed_C_idx = np.argmin(np.abs(C_range - FIXED_C_VAL))
FIXED_GAMMA_VAL = gamma_range[fixed_gamma_idx] # get gamma value (as selected one might be not exactly present in the gamma range)
FIXED_C_VAL = C_range[fixed_C_idx] # get C value (as selected one might be not exactly present in the gamma range)

# Initialize variables
protocol = []
all_acc_matrix_early, all_acc_matrix_late = [], []
all_best_gamma_early, all_best_gamma_late = [], []
all_best_c_early, all_best_c_late = [], []
all_sig_thresh = []
all_n_trials = []

# Function: Compute statistical threshold (chance level) for classification accuracy
def compute_statistical_threshold(y, alpha=0.05):
    """
    Compute statistical threshold (chance level) for classification accuracy using inverse binomial CDF (Combrisson & Jerbi, 2015).
    
    Parameters:
        y (array-like): Class labels
        alpha (float): Significance level (e.g., 0.05)

    Returns:
        float: Threshold accuracy (e.g., 0.5 (not 50%!)) that is statistically above chance

    Literature
        Combrisson, E., & Jerbi, K. (2015). Exceeding chance level by chance: The caveat
            of theoretical chance levels in brain signal classification and statistical
            assessment of decoding accuracy. Journal of Neuroscience Methods, 250,
            126–136. https://doi.org/10.1016/j.jneumeth.2015.01.010
    """
    n_trials = len(y)
    n_classes = len(np.unique(y))
    threshold_count = binom.ppf(1 - alpha, n_trials, 1 / n_classes)
    threshold_percent = threshold_count / n_trials
    return threshold_percent

# Function: Performs grid search using stratified k-fold CV, PCA and fits SVMs for each subject and condition
def process_subject(file, condition):
    """
    Performs grid search using stratified k-fold CV, PCA and fits SVMs for each subject and condition.

    Parameters:
        file (Path): Path to the subject's CSV data file containing features and labels.
        condition (str): Condition label (e.g., 'early' or 'late') used for file naming and output.

    Returns:
        acc_matrix (np.ndarray): Grid of mean classification accuracies (gamma × C).
        best_idx (tuple): Indices of best-performing gamma and C in the accuracy matrix.
        subj_outpath (str): Path to the subject-specific output folder.
        best_gamma (float): Gamma value achieving the highest accuracy.
        best_C (float): C value achieving the highest accuracy.
        sig_thresh (float): Subject-specific statistical threshold (accuracy %) for above-chance classification (see compute_statistical_threshold).
        n_trials (float): Subject-specific number of included trials
    """
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

    # Get number of trials
    n_trials = len(y)

    # Get individual threshold for above-chance classification
    sig_thresh = compute_statistical_threshold(y, alpha=0.05)

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
    skf = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=123)

    # Start grid search
    grid_iter = list(product(enumerate(gamma_range), enumerate(C_range)))
    for (i, gamma_val), (j, C_val) in tqdm(grid_iter, desc=f"{subj} - {condition} Grid Search", leave=False):
        fold_accuracies = []
        for train_idx, test_idx in skf.split(X, y):
            # Scale features
            scaler = StandardScaler()
            X_train_scaled = scaler.fit_transform(X[train_idx])
            X_test_scaled = scaler.transform(X[test_idx])
            # PCA
            pca = PCA(n_components=N_COMPS)
            X_train_pca = pca.fit_transform(X_train_scaled)
            X_test_pca = pca.transform(X_test_scaled)
            # Train SVM
            clf = SVC(C=C_val, kernel='rbf', gamma=gamma_val)
            clf.fit(X_train_pca, y[train_idx])
            acc = accuracy_score(y[test_idx], clf.predict(X_test_pca))
            fold_accuracies.append(acc)
        acc_matrix[i, j] = np.mean(fold_accuracies)
        #print('Fold Accuracy: ' + str(np.mean(fold_accuracies)))

    # Save accuracy grid
    acc_df = pd.DataFrame(acc_matrix,
                          index=[f"{g:.0e}" for g in gamma_range],
                          columns=[f"{c:.5f}" for c in C_range])
    acc_df.to_csv(os.path.join(subj_outpath, f"{subj}_{condition}_accuracy_grid.csv"))

    # Get individual best gamma and best C value (not further used)
    best_idx = np.unravel_index(np.argmax(acc_matrix), acc_matrix.shape)
    best_gamma = gamma_range[best_idx[0]]
    best_C = C_range[best_idx[1]]

    return acc_matrix, best_idx, subj_outpath, best_gamma, best_C, sig_thresh, n_trials

# Start processing
for file_early, file_late in tqdm(zip(dircont_subj_early, dircont_subj_late), total=len(dircont_subj_early), desc="Processing subjects"):
    # Get current ID
    subj = re.search(r'sub-\d+', file_early.name).group(0)
    tic = time.time()

    # Run process_subject function (see above)
    acc_early, best_idx_early, outpath, best_gamma_early, best_C_early, sig_thresh_early, n_trials = process_subject(file_early, 'early')
    acc_late, best_idx_late, _, best_gamma_late, best_C_late, sig_thresh_late, n_trials = process_subject(file_late, 'late')

    # Sanity Check: Same threshold for early and late time-window (has to be the case since n_trials has to be the same)
    if sig_thresh_early == sig_thresh_late:
        sig_thresh = sig_thresh_early
    else:
        raise ValueError('Non-matching significance thresholds!')

    # Store values
    all_acc_matrix_early.append(acc_early)
    all_acc_matrix_late.append(acc_late)
    all_best_gamma_early.append(best_gamma_early)
    all_best_gamma_late.append(best_gamma_late)
    all_best_c_early.append(best_C_early)
    all_best_c_late.append(best_C_late)
    all_sig_thresh.append(sig_thresh)
    all_n_trials.append(n_trials)
    protocol.append([subj, time.time() - tic, 'OK'])

    # Plot: Individual accuracy heatmaps
    fig, axs = plt.subplots(1, 2, figsize=(16, 6), constrained_layout=True)
    for i, (acc_matrix, label) in enumerate(zip([acc_early, acc_late], ['early', 'late'])):
        acc_df = pd.DataFrame(acc_matrix, index=[f"{g:.0e}" for g in gamma_range], columns=[f"{c:.5f}" for c in C_range])
        min_val = acc_matrix.min()
        max_val = acc_matrix.max()
        above_chance_prop = np.mean(acc_matrix > sig_thresh) 

        sns.heatmap(acc_df, ax=axs[i], cmap='magma', vmin=0.5, vmax=0.85,
                    cbar=(i == 1), cbar_kws={"label": "Validation Accuracy"}, xticklabels=True, yticklabels=True)
        axs[i].set_title(f"{subj} - {label.capitalize()}\n"
                        f"Min: {min_val:.3f}, Max: {max_val:.3f}, >Chance: {above_chance_prop:.2%}")
        axs[i].set_xlabel("C")
        axs[i].set_ylabel("Gamma")
        axs[i].tick_params(axis='x', rotation=90, labelsize=8)
        axs[i].tick_params(axis='y', labelsize=8)
    plt.savefig(os.path.join(outpath, f"{subj}_accuracy_heatmaps.png"), dpi=1200)
    plt.close()


    # Plot: Accuracy vs C (with fixed Gamma) for early and late
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

    # Plot: Accuracy vs Gamma (with fixed C) for early and late
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
ga_early = np.mean(all_acc_matrix_early, axis=0)
ga_late = np.mean(all_acc_matrix_late, axis=0)

# Plot: Grandaverage accuracy heatmaps
fig, axs = plt.subplots(1, 2, figsize=(16, 6), constrained_layout=True)
for i, (ga_matrix, label) in enumerate(zip([ga_early, ga_late], ['early', 'late'])):
    ga_df = pd.DataFrame(ga_matrix, index=[f"{g:.0e}" for g in gamma_range], columns=[f"{c:.5f}" for c in C_range])
    sns.heatmap(ga_df, ax=axs[i], cmap='magma', vmin=0.5, vmax=0.85,
                cbar=(i == 1), cbar_kws={"label": "Validation Accuracy"}, xticklabels=True, yticklabels=True)
    axs[i].set_title(f"Grand Average - {label.capitalize()}")
    axs[i].set_xlabel("C")
    axs[i].set_ylabel("Gamma")
    axs[i].tick_params(axis='x', rotation=90, labelsize=8)
    axs[i].tick_params(axis='y', labelsize=8)
plt.savefig(os.path.join(OUTPATH, "tid_psam_grandaverage_accuracy_heatmaps.png"), dpi=1200)
plt.close()

# Plot: Grandaverage accuracy as a function of C (with fixed Gamma)
plt.figure(figsize=(10, 5))
plt.plot(C_range, ga_early[fixed_gamma_idx, :], label='Early', marker='o')
plt.plot(C_range, ga_late[fixed_gamma_idx, :], label='Late', marker='s')
plt.xscale('log')
plt.xlabel("C")
plt.ylabel("Validation Accuracy")
plt.title(f"GA - Accuracy vs C (Fixed Gamma: {gamma_range[fixed_gamma_idx]:.5e})")
plt.xticks(C_range, [f"{c:.5f}" for c in C_range], rotation=90)
plt.tick_params(axis='both', labelsize=8)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(OUTPATH, "tid_psam_grandaverage_c_accuracy.png"), dpi=900)
plt.close()

# Plot: Grandaverage accuracy as a function of Gamma (with fixed C)
plt.figure(figsize=(10, 5))
plt.plot(gamma_range, ga_early[:, fixed_C_idx], label='Early', marker='o')
plt.plot(gamma_range, ga_late[:, fixed_C_idx], label='Late', marker='s')
plt.xscale('log')
plt.xlabel("Gamma")
plt.ylabel("Validation Accuracy")
plt.title(f"GA - Accuracy vs Gamma (Fixed C: {C_range[fixed_C_idx]:.2f})")
plt.xticks(gamma_range, [f"{g:.0e}" for g in gamma_range], rotation=90)
plt.tick_params(axis='both', labelsize=8)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(OUTPATH, "tid_psam_grandaverage_gamma_accuracy.png"), dpi=900)
plt.close()

# Save accuracy data
pd.DataFrame(ga_early, index=[f"{g:.0e}" for g in gamma_range], columns=[f"{c:.5f}" for c in C_range]) \
    .to_excel(os.path.join(OUTPATH, 'grandaverage_accuracy_grid_early.xlsx'))
pd.DataFrame(ga_late, index=[f"{g:.0e}" for g in gamma_range], columns=[f"{c:.5f}" for c in C_range]) \
    .to_excel(os.path.join(OUTPATH, 'grandaverage_accuracy_grid_late.xlsx'))

# Get proportion of the possible hyperparameter pairs that resulted in a above-chance-level classification for both early and late conditions
sig_props = []
for subj_idx, subj in enumerate([re.search(r'sub-\d+', f.name).group(0) for f in dircont_subj_early]):
    sig_thresh = all_sig_thresh[subj_idx]
    acc_early = all_acc_matrix_early[subj_idx]
    acc_late = all_acc_matrix_late[subj_idx]
    prop_early = np.mean(acc_early > sig_thresh)
    prop_late = np.mean(acc_late > sig_thresh)
    n_trials = all_n_trials[subj_idx]
    sig_props.append([subj, prop_early, prop_late, sig_thresh, n_trials])

# Save proportion data
sig_prop_df = pd.DataFrame(sig_props, columns=['subj', 'prop_sig_early', 'prop_sig_late', 'sig_thresh', 'n_trials'])
sig_prop_df.to_excel(os.path.join(OUTPATH, 'all_subj_accuracy_proportions.xlsx'), index=False)

# End of processing
protocol_df = pd.DataFrame(protocol, columns=['subj', 'time', 'status'])
protocol_df.to_excel(os.path.join(OUTPATH, 'tid_psam_svm_analysis_protocol.xlsx'), index=False)

check_done = "tid_psam_svm_analysis_DONE"
print(check_done)
