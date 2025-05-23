import os
import re
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from tqdm import tqdm

from sklearn.model_selection import StratifiedKFold
from sklearn.svm import SVC
from sklearn.decomposition import KernelPCA, PCA
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import accuracy_score

# Set up and validate script path
SCRIPTPATH = os.path.dirname(os.path.abspath(__file__))
expected_subpath = os.path.join('psam', 'analysis_script', 'python')
if re.search(re.escape(expected_subpath) + r'$', SCRIPTPATH):
    print('Path OK')
else:
    raise EnvironmentError('Path not OK')

# Define input and output paths
MAINPATH = os.path.abspath(os.path.join(SCRIPTPATH, '..', '..'))
INPATH = os.path.join(MAINPATH, 'data', 'processed_data', 'svm_prepared_clean')
OUTPATH = os.path.join(MAINPATH, 'data', 'analysis_data', 'svm_analysis')
os.makedirs(OUTPATH, exist_ok=True)

# Subject files
dircont_subj_early = [f for f in Path(INPATH).glob("sub-*early.csv")]

# Protocol
protocol = []

# Hyperparameter ranges (narrowed)
C_range = np.logspace(-2, 2, 10)
gamma_range = np.logspace(-3, 1, 10)
n_splits = 5  # Stratified CV

for subj_idx, file in enumerate(tqdm(dircont_subj_early, desc="SVM Analysis")):
    subj = file.name
    subj = re.search(r'sub-\d+', subj).group(0)
    tic = time.time()
    print(subj)

    df = pd.read_csv(file)
    X = df.iloc[:, :-1].values
    y = df.iloc[:, -1].values
    print('read file...')

    ###### Classwise mean/std print
    label_col = df.columns[-1]
    means_by_class = df.groupby(label_col).mean()
    sd_by_class = df.groupby(label_col).std()
    print(means_by_class)
    print(sd_by_class)
    ######

    # Sanity Check: Label distribution
    unique_labels, label_counts = np.unique(y, return_counts=True)
    print(f"Sanity Check: Label counts: {dict(zip(unique_labels, label_counts))}")

    # === PCA Visualization ===
    scaler_vis = StandardScaler()
    X_scaled_vis = scaler_vis.fit_transform(X)

    pca_vis = PCA(n_components=2)
    X_pca_vis = pca_vis.fit_transform(X_scaled_vis)

    # Convert string labels to numeric for coloring
    y_encoded = LabelEncoder().fit_transform(y)

    plt.figure(figsize=(6, 5))
    plt.scatter(X_pca_vis[:, 0], X_pca_vis[:, 1], c=y_encoded, cmap='coolwarm', edgecolors='k', alpha=0.7)
    plt.title(f"{subj} - PCA Projection", fontsize=12)
    plt.xlabel("PC 1")
    plt.ylabel("PC 2")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPATH, f"{subj}_pca_projection.png"), dpi=300)
    plt.close()

    # === Manual Feature Plot: mean_E01 vs band_power_14_30_E01 ===
    if 'mean_E01' in df.columns and 'band_power_14_30_E01' in df.columns:
        x_manual = df['mean_E01'].values
        y_manual = df['band_power_14_30_E01'].values

        plt.figure(figsize=(6, 5))
        plt.scatter(x_manual, y_manual, c=y_encoded, cmap='coolwarm', edgecolors='k', alpha=0.7)
        plt.xlabel("mean_E01")
        plt.ylabel("band_power_14_30_E01")
        plt.title(f"{subj} - Feature Space: mean_E01 vs band_power_14_30_E01", fontsize=12)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(OUTPATH, f"{subj}_manual_features_plot.png"), dpi=300)
        plt.close()
    else:
        print(f"Warning: One or both of mean_E01 and band_power_14_30_E01 not found in {subj}!")

    acc_matrix = np.zeros((len(gamma_range), len(C_range)))
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=123)

    print('start grid search...')
    for i, gamma_val in enumerate(gamma_range):
        for j, C_val in enumerate(C_range):
            fold_accuracies = []

            for fold_idx, (train_idx, test_idx) in enumerate(skf.split(X, y)):
                X_train, X_test = X[train_idx], X[test_idx]
                y_train, y_test = y[train_idx], y[test_idx]

                # Sanity Check: Per-fold label distribution
                u_train, c_train = np.unique(y_train, return_counts=True)
                u_test, c_test = np.unique(y_test, return_counts=True)
                print(f"Sanity Check: Fold {fold_idx} train labels: {dict(zip(u_train, c_train))}")
                print(f"Sanity Check: Fold {fold_idx} test labels: {dict(zip(u_test, c_test))}")

                # Scale and apply PCA
                scaler = StandardScaler()
                X_train_scaled = scaler.fit_transform(X_train)
                X_test_scaled = scaler.transform(X_test)

                # PCA (or KernelPCA alternative)
                pca = PCA(n_components=180)
                X_train_pca = pca.fit_transform(X_train_scaled)
                X_test_pca = pca.transform(X_test_scaled)

                explained_var = pca.explained_variance_ratio_.sum()
                print(f"Total variance retained: {explained_var:.2f}")

                # SVM training and evaluation
                clf = SVC(C=C_val, kernel='rbf', gamma=gamma_val)
                clf.fit(X_train_pca, y_train)
                y_pred = clf.predict(X_test_pca)

                # Sanity Check: Prediction distribution
                u_pred, c_pred = np.unique(y_pred, return_counts=True)
                print(f"Sanity Check: Fold {fold_idx} predictions: {dict(zip(u_pred, c_pred))}")

                acc = accuracy_score(y_test, y_pred)
                fold_accuracies.append(acc)

            acc_matrix[i, j] = np.mean(fold_accuracies)
            print(f"Mean accuracy (γ={gamma_val:.4f}, C={C_val:.4f}): {acc_matrix[i, j]:.4f}")

    subj_time = time.time() - tic
    protocol.append([subj, subj_time, 'OK', np.mean(acc_matrix)])

    # Save accuracy grid
    acc_df = pd.DataFrame(acc_matrix,
                          index=[f"{g:.5f}" for g in gamma_range],
                          columns=[f"{c:.5f}" for c in C_range])
    acc_path = os.path.join(OUTPATH, f"{subj}_acc_grid.csv")
    acc_df.to_csv(acc_path)

    # Heatmap plot
    plt.figure(figsize=(8, 6))
    im = plt.imshow(acc_matrix, origin='lower', aspect='auto',
                    extent=[C_range[0], C_range[-1], gamma_range[0], gamma_range[-1]],
                    cmap='magma', vmin=0.5, vmax=0.95)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("C", fontsize=12)
    plt.ylabel("Gamma", fontsize=12)
    plt.title(f"Validation Accuracy Grid: {subj}", fontsize=14)
    cbar = plt.colorbar(im)
    cbar.set_label('Mean Accuracy', fontsize=12)

    # Add min and max accuracy annotation
    min_acc = acc_matrix.min()
    max_acc = acc_matrix.max()
    text_str = f"Min Acc: {min_acc:.3f}\nMax Acc: {max_acc:.3f}"
    plt.text(0.05, 0.95, text_str, transform=plt.gca().transAxes,
             fontsize=10, color='white', ha='left', va='top',
             bbox=dict(boxstyle='round,pad=0.3', facecolor='black', alpha=0.5))

    plt.tight_layout()
    plot_path = os.path.join(OUTPATH, f"{subj}_heatmap.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()

# Save protocol summary
protocol_df = pd.DataFrame(protocol, columns=['subj', 'time', 'status', 'mean_accuracy'])
protocol_path = os.path.join(OUTPATH, 'svm_analysis_protocol.xlsx')
protocol_df.to_excel(protocol_path, index=False)
print(f"Protocol saved to: {protocol_path}")
print("tid_psam_svm_analysis_DONE")
