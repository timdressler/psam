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
from sklearn.decomposition import PCA
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
dircont_subj_early = [f for f in Path(INPATH).glob("sub-96*late.csv")]

# Protocol
protocol = []

# Hyperparameter ranges (narrowed)
C_range = np.logspace(-2, 4, 10)
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

    # === Manual Feature Plot: mean_E01 vs mean_E02 ===
    if 'mean_E01' in df.columns and 'mean_E02' in df.columns:
        x_manual = df['mean_E01'].values
        y_manual = df['mean_E02'].values

        plt.figure(figsize=(6, 5))
        plt.scatter(x_manual, y_manual, c=y_encoded, cmap='coolwarm', edgecolors='k', alpha=0.7)
        plt.xlabel("mean_E01")
        plt.ylabel("mean_E02")
        plt.title(f"{subj} - Feature Space: mean_E01 vs mean_E02", fontsize=12)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(OUTPATH, f"{subj}_manual_features_plot.png"), dpi=300)
        plt.close()
    else:
        print(f"Warning: One or both of mean_E01 and mean_E02 not found in {subj}!")

    acc_matrix = np.zeros((len(gamma_range), len(C_range)))
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=123)

    print('start grid search...')
    for i, gamma_val in enumerate(gamma_range):
        for j, C_val in enumerate(C_range):
            fold_accuracies = []

            for fold_idx, (train_idx, test_idx) in enumerate(skf.split(X, y)):
                X_train, X_test = X[train_idx], X[test_idx]
                y_train, y_test = y[train_idx], y[test_idx]

                scaler = StandardScaler()
                X_train_scaled = scaler.fit_transform(X_train)
                X_test_scaled = scaler.transform(X_test)

                pca = PCA(n_components=100)
                X_train_pca = pca.fit_transform(X_train_scaled)
                X_test_pca = pca.transform(X_test_scaled)

                explained_var = pca.explained_variance_ratio_.sum()
                print(f"Total variance retained: {explained_var:.2f}")

                clf = SVC(C=C_val, kernel='rbf', gamma=gamma_val)
                clf.fit(X_train_pca, y_train)
                y_pred = clf.predict(X_test_pca)

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

    # === NEW: Plot accuracy vs C and Gamma for best hyperparameters ===
    best_idx = np.unravel_index(np.argmax(acc_matrix), acc_matrix.shape)
    best_gamma = gamma_range[best_idx[0]]
    best_C = C_range[best_idx[1]]
    print(f"Best gamma: {best_gamma}, Best C: {best_C}")

    # Plot accuracy vs C
    plt.figure(figsize=(8, 5))
    plt.plot(C_range, acc_matrix[best_idx[0], :], marker='o', label=f'Gamma = {best_gamma:.4f}')
    plt.fill_between(C_range,
                     acc_matrix[best_idx[0], :] - 0.01,
                     acc_matrix[best_idx[0], :] + 0.01,
                     alpha=0.2)
    plt.xscale('log')
    plt.ylim(0.4, 1.05)
    plt.xlabel("C")
    plt.ylabel("Accuracy")
    plt.title(f"{subj} - Accuracy vs C (Gamma = {best_gamma:.4f})")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPATH, f"{subj}_accuracy_vs_C.png"), dpi=300)
    plt.close()

    # Plot accuracy vs Gamma
    plt.figure(figsize=(8, 5))
    plt.plot(gamma_range, acc_matrix[:, best_idx[1]], marker='o', label=f'C = {best_C:.4f}')
    plt.fill_between(gamma_range,
                     acc_matrix[:, best_idx[1]] - 0.01,
                     acc_matrix[:, best_idx[1]] + 0.01,
                     alpha=0.2)
    plt.xscale('log')
    plt.ylim(0.4, 1.05)
    plt.xlabel("Gamma")
    plt.ylabel("Accuracy")
    plt.title(f"{subj} - Accuracy vs Gamma (C = {best_C:.4f})")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPATH, f"{subj}_accuracy_vs_Gamma.png"), dpi=300)
    plt.close()
    # === END NEW PLOTS ===

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
