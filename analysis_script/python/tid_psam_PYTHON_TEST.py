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

SCRIPTPATH = os.path.dirname(os.path.abspath(__file__))
expected_subpath = os.path.join('psam', 'analysis_script', 'python')
if not re.search(re.escape(expected_subpath) + r'$', SCRIPTPATH):
    raise EnvironmentError('Path not OK')

MAINPATH = os.path.abspath(os.path.join(SCRIPTPATH, '..', '..'))
INPATH = os.path.join(MAINPATH, 'data', 'processed_data', 'svm_prepared_clean')
OUTPATH = os.path.join(MAINPATH, 'data', 'analysis_data', 'svm_analysis2')
os.makedirs(OUTPATH, exist_ok=True)

dircont_subj_early = sorted([f for f in Path(INPATH).glob("sub-*early.csv")])
dircont_subj_late = sorted([f for f in Path(INPATH).glob("sub-*late.csv")])

C_range = np.logspace(-2, 4, 10)
gamma_range = np.logspace(-3, 1, 10)
n_splits = 5
protocol = []
acc_collection_early, acc_collection_late = [], []

def process_subject(file, condition):
    subj = re.search(r'sub-\d+', file.name).group(0)
    subj_outpath = os.path.join(OUTPATH, subj)
    os.makedirs(subj_outpath, exist_ok=True)
    df = pd.read_csv(file)
    X = df.iloc[:, :-1].values
    y = df.iloc[:, -1].values

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)

    y_encoded = LabelEncoder().fit_transform(y)
    plt.figure(figsize=(6, 5))
    plt.scatter(X_pca[:, 0], X_pca[:, 1], c=y_encoded, cmap='coolwarm', edgecolors='k', alpha=0.7)
    plt.title(f"{subj} - {condition} - PCA Projection")
    plt.xlabel("PC 1")
    plt.ylabel("PC 2")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(subj_outpath, f"{condition}_pca_projection.png"), dpi=300)
    plt.close()

    acc_matrix = np.zeros((len(gamma_range), len(C_range)))
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=123)

    for i, gamma_val in enumerate(gamma_range):
        for j, C_val in enumerate(C_range):
            fold_accuracies = []
            for train_idx, test_idx in skf.split(X, y):
                scaler = StandardScaler()
                X_train_scaled = scaler.fit_transform(X[train_idx])
                X_test_scaled = scaler.transform(X[test_idx])
                pca = PCA(n_components=100)
                X_train_pca = pca.fit_transform(X_train_scaled)
                X_test_pca = pca.transform(X_test_scaled)
                clf = SVC(C=C_val, kernel='rbf', gamma=gamma_val)
                clf.fit(X_train_pca, y[train_idx])
                acc = accuracy_score(y[test_idx], clf.predict(X_test_pca))
                print('fold acc:' + str(acc))
                fold_accuracies.append(acc)
            acc_matrix[i, j] = np.mean(fold_accuracies)

    acc_df = pd.DataFrame(acc_matrix,
                          index=[f"{g:.5f}" for g in gamma_range],
                          columns=[f"{c:.5f}" for c in C_range])
    acc_df.to_csv(os.path.join(subj_outpath, f"{condition}_acc_grid.csv"))

    best_idx = np.unravel_index(np.argmax(acc_matrix), acc_matrix.shape)
    best_gamma = gamma_range[best_idx[0]]
    best_C = C_range[best_idx[1]]

    return acc_matrix, best_idx, subj_outpath

for file_early, file_late in tqdm(zip(dircont_subj_early, dircont_subj_late), total=len(dircont_subj_early), desc="Processing subjects"):
    subj = re.search(r'sub-\d+', file_early.name).group(0)
    tic = time.time()
    acc_early, best_idx_early, outpath = process_subject(file_early, 'early')
    acc_late, best_idx_late, _ = process_subject(file_late, 'late')

    acc_collection_early.append(acc_early)
    acc_collection_late.append(acc_late)
    protocol.append([subj, time.time() - tic, 'OK', np.mean(acc_early), np.mean(acc_late)])

    # Combined individual line plots
    plt.figure(figsize=(10, 5))
    plt.plot(C_range, acc_early[best_idx_early[0], :], label='Early', marker='o')
    plt.plot(C_range, acc_late[best_idx_late[0], :], label='Late', marker='s')
    plt.xscale('log')
    plt.xlabel("C")
    plt.ylabel("Accuracy")
    plt.title(f"{subj} - Accuracy vs C (Best Gamma)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(outpath, f"combined_accuracy_vs_C.png"), dpi=300)
    plt.close()

    plt.figure(figsize=(10, 5))
    plt.plot(gamma_range, acc_early[:, best_idx_early[1]], label='Early', marker='o')
    plt.plot(gamma_range, acc_late[:, best_idx_late[1]], label='Late', marker='s')
    plt.xscale('log')
    plt.xlabel("Gamma")
    plt.ylabel("Accuracy")
    plt.title(f"{subj} - Accuracy vs Gamma (Best C)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(outpath, f"combined_accuracy_vs_Gamma.png"), dpi=300)
    plt.close()

    # Combined individual heatmaps with min/max
    for acc_matrix, label in zip([acc_early, acc_late], ['early', 'late']):
        plt.figure(figsize=(8, 6))
        im = plt.imshow(acc_matrix, origin='lower', aspect='auto',
                        extent=[C_range[0], C_range[-1], gamma_range[0], gamma_range[-1]],
                        cmap='magma', vmin=0.5, vmax=0.95)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel("C")
        plt.ylabel("Gamma")
        plt.title(f"{subj} - {label} - Accuracy Grid")
        cbar = plt.colorbar(im, label='Mean Accuracy')
        min_acc = acc_matrix.min()
        max_acc = acc_matrix.max()
        plt.text(0.05, 0.95, f"Min Acc: {min_acc:.3f}\nMax Acc: {max_acc:.3f}", transform=plt.gca().transAxes,
                 fontsize=10, color='white', ha='left', va='top',
                 bbox=dict(boxstyle='round,pad=0.3', facecolor='black', alpha=0.5))
        plt.tight_layout()
        plt.savefig(os.path.join(outpath, f"{label}_heatmap.png"), dpi=300)
        plt.close()

# Grand average plots
ga_early = np.mean(acc_collection_early, axis=0)
ga_late = np.mean(acc_collection_late, axis=0)

plt.figure(figsize=(16, 6))
plt.subplot(1, 2, 1)
im1 = plt.imshow(ga_early, origin='lower', aspect='auto', extent=[C_range[0], C_range[-1], gamma_range[0], gamma_range[-1]], cmap='magma', vmin=0.5, vmax=0.95)
plt.xscale('log')
plt.yscale('log')
plt.title("Grand Average - Early")
plt.xlabel("C")
plt.ylabel("Gamma")
plt.colorbar(im1, label='Accuracy')
min_acc = ga_early.min()
max_acc = ga_early.max()
plt.text(0.05, 0.95, f"Min Acc: {min_acc:.3f}\nMax Acc: {max_acc:.3f}", transform=plt.gca().transAxes,
         fontsize=10, color='white', ha='left', va='top',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='black', alpha=0.5))

plt.subplot(1, 2, 2)
im2 = plt.imshow(ga_late, origin='lower', aspect='auto', extent=[C_range[0], C_range[-1], gamma_range[0], gamma_range[-1]], cmap='magma', vmin=0.5, vmax=0.95)
plt.xscale('log')
plt.yscale('log')
plt.title("Grand Average - Late")
plt.xlabel("C")
plt.ylabel("Gamma")
plt.colorbar(im2, label='Accuracy')
min_acc = ga_late.min()
max_acc = ga_late.max()
plt.text(0.05, 0.95, f"Min Acc: {min_acc:.3f}\nMax Acc: {max_acc:.3f}", transform=plt.gca().transAxes,
         fontsize=10, color='white', ha='left', va='top',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='black', alpha=0.5))

plt.tight_layout()
plt.savefig(os.path.join(OUTPATH, "grand_average_heatmaps.png"), dpi=300)
plt.close()

best_gamma_idx = np.unravel_index(np.argmax(ga_early), ga_early.shape)[0]
best_C_idx = np.unravel_index(np.argmax(ga_early), ga_early.shape)[1]

plt.figure(figsize=(10, 5))
plt.plot(C_range, ga_early[best_gamma_idx, :], label='Early', marker='o')
plt.plot(C_range, ga_late[best_gamma_idx, :], label='Late', marker='s')
plt.xscale('log')
plt.xlabel("C")
plt.ylabel("Accuracy")
plt.title("GA - Accuracy vs C (Best Gamma)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(OUTPATH, "GA_accuracy_vs_C.png"), dpi=300)
plt.close()

plt.figure(figsize=(10, 5))
plt.plot(gamma_range, ga_early[:, best_C_idx], label='Early', marker='o')
plt.plot(gamma_range, ga_late[:, best_C_idx], label='Late', marker='s')
plt.xscale('log')
plt.xlabel("Gamma")
plt.ylabel("Accuracy")
plt.title("GA - Accuracy vs Gamma (Best C)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(OUTPATH, "GA_accuracy_vs_Gamma.png"), dpi=300)
plt.close()

protocol_df = pd.DataFrame(protocol, columns=['subj', 'time', 'status', 'mean_accuracy_early', 'mean_accuracy_late'])
protocol_df.to_excel(os.path.join(OUTPATH, 'svm_analysis_protocol.xlsx'), index=False)
print("tid_psam_svm_analysis_DONE")
