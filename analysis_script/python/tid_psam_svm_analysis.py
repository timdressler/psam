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
from sklearn.decomposition import KernelPCA
from sklearn.decomposition import PCA
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

# Color definitions
def hex_to_rgb(color):
    color = color.lstrip('#')
    return tuple(int(color[i:i+2], 16) / 255 for i in (0, 2, 4))

main_blue = hex_to_rgb('#004F9F')
main_red = hex_to_rgb('#D53D0E')
main_green = hex_to_rgb('#00786B')

# Subject files
dircont_subj_early = [f for f in Path(INPATH).glob("sub-95*early.csv")]

# Protocol
protocol = []

# Hyperparameter ranges (log-spaced)
lower_C = -5
upper_C = 5
lower_gamma = 0
upper_gamma = 100
C_range = np.logspace(lower_C, upper_C, 32)
#C_range = np.linspace(lower_C, upper_C, 32)
gamma_range = np.logspace(lower_gamma, upper_gamma, 32)

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

    ######
    df2 = pd.read_csv(file)
    # Assume the last column is the label
    label_col = df2.columns[-1]
    # Group by the label and compute the mean of each feature
    means_by_class = df2.groupby(label_col).mean()
    print(means_by_class)
    sd_by_class = df2.groupby(label_col).std()
    print(sd_by_class)
    ######

    acc_matrix = np.zeros((len(gamma_range), len(C_range)))
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=123)

    print('start grid search...')
    for i, gamma_val in enumerate(gamma_range):
        for j, C_val in enumerate(C_range):
            fold_accuracies = []

            for train_idx, test_idx in skf.split(X, y):
                X_train, X_test = X[train_idx], X[test_idx]
                y_train, y_test = y[train_idx], y[test_idx]

                kpca = KernelPCA(n_components=50,
                                 kernel='poly')
                print('pca')
                X_train_kpca = kpca.fit_transform(X_train)
                X_test_kpca = kpca.transform(X_test)

                # Use standard PCA
                #kpca = PCA(n_components=120)
                #X_train_kpca = kpca.fit_transform(X_train)
                #X_test_kpca = kpca.transform(X_test)

                #explained_var = kpca.explained_variance_ratio_.sum()
                #print(f"Total variance retained: {explained_var:.2f}")

                print("svm")
                clf = SVC(C=C_val, kernel='rbf', gamma=gamma_val) # class_weight = 'balanced'
                clf.fit(X_train, y_train)
                y_pred = clf.predict(X_test)
                #print(y_pred)
                acc = accuracy_score(y_test, y_pred)
                fold_accuracies.append(acc)

            acc_matrix[i, j] = np.mean(fold_accuracies)
            print(str(np.mean(fold_accuracies)))

    subj_time = time.time() - tic
    protocol.append([subj, subj_time, 'OK', np.mean(acc_matrix)])

    # Save accuracy grid
    acc_df = pd.DataFrame(acc_matrix,
                          index=[f"{g:.5f}" for g in gamma_range],
                          columns=[f"{c:.5f}" for c in C_range])
    acc_path = os.path.join(OUTPATH, f"{subj}_acc_grid.csv")
    acc_df.to_csv(acc_path)

    # Plot heatmap
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
    plt.tight_layout()

    plot_path = os.path.join(OUTPATH, f"{subj}_heatmap.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()

# Save protocol
protocol_df = pd.DataFrame(protocol, columns=['subj', 'time', 'status', 'mean_accuracy'])
protocol_path = os.path.join(OUTPATH, 'svm_analysis_protocol.xlsx')
protocol_df.to_excel(protocol_path, index=False)
print(f"Protocol saved to: {protocol_path}")
print("tid_psam_svm_analysis_DONE")
