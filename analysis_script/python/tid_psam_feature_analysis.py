import os
import sys
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import pointbiserialr
import numpy as np
import math

# Define paths
SCRIPTPATH = os.path.dirname(os.path.abspath(__file__))
MAINPATH = os.path.abspath(os.path.join(SCRIPTPATH, '..', '..'))
INPATH = os.path.join(MAINPATH, 'data', 'processed_data', 'svm_prepared_clean')
OUTPATH = os.path.join(MAINPATH, 'data', 'analysis_data', 'feature_analysis')

FUNPATH = os.path.join(MAINPATH, 'functions')
sys.path.append(FUNPATH)

# Load costum functions
from tid_psam_check_folder_clean_up_folder_TD import tid_psam_check_folder_TD, tid_psam_clean_up_folder_TD

tid_psam_check_folder_TD(MAINPATH, INPATH, OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH)

EPOCHS = ['early', 'late']
group_rvals_by_epoch = {epoch: {} for epoch in EPOCHS}

def chunk_list(data, n_chunks):
    if not data:
        return []
    avg = max(1, math.ceil(len(data) / n_chunks))
    return [data[i:i + avg] for i in range(0, len(data), avg)]

def plot_horizontal_bar_columns(title, data_dict, color_dict, outpath, n_cols=10, format='png'):
    items = sorted(data_dict.items(), key=lambda x: x[1])
    chunks = chunk_list(items, n_cols)

    fig_width = 4.5 * n_cols
    fig_height = max(4, len(chunks[0]) * 0.4) if chunks else 4

    fig, axs = plt.subplots(1, len(chunks), figsize=(fig_width, fig_height), squeeze=False)

    for ax, chunk in zip(axs[0], chunks):
        names = [k for k, _ in chunk]
        values = [v for _, v in chunk]
        colors = [color_dict.get(k, 'blue') for k in names]

        y_pos = np.arange(len(names))
        bars = ax.barh(y_pos, values, color=colors)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(names, fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel('Effect Size (r)')
        ax.set_xlim(min(values) - 0.1, max(values) + 0.1)

        for i, (bar, val) in enumerate(zip(bars, values)):
            ax.text(val + 0.02 * np.sign(val), i, f'{val:.2f}', va='center', fontsize=6)

    fig.suptitle(title, fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(outpath, dpi=300)
    plt.close()

# === MAIN PROCESSING LOOP ===
for epoch in EPOCHS:
    print(f"\n=== Processing {epoch.upper()} epoch ===")
    filelist = list(Path(INPATH).glob(f"sub-*{epoch}.csv"))

    for file in filelist:
        subj_name = re.search(r'sub-\d+', file.name).group(0)
        print(f"\nProcessing {subj_name}...")

        df = pd.read_csv(file)

        cond_col = 'labels'
        original_labels = df[cond_col].unique()
        label_map = {label: idx for idx, label in enumerate(sorted(original_labels))}
        reverse_label_map = {v: k for k, v in label_map.items()}
        df['labels_numeric'] = df[cond_col].map(label_map)
        print(f"Label mapping: {label_map}")

        feature_cols = [col for col in df.columns if col not in ['labels', 'trial', 'labels_numeric']]
        df['trial'] = df.groupby('labels_numeric').cumcount()

        subj_rvals = {}
        subj_signif = {}

        # Output path format: OUTPATH/sub-X/epoch
        subj_outpath = os.path.join(OUTPATH, subj_name, epoch)
        os.makedirs(subj_outpath, exist_ok=True)

        for val_col in feature_cols:
            x = pd.to_numeric(df[val_col], errors='coerce').values
            y = pd.to_numeric(df['labels_numeric'], errors='coerce').values

            valid_mask = ~np.isnan(x) & ~np.isnan(y)
            x_valid = x[valid_mask]
            y_valid = y[valid_mask]

            if len(x_valid) < 2 or len(np.unique(y_valid)) < 2:
                print(f"Skipping feature '{val_col}': not enough valid or variable data.")
                continue

            try:
                r, p_val = pointbiserialr(y_valid, x_valid)
            except Exception as e:
                print(f"Skipping feature '{val_col}': {e}")
                continue

            subj_rvals[val_col] = r
            subj_signif[val_col] = p_val < 0.05

            print(f"Feature '{val_col}': r = {r:.4f}, p = {p_val:.4g}")

            if p_val < 0.10:
                plt.figure(figsize=(8, 6))
                temp_df = df[[cond_col, val_col]].copy()
                temp_df[cond_col] = temp_df[cond_col].astype(str)

                sns.boxplot(data=temp_df, x=cond_col, y=val_col, whis=[5, 95], width=0.5, fliersize=0)
                sns.stripplot(data=temp_df, x=cond_col, y=val_col, color='black', alpha=0.6, jitter=True)
                plt.title(f"{subj_name} - {val_col} ({epoch})\nr = {r:.3f}, p = {p_val:.3f}")
                plt.ylabel(val_col)
                plt.xlabel("Condition")
                plt.grid(True, axis='y')
                plt.tight_layout()

                # Save as sub-XX_feature_epoch.png
                plot_path = os.path.join(subj_outpath, f"{subj_name}_{val_col}_{epoch}.png")
                plt.savefig(plot_path, dpi=300)
                plt.close()

        # Store Fisher-Z transformed r-values
        for feat, r in subj_rvals.items():
            if np.abs(r) < 1:
                z = np.arctanh(r)
                group_rvals_by_epoch[epoch].setdefault(feat, []).append(z)

        # Per-subject bar plot
        color_dict = {feat: 'red' if subj_signif[feat] else 'blue' for feat in subj_rvals}
        barplot_path = os.path.join(subj_outpath, f"{subj_name}_{epoch}_rvals_barplot.png")
        plot_horizontal_bar_columns(
            title=f"{subj_name} - Point-Biserial Correlation by Feature ({epoch})",
            data_dict=subj_rvals,
            color_dict=color_dict,
            outpath=barplot_path,
            n_cols=10,
            format='png'
        )

# === GROUP PLOTS BY EPOCH ===
for epoch, zvals_dict in group_rvals_by_epoch.items():
    mean_group_rvals = {}
    for feat, zvals in zvals_dict.items():
        mean_z = np.mean(zvals)
        mean_r = np.tanh(mean_z)
        mean_group_rvals[feat] = mean_r

    color_dict_group = {feat: 'blue' for feat in mean_group_rvals}
    group_barplot_path = os.path.join(OUTPATH, f"group_mean_rvals_barplot_{epoch}.png")
    plot_horizontal_bar_columns(
        title=f"Group-Averaged r-values (Fisher Z-transformed) ({epoch})",
        data_dict=mean_group_rvals,
        color_dict=color_dict_group,
        outpath=group_barplot_path,
        n_cols=10,
        format='png'
    )
