import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import ttest_rel
import numpy as np
import math

# Define paths
SCRIPTPATH = os.path.dirname(os.path.abspath(__file__))
MAINPATH = os.path.abspath(os.path.join(SCRIPTPATH, '..', '..'))
INPATH = os.path.join(MAINPATH, 'data', 'processed_data', 'svm_prepared_clean')
OUTPATH = os.path.join(MAINPATH, 'data', 'analysis_data', 'svm_analysis_TEST')
os.makedirs(OUTPATH, exist_ok=True)

# Define both epochs
EPOCHS = ['early', 'late']

# Store t-values for group-level plot (combined across epochs)
group_tvals_by_epoch = {epoch: {} for epoch in EPOCHS}

def chunk_list(data, n_chunks):
    avg = math.ceil(len(data) / n_chunks)
    return [data[i:i + avg] for i in range(0, len(data), avg)]

def plot_horizontal_bar_columns(title, data_dict, color_dict, outpath, n_cols=10, format='png'):
    items = sorted(data_dict.items(), key=lambda x: x[1])
    chunks = chunk_list(items, n_cols)

    fig_width = 4.5 * n_cols
    fig_height = max(4, len(chunks[0]) * 0.4)

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
        ax.set_xlabel('t-value')
        ax.set_xlim(min(values) - 0.5, max(values) + 0.5)

        for i, (bar, val) in enumerate(zip(bars, values)):
            ax.text(val + 0.05 * np.sign(val), i, f'{val:.2f}', va='center', fontsize=6)

    fig.suptitle(title, fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(outpath, dpi=600, format=format)
    plt.close()

# === MAIN PROCESSING LOOP ===
for epoch in EPOCHS:
    print(f"\n=== Processing {epoch.upper()} epoch ===")
    filelist = list(Path(INPATH).glob(f"sub-*{epoch}.csv"))

    for file in filelist:
        subj_name = re.search(r'sub-\d+', file.name).group(0)
        print(f"\nProcessing {subj_name}...")

        subj_outpath = os.path.join(OUTPATH, epoch, subj_name)
        os.makedirs(subj_outpath, exist_ok=True)

        df = pd.read_csv(file)

        cond_col = 'labels'
        feature_cols = [col for col in df.columns if col not in ['labels', 'trial']]
        df['trial'] = df.groupby(cond_col).cumcount()

        subj_tvals = {}
        subj_signif = {}

        for val_col in feature_cols:
            df_wide = df.pivot(index='trial', columns=cond_col, values=val_col)

            if df_wide.shape[1] != 2:
                print(f"Warning: Expected 2 conditions for '{val_col}', found {df_wide.shape[1]}. Skipping.")
                continue

            conds = df_wide.columns.tolist()
            x = df_wide[conds[0]].values
            y = df_wide[conds[1]].values

            t_stat, p_val = ttest_rel(x, y, nan_policy='omit')
            subj_tvals[val_col] = t_stat
            subj_signif[val_col] = p_val < 0.05

            print(f"Feature '{val_col}': t = {t_stat:.4f}, p = {p_val:.4g}")

            if p_val < 0.10:
                plt.figure(figsize=(8, 6))
                sns.boxplot(data=df, x=cond_col, y=val_col, whis=[5, 95], width=0.5, fliersize=0)
                sns.stripplot(data=df, x=cond_col, y=val_col, color='black', alpha=0.6, jitter=True)
                plt.title(f"{subj_name} - {val_col} ({epoch}, p = {p_val:.3f})")
                plt.ylabel(val_col)
                plt.xlabel("Condition")
                plt.grid(True, axis='y')
                plt.tight_layout()

                plot_filename = f"{val_col}_p{p_val:.3f}.png".replace('.', '_')
                plot_path = os.path.join(subj_outpath, plot_filename)
                plt.savefig(plot_path, dpi=600, format='png')
                plt.close()

        # Store for group-level analysis
        for feat, tval in subj_tvals.items():
            group_tvals_by_epoch[epoch].setdefault(feat, []).append(tval)

        # Per-subject bar plot
        color_dict = {feat: 'red' if subj_signif[feat] else 'blue' for feat in subj_tvals}
        barplot_path = os.path.join(subj_outpath, "tvals_barplot.png")
        plot_horizontal_bar_columns(
            title=f"{subj_name} - Paired t-test by Feature ({epoch})",
            data_dict=subj_tvals,
            color_dict=color_dict,
            outpath=barplot_path,
            n_cols=10,
            format='png'
        )

# === GROUP PLOTS BY EPOCH ===
for epoch, tvals_dict in group_tvals_by_epoch.items():
    mean_group_tvals = {feat: np.mean(tvals) for feat, tvals in tvals_dict.items()}
    color_dict_group = {feat: 'blue' for feat in mean_group_tvals}
    group_barplot_path = os.path.join(OUTPATH, f"group_mean_tvals_barplot_{epoch}.pdf")
    plot_horizontal_bar_columns(
        title=f"Group-Averaged t-values ({epoch})",
        data_dict=mean_group_tvals,
        color_dict=color_dict_group,
        outpath=group_barplot_path,
        n_cols=10,
        format='pdf'
    )
