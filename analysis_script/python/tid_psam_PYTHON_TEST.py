import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import ttest_rel

# Define paths
SCRIPTPATH = os.path.dirname(os.path.abspath(__file__))
MAINPATH = os.path.abspath(os.path.join(SCRIPTPATH, '..', '..'))
INPATH = os.path.join(MAINPATH, 'data', 'processed_data', 'svm_prepared_clean')

files = list(Path(INPATH).glob("sub-96*late.csv"))

for file in files:
    subj_name = re.search(r'sub-\d+', file.name).group(0)
    print(f"\nProcessing {subj_name}...")

    df = pd.read_csv(file)

    print(df.head())

    # Replace with your actual column names
    cond_col = 'labels'
    val_col = 'mean_E01'

    # Create trial index within each condition
    df['trial'] = df.groupby(cond_col).cumcount()

    # Pivot for stats
    df_wide = df.pivot(index='trial', columns=cond_col, values=val_col)

    print("Pivoted DataFrame head:")
    print(df_wide.head())

    if df_wide.shape[1] != 2:
        print(f"Warning: Expected 2 conditions, found {df_wide.shape[1]}. Skipping {subj_name}.")
        continue

    conds = df_wide.columns.tolist()
    x = df_wide[conds[0]].values
    y = df_wide[conds[1]].values

    t_stat, p_val = ttest_rel(x, y, nan_policy='omit')
    print(f"Paired t-test between '{conds[0]}' and '{conds[1]}':")
    print(f"  t-statistic = {t_stat:.4f}, p-value = {p_val:.4g}")

    # Boxplot with individual points
    plt.figure(figsize=(8, 6))
    sns.boxplot(data=df, x=cond_col, y=val_col, whis=[5, 95], width=0.5, fliersize=0)
    sns.stripplot(data=df, x=cond_col, y=val_col, color='black', alpha=0.6, jitter=True)
    plt.title(f"{subj_name} - {val_col} by Condition")
    plt.ylabel(val_col)
    plt.xlabel("Condition")
    plt.grid(True, axis='y')
    plt.tight_layout()
    plt.show()
