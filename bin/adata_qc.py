# Ploting for QC
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def plot_qc_density(adata, qc_metric, sample_col="sample_name", cutoff=None, n_cols=4, figsize=(4, 4)):
    
    # Check if sample col exists
    if sample_col not in adata.obs:
        raise ValueError(f"Sample column '{sample_col}' is missing in adata.obs.")

    samples = adata.obs[sample_col].unique()
    n_rows = int(np.ceil(len(samples) / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(figsize[0] * n_cols, figsize[1] * n_rows), constrained_layout=True)

    for ax, sample in zip(axes.flatten(), samples):
        subset = adata.obs[adata.obs[sample_col] == sample]
        sns.kdeplot(subset[qc_metric], ax=ax, fill=True)
        if cutoff is not None:
            ax.axvline(cutoff, color="red", linestyle="--", label=f"Cutoff: {cutoff}")
        ax.set_title(sample)
        ax.set_xlabel(qc_metric)
        ax.set_ylabel("Density")
        ax.legend()

    for ax in axes.flatten()[len(samples):]:  # Hide unused subplots
        ax.set_visible(False)

    plt.show()
    
    
def plot_qc_scatter(adata, x_metric, y_metric, sample_col="sample_name", dot_size=5, x_cutoff=None, y_cutoff=None, n_cols=4, figsize=(4, 4)):
    
    if sample_col not in adata.obs:
        raise ValueError(f"Sample column '{sample_col}' is missing in adata.obs.")

    samples = adata.obs[sample_col].unique()
    n_rows = int(np.ceil(len(samples) / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(figsize[0] * n_cols, figsize[1] * n_rows), constrained_layout=True)

    for ax, sample in zip(axes.flatten(), samples):
        subset = adata.obs[adata.obs[sample_col] == sample]
        sns.scatterplot(x=subset[x_metric], y=subset[y_metric], ax=ax, s=dot_size, alpha=0.5)

        # Add cutoff lines if specified
        if x_cutoff is not None:
            ax.axvline(x_cutoff, color="red", linestyle="--", label=f"X Cutoff: {x_cutoff}")
        if y_cutoff is not None:
            ax.axhline(y_cutoff, color="red", linestyle="--", label=f"Y Cutoff: {y_cutoff}")

        ax.set_title(f"Sample {sample}")
        ax.set_xlabel(x_metric)
        ax.set_ylabel(y_metric)
        # ax.legend()

    for ax in axes.flatten()[len(samples):]:  # Hide unused subplots
        ax.set_visible(False)

    plt.show()
    
def plot_obs_barplot(adata, x_col, fill_col, plot_fractions=False, figsize=(6, 6)):
    
    if x_col not in adata.obs or fill_col not in adata.obs:
        raise ValueError(f"One of the columns '{x_col}' or '{fill_col}' is missing in adata.obs.")
    
    # Add observed=False to silence the warning
    count_df = adata.obs.groupby([x_col, fill_col], observed=False).size().reset_index(name='count')
    
    # Calculate fractions if plot_fractions is True
    if plot_fractions:
        count_df['fraction'] = count_df.groupby(x_col)['count'].transform(lambda x: x / x.sum())
        y_col = 'fraction'
    else:
        y_col = 'count'
    
    # Create the stacked bar plot
    plt.figure(figsize=figsize)
    
    # Get the unique values for x_col
    x_vals = count_df[x_col].unique()
    
    # Create a new dataframe where each column corresponds to a different fill_col category
    stacked_df = count_df.pivot_table(index=x_col, columns=fill_col, values=y_col, aggfunc='sum', fill_value=0)
    
    # Plot the stacked bars
    stacked_df.plot(kind='bar', stacked=True, figsize=figsize)

    plt.xlabel(x_col)
    plt.ylabel("Fraction" if plot_fractions else "Count")
    plt.title(f"Stacked bar plot of {x_col} with {fill_col} ({'fractions' if plot_fractions else 'counts'})")
    plt.legend(title=fill_col)
    plt.show()
    
def plot_obs_violinplot(adata, x_col, y_col, figsize=(6, 6), width=1.0):
    
    if x_col not in adata.obs or y_col not in adata.obs:
        raise ValueError(f"One of the columns '{x_col}' or '{y_col}' is missing in adata.obs.")
    
    plt.figure(figsize=figsize)
    sns.violinplot(data=adata.obs, x=x_col, y=y_col, width=width)
    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.title(f"Violin plot of {y_col} by {x_col}")
    plt.show()