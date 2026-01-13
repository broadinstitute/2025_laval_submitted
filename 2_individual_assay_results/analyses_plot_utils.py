import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from itertools import combinations
from sklearn.metrics import precision_recall_curve, average_precision_score, roc_auc_score, roc_curve
from scipy.stats import fisher_exact, mannwhitneyu
from statannotations.Annotator import Annotator

set2_palette = sns.color_palette("Set2")


def plot_cat_count_perc(df, cat, title="", ax=None, palette="Dark2"):
    clin_counts = df[cat].value_counts().sort(cat)
    clin_counts = clin_counts.with_columns(
        (pl.col("count") / pl.col("count").sum() * 100).alias("percent")
    )
    clin_counts = clin_counts.to_pandas()
    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(6, clin_counts.shape[0]))
    sns.barplot(data=clin_counts,
               y=cat,
               x="count",
               hue=cat,
               palette=palette,
               ax=ax
    )
    # percentage annotations
    max_tot = clin_counts["count"].max()
    for idx, row in clin_counts.iterrows():
        ax.text(
            max_tot * 0.01,
            idx,
            f"{row['percent']:.1f}%",
            va="center",
            ha="left",
            fontweight="bold",
            fontsize=12,
            color="black"
        )
        ax.grid(alpha=.2, axis="x")
    return clin_counts
    

def plot_assay_hit_by_category(df, hit_col, cat_cols, title="", ax=None):
    """
    plot_assay_hit_by_category
    ==========================
    
    Create a paired horizontal bar chart that shows, for every category in
    `cat_cols`, both the **total number of observations** and the **number of
    positive hits** (as defined by `hit_col`).  The chart also annotates the
    hit-rate as a percentage.
    
    Parameters
    ----------
    df : polars.DataFrame
        Input table containing at least the columns listed in `cat_cols` and a
        boolean column indicating hits.
    hit_col : str
        Name of the boolean column that flags hits (e.g. `"is_hit"`).
    cat_cols : str or list[str]
        Column(s) that define the categorical grouping(s).  If a single string is
        supplied it is handled automatically.
    title : str, optional
        Title placed above the plot.
    ax : matplotlib.axes.Axes, optional
        Pre-existing axes to draw on.  If None, a new figure and axes are created.
    
    Returns
    -------
    polars.DataFrame
        A tidy frame with one row per category, containing:
            - the grouping key(s)
            - total_counts   : total observations in that category
            - hit_counts     : observations that were hits
            - hit_perc       : hit rate in percent (hit_counts / total_counts * 100)
    
    Example
    -------
    >>> plot_assay_hit_by_category(
    ...     df=my_data,
    ...     hit_col="is_hit",
    ...     cat_cols="clinvar_clnsig_clean",
    ...     title="Hit rate by ClinVar significance"
    ... )
    """
    # --- compute hit counts and percentages ---
    df_tot_cnt_per_cat = (
        df.group_by(cat_cols)
        .len()
        .rename({"len": "total_counts"})
    )
    df_hit_cnt_per_cat = (
        df.unique("gene_variant")
        .filter(pl.col(hit_col))
        .group_by(cat_cols)
        .len()
        .rename({"len": "hit_counts"})
    )

    img_hits_perc_df = (
        df_tot_cnt_per_cat.join(
            df_hit_cnt_per_cat,
            on=cat_cols,
            how="left"
        )
        .fill_null(0)
        .with_columns(
            (pl.col("hit_counts") / pl.col("total_counts") * 100).alias("hit_perc")
        )
        .sort(cat_cols)
    )

    # --- plotting ---
    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(6, img_hits_perc_df.shape[0]))

    plot_df = img_hits_perc_df.to_pandas()

    # total bars (blue)
    sns.barplot(
        data=plot_df,
        y=cat_cols if isinstance(cat_cols, str) else cat_cols[-1],
        x="total_counts",
        color="#1f77b4",
        label="Total count",
        ax=ax
    )

    # hit bars (red)
    sns.barplot(
        data=plot_df,
        y=cat_cols if isinstance(cat_cols, str) else cat_cols[-1],
        x="hit_counts",
        color="#d62728",
        label="Hit count",
        ax=ax
    )

    ax.set_title(title)
    ax.grid(alpha=.2)

    # percentage annotations
    max_tot = plot_df["total_counts"].max()
    for idx, row in plot_df.iterrows():
        ax.text(
            row["hit_counts"] + max_tot * 0.01,
            idx,
            f"{row['hit_perc']:.1f}%",
            va="center",
            ha="left",
            fontweight="bold",
            fontsize=12,
            color="black"
        )

    return img_hits_perc_df


def plot_assay_hit_by_category_perc(df, hit_col, cat_cols, title="", ax=None):
    """
    plot_assay_hit_by_category
    ==========================
    
    Create a stacked horizontal bar chart showing hit rate percentage per category.
    Each bar represents 100% of the category, with the hit portion highlighted.
    
    Parameters
    ----------
    df : polars.DataFrame
        Input table containing at least the columns listed in `cat_cols` and a
        boolean column indicating hits.
    hit_col : str
        Name of the boolean column that flags hits (e.g. `"is_hit"`).
    cat_cols : str or list[str]
        Column(s) that define the categorical grouping(s).  If a single string is
        supplied it is handled automatically.
    title : str, optional
        Title placed above the plot.
    ax : matplotlib.axes.Axes, optional
        Pre-existing axes to draw on.  If None, a new figure and axes are created.
    
    Returns
    -------
    polars.DataFrame
        A tidy frame with one row per category, containing:
            - the grouping key(s)
            - total_counts   : total observations in that category
            - hit_counts     : observations that were hits
            - hit_perc       : hit rate in percent (hit_counts / total_counts * 100)
    """
    # --- compute hit counts and percentages ---
    df_tot_cnt_per_cat = (
        df.group_by(cat_cols)
        .len()
        .rename({"len": "total_counts"})
    )
    df_hit_cnt_per_cat = (
        df.unique("gene_variant")
        .filter(pl.col(hit_col))
        .group_by(cat_cols)
        .len()
        .rename({"len": "hit_counts"})
    )

    cat_col_name = cat_cols if isinstance(cat_cols, str) else cat_cols[-1]
    
    img_hits_perc_df = (
        df_tot_cnt_per_cat.join(
            df_hit_cnt_per_cat,
            on=cat_cols,
            how="left"
        )
        .fill_null(0)
        .with_columns(
            (pl.col("hit_counts") / pl.col("total_counts") * 100).alias("hit_perc")
        )
        .sort(cat_col_name, descending=True)
    )

    # --- plotting ---
    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(5, img_hits_perc_df.shape[0]))

    plot_df = img_hits_perc_df.to_pandas()
    y_positions = range(len(plot_df))

    # background bars (full 100% - gray)
    ax.barh(y_positions, [100] * len(plot_df), color='#00A3A3', label="Non-hit")
    # hit bars (foreground - red)
    ax.barh(y_positions, plot_df["hit_perc"], color='#E67300', label="Hit")

    # set labels
    ax.set_yticks(y_positions)
    ax.set_yticklabels([f"{cat}\n(n={n})" for cat, n in zip(plot_df[cat_col_name], plot_df["total_counts"])]) ##.split('_')[1]
    ax.set_xlabel("Percentage (%)")
    ax.set_xlim(0, 100)
    ax.set_title(title)
    ax.grid(axis='x', alpha=.2)
    ax.legend(loc='lower right')

    # percentage annotations - positioned at the end of hit bar
    for idx, row in plot_df.iterrows():
        ax.text(
            row["hit_perc"] + 1,
            idx,
            f"{row['hit_perc']:.1f}%",
            va="center",
            ha="left",
            fontweight="bold",
            fontsize=11,
            color="white"
        )
    return img_hits_perc_df


def plot_assay_scores(df, assay, category, ax=None, type="box", pairs=None, order=None, palette=None):
    """
    plot_assay_scores
    =================
    
    Create a box/boxen/violin plot of assay scores by category with statistical
    annotations using Mann-Whitney U test.
    
    Parameters
    ----------
    df : polars.DataFrame
        Input table containing at least the columns `assay` and `category`.
    assay : str
        Name of the column containing assay scores.
    category : str
        Name of the categorical column to group by.
    ax : matplotlib.axes.Axes, optional
        Pre-existing axes to draw on. If None, a new figure and axes are created.
    type : str, optional
        Type of plot: "box", "boxen", or "violin". Default is "box".
    pairs : list of tuples, optional
        List of category pairs to compare (e.g., [("PLP", "BLB")]).
        If None, will automatically create pairs from unique categories.
        If exactly 2 categories exist, creates one pair. If more, compares all pairs.
    order : list, optional
        Order of categories on x-axis. If None, uses sorted unique categories.
    palette : dict or list, optional
        Color palette for the plot.
    
    Returns
    -------
    matplotlib.axes.Axes
        The axes object with the plot.
    """
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(3, 3.5))
    
    # Filter out null values in assay column
    df_filtered = df.filter(pl.col(assay).is_not_null())
    
    # Get unique categories from filtered data - preserve original order
    unique_cats = df_filtered[category].unique().to_list()
    
    if len(unique_cats) == 0:
        raise ValueError(f"No valid data found for {assay} and {category}")
    
    # Determine order - preserve original order if not specified
    if order is None:
        order = unique_cats
    else:
        # Filter order to only include categories that exist in the data
        order = [cat for cat in order if cat in unique_cats]
        # Add any missing categories from unique_cats that aren't in order
        for cat in unique_cats:
            if cat not in order:
                order.append(cat)
    
    # Determine pairs for statistical testing
    if pairs is None:
        if len(unique_cats) == 2:
            # If exactly 2 categories, create one pair
            pairs = [(unique_cats[0], unique_cats[1])]
        elif len(unique_cats) > 2:
            # If more than 2, compare all pairs
            pairs = list(combinations(unique_cats, 2))
        else:
            # Only one category, no pairs to compare
            pairs = []
    
    # Convert to pandas for plotting
    plot_df = df_filtered.to_pandas()
    
    # Create the plot
    if type == "box":
        sns.boxplot(
            data=plot_df,
            x=category,
            y=assay,
            palette=palette,
            width=0.4,
            notch=True,
            order=order,
            ax=ax,
            legend=False
        )
    elif type == "boxen":
        sns.boxenplot(
            data=plot_df,
            x=category,
            y=assay,
            palette=palette,
            width=0.6,
            order=order,
            ax=ax,
            legend=False
        )
    elif type == "violin":
        sns.violinplot(
            data=plot_df,
            x=category,
            y=assay,
            palette=palette,
            cut=0,
            order=order,
            ax=ax,
            legend=False
        )

    plot_count = df_filtered.group_by(category).agg(pl.len())
    plot_count_dict = dict(zip(plot_count[category], plot_count["len"]))

    x_labels = [
        f"{cat}\n(n={n})" 
        for cat, n in plot_count_dict.items()
    ]
    ax.set_xticks(range(len(plot_count)))
    ax.set_xticklabels(x_labels, fontsize=11)
    
    # Add statistical annotations if pairs exist
    if len(pairs) > 0:
        # Custom p-value formatter function
        def format_pvalue(pvalue):
            """Format p-value as 'p = xxx' or 'p < 0.0001' if below threshold"""
            if pvalue < 0.0001:
                return "p < 0.0001"
            else:
                return f"p = {pvalue:.4f}"
        
        # Compute p-values manually for custom formatting
        pvalues_dict = {}
        for pair in pairs:
            cat1, cat2 = pair
            data1 = plot_df[plot_df[category] == cat1][assay].dropna()
            data2 = plot_df[plot_df[category] == cat2][assay].dropna()
            
            if len(data1) > 0 and len(data2) > 0:
                try:
                    _, pval = mannwhitneyu(data1, data2, alternative='two-sided')
                    pvalues_dict[pair] = pval
                except Exception:
                    pvalues_dict[pair] = np.nan
            else:
                pvalues_dict[pair] = np.nan
        
        annotator = Annotator(
            ax, 
            pairs,
            data=plot_df,
            x=category, 
            y=assay,
            order=order,
            verbose=False
        )
        annotator.configure(
            test="Mann-Whitney",
            text_format='simple', 
            loc='inside', 
            text_offset=4,
            verbose=False
        )
        annotator.apply_test()
        annotator.annotate(line_offset_to_group=0.05, line_offset=0.15)
        
        # Add custom formatted p-value text annotations
        # Get the y-axis limits to position text
        y_max = ax.get_ylim()[1]
        
        # for pair in pairs:
        #     if pair in pvalues_dict and not np.isnan(pvalues_dict[pair]):
        #         cat1, cat2 = pair
        #         pval = pvalues_dict[pair]
        #         formatted_text = format_pvalue(pval)
                
        #         # Get x positions of the categories
        #         if cat1 in order and cat2 in order:
        #             x1 = order.index(cat1)
        #             x2 = order.index(cat2)
        #             x_pos = (x1 + x2) / 2
                    
        #             # Get max y value for positioning
        #             data1 = plot_df[plot_df[category] == cat1][assay].dropna()
        #             data2 = plot_df[plot_df[category] == cat2][assay].dropna()
        #             if len(data1) > 0 and len(data2) > 0:
        #                 max_y = max(data1.max(), data2.max())
        #                 # Position text above the comparison line
        #                 text_y = max_y + (y_max - max_y) * 0.1
                        
        #                 # Add text annotation
        #                 ax.text(x_pos, text_y, formatted_text,
        #                        ha='center', va='bottom', fontsize=9,
        #                        fontweight='bold', color='black')

    ax.set_xlabel("")
    ax.set_xlim(-.5, 1.5)
    plt.tight_layout()
    plt.show()
    
    return ax


def plot_assay_hit_fisher_test(df, hit_col, cat_col, alternative='two-sided', plot=False, ax=None, title="", palette=None):
    """
    plot_assay_hit_fisher_test
    ==========================
    
    Perform Fisher's exact test on a contingency table comparing assay hits
    across categories. Creates a contingency table and computes the odds ratio
    and p-value using Fisher's exact test. Optionally creates a bar plot showing
    hit percentages with sample counts.
    
    Parameters
    ----------
    df : polars.DataFrame
        Input table containing at least the columns `hit_col` and `cat_col`.
    hit_col : str
        Name of the boolean column that flags hits (e.g. `"dualipa_hit"`).
    cat_col : str
        Name of the categorical column to compare (e.g. `"clinvar_clnsig_clean"`).
        Must have exactly 2 unique categories for Fisher's exact test.
    alternative : str, optional
        Defines the alternative hypothesis. Options are:
        - 'two-sided' (default): two-tailed test
        - 'less': one-tailed test, odds ratio < 1
        - 'greater': one-tailed test, odds ratio > 1
    plot : bool, optional
        If True, create a bar plot showing hit percentages. Default is False.
    ax : matplotlib.axes.Axes, optional
        Pre-existing axes to draw on. If None and plot=True, a new figure and axes are created.
    title : str, optional
        Title for the plot. If empty and plot=True, a default title will be generated.
    palette : dict or list, optional
        Color palette for the plot. If None, uses default colors.
    
    Returns
    -------
    dict
        Dictionary containing:
            - 'contingency_table': polars.DataFrame
                The contingency table with hit/non-hit counts by category
            - 'odds_ratio': float
                The odds ratio from Fisher's exact test
            - 'p_value': float
                The p-value from Fisher's exact test
            - 'hit_rates': polars.DataFrame
                Hit rates (percentages) for each category
            - 'table_array': numpy.ndarray
                The 2x2 numpy array used for Fisher's exact test in format:
                [[cat1_hits, cat1_non_hits], [cat2_hits, cat2_non_hits]]
    
    Example
    -------
    >>> result = plot_assay_hit_fisher_test(
    ...     df=my_data,
    ...     hit_col="dualipa_hit",
    ...     cat_col="clinvar_clnsig_clean",
    ...     plot=True
    ... )
    >>> print(f"Odds ratio: {result['odds_ratio']:.2f}, p-value: {result['p_value']:.2e}")
    """
    # Filter out null values in both columns
    df_filtered = df.filter(
        pl.col(hit_col).is_not_null(),
        pl.col(cat_col).is_not_null()
    )
    
    # Get unique categories - must be exactly 2 for Fisher's exact test
    categories = sorted(df_filtered[cat_col].unique().to_list())
    
    if len(categories) != 2:
        raise ValueError(f"Fisher's exact test requires exactly 2 categories, but found {len(categories)}: {categories}")
    
    cat1, cat2 = categories
    
    # Calculate counts for each category
    # Total counts per category
    df_tot_cnt_per_cat = (
        df_filtered
        .group_by(cat_col)
        .len()
        .rename({"len": "total_counts"})
    )
    
    # Hit counts per category
    df_hit_cnt_per_cat = (
        df_filtered
        .filter(pl.col(hit_col))
        .group_by(cat_col)
        .len()
        .rename({"len": "hit_counts"})
    )
    
    # Join to get hit and total counts
    counts_df = (
        df_tot_cnt_per_cat
        .join(df_hit_cnt_per_cat, on=cat_col, how="left")
        .fill_null(0)
        .with_columns(
            (pl.col("total_counts") - pl.col("hit_counts")).alias("non_hit_counts")
        )
    )
    
    # Extract counts for each category
    cat1_row = counts_df.filter(pl.col(cat_col) == cat1)
    cat2_row = counts_df.filter(pl.col(cat_col) == cat2)
    
    cat1_hits = cat1_row["hit_counts"][0] if cat1_row.height > 0 else 0
    cat1_non_hits = cat1_row["non_hit_counts"][0] if cat1_row.height > 0 else 0
    cat2_hits = cat2_row["hit_counts"][0] if cat2_row.height > 0 else 0
    cat2_non_hits = cat2_row["non_hit_counts"][0] if cat2_row.height > 0 else 0
    
    # Create contingency table in format: [[cat1_hits, cat1_non_hits], [cat2_hits, cat2_non_hits]]
    table_for_test = np.array([
        [cat1_hits, cat1_non_hits],
        [cat2_hits, cat2_non_hits]
    ])
    
    # Perform Fisher's exact test
    try:
        odds_ratio, p_value = fisher_exact(table_for_test, alternative=alternative)
    except Exception as e:
        odds_ratio = np.nan
        p_value = np.nan
    
    # Create a readable contingency table DataFrame
    contingency_table = pl.DataFrame({
        hit_col: [True, False],
        cat1: [cat1_hits, cat1_non_hits],
        cat2: [cat2_hits, cat2_non_hits]
    })
    
    # Calculate hit rates for each category
    hit_rates = counts_df.with_columns(
        (pl.col("hit_counts") / pl.col("total_counts") * 100).alias("hit_perc")
    ).sort(cat_col)
    
    # Plotting functionality
    if plot:
        if ax is None:
            _, ax = plt.subplots(1, 1, figsize=(3,3.5))
        
        # Convert to pandas for plotting
        plot_df = hit_rates.to_pandas()
        
        # Create bar plot
        if palette is None:
            # Default colors
            colors = ['#CA7682', '#1D7AAB']  # Red and blue
        elif isinstance(palette, dict):
            colors = [palette.get(cat, '#808080') for cat in plot_df[cat_col]]
        else:
            colors = palette[:len(plot_df)]
        
        bars = ax.bar(
            range(len(plot_df)),
            plot_df["hit_perc"],
            color=colors,
            alpha=0.7,
            edgecolor='black',
            linewidth=1
        )
        
        # Set x-axis labels with category names and sample counts
        x_labels = [
            f"{cat}\n(n={n})" 
            for cat, n in zip(plot_df[cat_col], plot_df["total_counts"])
        ]
        ax.set_xticks(range(len(plot_df)))
        ax.set_xticklabels(x_labels, fontsize=11)
        
        # Add percentage labels on top of bars
        for i, (bar, perc) in enumerate(zip(bars, plot_df["hit_perc"])):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                perc + 1,  # Small offset above bar
                f"{perc:.1f}%",
                ha="center",
                va="bottom",
                fontweight="bold",
                fontsize=11
            )
        
        # Add comparison bracket and p-value annotation (similar to Annotator)
        if not np.isnan(p_value):
            # Format p-value for annotation
            def format_pvalue(pvalue):
                """Format p-value as 'p = xxx' or 'p < 0.0001' if below threshold"""
                if pvalue < 0.0001:
                    return "p < 0.0001"
                else:
                    return f"p = {pvalue:.4f}"
            
            formatted_p = format_pvalue(p_value)
            
            # Get bar positions and heights
            bar1_height = plot_df[plot_df[cat_col] == cat1]["hit_perc"].iloc[0]
            bar2_height = plot_df[plot_df[cat_col] == cat2]["hit_perc"].iloc[0]
            max_bar_height = max(bar1_height, bar2_height)
            
            # Calculate positions (bars are at x=0 and x=1)
            x1, x2 = 0, 1
            bracket_y = max_bar_height + 10  # Position bracket above bars
            text_y = max_bar_height + 12     # Position text above bracket
            
            # Draw comparison bracket
            # Horizontal line
            ax.plot([x1, x2], [bracket_y, bracket_y], 'k-', linewidth=1.5)
            # Left vertical line
            ax.plot([x1, x1], [bracket_y - 2, bracket_y], 'k-', linewidth=1.5)
            # Right vertical line
            ax.plot([x2, x2], [bracket_y - 2, bracket_y], 'k-', linewidth=1.5)
            
            # Add p-value text
            ax.text((x1 + x2) / 2, text_y, formatted_p,
                   ha='center', va='bottom', fontsize=12,
                   fontweight='bold', color='black')
        
        # Set title
        if title:
            plot_title = title
        else:
            plot_title = f"Hit Rate by {cat_col}"
        
        # Add p-value and odds ratio to title if available
        # if not np.isnan(p_value) and not np.isnan(odds_ratio):
        #     # Format p-value as 'p = xxx' or 'p < 0.0001' if below threshold
        #     if p_value < 0.0001:
        #         p_str = "p < 0.0001"
        #     else:
        #         p_str = f"p = {p_value:.4f}"
        #     plot_title += f"\nOR={odds_ratio:.2f}, {p_str}"
        
        ax.set_title(plot_title, fontsize=12, fontweight='bold')
        ax.set_ylabel("Hit Rate (%)", fontsize=11)
        max_perc = max(plot_df["hit_perc"]) if len(plot_df) > 0 and max(plot_df["hit_perc"]) > 0 else 10
        ax.set_ylim(0, 100)  # Add some space for labels
        ax.grid(alpha=0.2, axis='y')
        
        plt.tight_layout()
        plt.show()
    
    return {
        'contingency_table': contingency_table,
        'odds_ratio': odds_ratio,
        'p_value': p_value,
        'hit_rates': hit_rates,
        'table_array': table_for_test  # The numpy array used for the test
    }


def compute_aubprc(auprc, prior):
    return (auprc*(1-prior))/((auprc*(1-prior)) + ((1-auprc)*prior))


def plot_auroc_curves(df_clinvar_w_scores, label, pos_labels, methods, set2_palette=set2_palette):
    """
    For each method, plot the adjusted PR curve (with both balanced and monotonized adjustments)
    and compute the adjusted AUPRC via bootstrapping. Also perform pairwise significance testing 
    of the AUPRC across different RSA categories.
    
    Parameters
    ----------
    df_clinvar_w_scores : Polars DataFrame
        Input data including columns "clinvar_clnsig_clean", category and predictor scores.
    methods : list of str
        List of method names corresponding to score columns in the data.
    set2_palette : list of str
        List of colors (e.g., from a color palette) to use for different RSA categories.
    n_bootstrap : int, optional
        Number of bootstrap iterations for significance testing (default is 1000).
    random_state : int or None, optional
        Seed for reproducibility.
    
    Returns
    -------
    ap_results : dict
        Dictionary of adjusted AUPRC values and bootstrap distributions, structured as:
          { method: { category: {"AP": <adjusted AUPRC>, "boot_ap": <bootstrap array> } } }
    significance : dict
        Dictionary mapping each method to a DataFrame of pairwise p-values comparing AUPRC between categories.
    """
    fig, axes = plt.subplots(1, 1, figsize=(4, 4))
        
    ap_results = {}  # to store adjusted AUPRC and bootstrap distributions
    # significance = {}  # to store pairwise significance p-values per method
    
    for m_idx, met in enumerate(methods):
        ap_results[met] = {}
        # Filter data for current category and valid clinical labels.
        df_clin_per_method = df_clinvar_w_scores.with_columns(
            (
                pl.when(pl.col(label).is_in(pos_labels))
                .then(1)
                .otherwise(0)
            ).alias("clinvar_label")
        ).drop_nulls(subset=[label, met])
        
        # Convert to pandas DataFrame.
        df_pd = df_clin_per_method.to_pandas()
        df_pd = df_pd[["clinvar_label", met]].dropna()
        # Convert clinvar_label to boolean.
        truth = df_pd["clinvar_label"].astype(bool)
        # Create a scores DataFrame (with one column).
        scores = df_pd[[met]]
        # Create the YogiROC2 object using the new YogiROC2() constructor.
        fpr, tpr, thresholds = roc_curve(truth, scores)
        # Get the performance table from the YogiROC2 object.
        auc_score = roc_auc_score(truth, scores)
        ap_results[met] = auc_score

        # Plot the adjusted PR curve.
        axes.plot(100 * fpr, 100 * tpr,
                label=f"{met}: {auc_score:.2f}",
                color=set2_palette[methods.index(met)])
        # Plot chance level as a horizontal dashed line.
        # axes[m_idx].axhline(y=100 * prior, color=set2_palette[plddt_cats.index(cat)], ls="--")
    # Customize the subplot.
    handles, labels = axes.get_legend_handles_labels()
    # order = [2,1,0] if "rsa" in category else [3,1,0,2]
    precision_label = "True Positive Rate (%)"
    # title = f"AUROC by {category}"
    
    axes.legend(loc="lower right",
                    title="Assay")
    axes.set_xlabel("False Positive Rate (%)")
    axes.set_ylabel(precision_label)
    axes.set_title("AUROC")
    axes.grid(alpha=0.2)
    plt.tight_layout()
    plt.show()
    # plt.subplots_adjust(wspace=0.04)
    
    return ap_results ## , sig_results_dict


def plot_gene_level_summary(assay_df, assay, plot=False, ax=None, xlim=None, cat="clinvar_clnsig_clean", palette=None):
    gene_allele_cnts = assay_df.group_by(
        ["symbol"]
    ).len().rename(
        {"len": "total_counts_all"}
    ).sort("total_counts_all", descending=True)
    
    gene_allele_cat_cnts = assay_df.group_by(
        ["symbol", cat]
    ).len().rename(
        {"len": "total_counts"}
    )
    gene_allele_cat_hits = assay_df.filter(
        (pl.col(assay))
    ).group_by(
        ["symbol", cat]
    ).len().rename(
        {"len": "hit_counts"}
    )
    total_allele_hit_sum_df = gene_allele_cat_cnts.join(
        gene_allele_cat_hits,
        on=["symbol", cat],
        how="left"
    ).fill_null(0).with_columns(
        (pl.col("hit_counts") / pl.col("total_counts") * 100).alias("hit_perc")
    ).join(gene_allele_cnts, on="symbol").sort(
        "total_counts_all", "symbol", cat, 
        descending=[True, True, False]
    )

    if plot:
        genes_with_hits = total_allele_hit_sum_df.filter(pl.col("hit_perc")>0)["symbol"].unique()
        # Convert to pandas and sort
        total_allele_hit_sum_df = total_allele_hit_sum_df.to_pandas().reset_index(drop=True)
        n_genes = len(total_allele_hit_sum_df)
        
        # Add percentage labels
        # Get unique categories and their positions
        categories = total_allele_hit_sum_df[cat].unique()
        n_categories = len(categories)
        # Calculate bar positions for each category
        bar_height = 0.8 / n_categories  # Default seaborn bar height divided by number of categories
        
        # Dynamically scale figure size
        if ax is None:
            fig_height = n_genes * bar_height * n_categories / 2.5
            fig_width = 6  # Slightly wider to accommodate labels
            fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        if palette is None:
            palette = "Set2"
            
        # Create the barplots
        # Total counts (background bars)
        sns.barplot(
            data=total_allele_hit_sum_df, 
            y="symbol", 
            x="total_counts", 
            hue=cat,  # "clinvar_clnsig_clean"
            palette=palette,
            alpha=0.4,  # Make background bars slightly transparent
            ax=ax,
            # gap=.1,
            order=total_allele_hit_sum_df["symbol"].drop_duplicates(),
        )
        # Hit counts (foreground bars)
        sns.barplot(
            data=total_allele_hit_sum_df, 
            y="symbol", 
            x="hit_counts", 
            hue=cat,#"clinvar_clnsig_clean", 
            palette=palette,
            ax=ax,
            # gap=.1,
            order=total_allele_hit_sum_df["symbol"].drop_duplicates()
        )
        for i, gene in enumerate(total_allele_hit_sum_df["symbol"].unique()):
            gene_data = total_allele_hit_sum_df[total_allele_hit_sum_df["symbol"] == gene]
            for j, category in enumerate(categories):
                cat_data = gene_data[gene_data[cat] == category]
                if not cat_data.empty:
                    # Calculate bar position
                    bar_center = i + (j - (n_categories - 1) / 2) * bar_height
                    # Get the hit percentage value
                    hit_perc = cat_data["hit_perc"].iloc[0]
                    hit_counts = cat_data["hit_counts"].iloc[0]
                    # Add percentage label at the end of the hit_counts bar
                    ax.text(
                        hit_counts + max(total_allele_hit_sum_df["hit_counts"]) * 0.01,  # Small offset
                        bar_center,
                        f"{hit_perc:.1f}%",
                        ha="left",
                        va="center",
                        fontsize=9 if n_genes > 20 else 10,
                        fontweight="bold"
                    )
        
        # # Formatting
        ax.grid(alpha=0.2, axis='x')
        ax.set_xlabel("Count", fontsize=11)
        ax.set_ylabel("Gene Symbol", fontsize=11)
        ax.set_title(f"{assay}", fontsize=13)
        
        # Improve legend
        handles, labels = ax.get_legend_handles_labels()
        # Remove duplicate legend entries (seaborn creates duplicates with multiple barplots)
        unique_labels = []
        unique_handles = []
        for handle, label in zip(handles, labels):
            if label not in unique_labels:
                unique_labels.append(label)
                unique_handles.append(handle)
        
        ax.legend(unique_handles, unique_labels, 
                title=cat, 
                fontsize=10, 
                loc="lower right",
                framealpha=0.9)
        if xlim:
            ax.set_xlim(xlim)
        
        plt.tight_layout()
        plt.show()
    
    return total_allele_hit_sum_df