def summarize_taxonomic_levels(df, column, level):
    """
    Summarizes the taxonomic composition of a dataset at a specified taxonomic level.

    Parameters
    ----------
    df : pd.DataFrame
        The input DataFrame, typically a normalized data format (ndf) containing taxonomic lineages.
    column : str
        The name of the column in the DataFrame that contains taxonomic lineage strings.
    level : int
        The taxonomic level to summarize (1 for kingdom/domain, 2 for phylum, etc.).
        Taxonomic levels are assumed to be separated by '>'.

    Returns
    -------
    pd.DataFrame
        A summary DataFrame with the count and percentage of entries at the specified taxonomic level.
        Example output:

                    Count  Percentage
        bacteria     11350   93.16
        archaea        658    5.40
        unclassified   175    1.44

    Notes
    -----
    It is recommended to run `epsoares.update_lineage` on your dataset before using this function.

    Example
    -------
    >>> summary = summarize_taxonomic_levels(rad_sam.ndf, 'lineage', 1)
    >>> summary.head()
    """
    import pandas as pd

    def extract_taxonomic_level(lineage, level):
        if pd.isna(lineage):  # Handle NaN values
            return 'unclassified'
        parts = lineage.split('>')
        if len(parts) >= level:
            return parts[level - 1].strip().lower()  # Standardize to lowercase for consistency
        else:
            return 'unclassified'

    # Extract the specified taxonomic level
    taxonomic_levels = df[column].apply(lambda x: extract_taxonomic_level(x, level))

    # Calculate absolute and normalized counts
    absolute_counts = taxonomic_levels.value_counts()
    normalized_counts = taxonomic_levels.value_counts(normalize=True) * 100

    # Combine results into a single DataFrame
    summary = pd.DataFrame({
        'Count': absolute_counts,
        'Percentage': normalized_counts
    })

    return summary

def compute_shannon_entropy(self, ignore_gaps=True):
    """
    Computes Shannon entropy for each column in the alignment.

    Parameters
    ----------
    ignore_gaps : bool
        Whether to exclude gaps ('-' or '.') from entropy calculation.

    Notes
    -----
    You need to run the following commands in order for this method to work properly:

    from rotifer.devel.beta import sequence as rdbs
    rdbs.sequence.compute_shannon_entropy = compute_shannon_entropy

    Returns
    -------
    entropy : pd.Series
        A pandas Series with entropy values indexed by column number.
    """

    from collections import Counter
    import numpy as np
    import pandas as pd


    # Make sure we only work with sequences
    sequences = self.df[self.df['type'] == 'sequence']['sequence'].tolist()

    if not sequences:
        raise ValueError("No sequences found to compute entropy.")

    # Transpose the alignment into columns
    alignment_length = len(sequences[0])
    columns = zip(*sequences)  # Converts list of sequences to columns

    entropy_values = []
    for col in columns:
        if ignore_gaps:
            col = [res for res in col if res not in ['-', '.']]
        freqs = Counter(col)
        total = sum(freqs.values())
        if total == 0:
            entropy = 0
        else:
            probs = [count / total for count in freqs.values()]
            entropy = -sum(p * np.log2(p) for p in probs if p > 0)
        entropy_values.append(entropy)

    return pd.Series(entropy_values, name='shannon_entropy')
