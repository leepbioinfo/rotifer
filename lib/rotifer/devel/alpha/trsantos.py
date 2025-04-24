def taxon_summary(
    df,
    level,
    column='classification',
    separator=';',
    unclassified_terms=[
        'unclassified',
        'environmental samples',
        'uncultured',
        'incertae sedis',
        'candidatus',
        'synthetic construct',
        'other']
):
    """
    Summarizes the taxonomic composition of a dataset at a specified taxonomic level.

    Parameters
    ----------
    df : pd.DataFrame
        Input table containing taxonomic lineage strings.
    level : int
        The taxonomic rank to extract (e.g., 1 = kingdom/domain, 2 = phylum, 3 = class, etc.).
    column : str, default 'classification'
        Name of the column in `df` that holds the full lineage strings.
    separator : str, default ';'
        Character or substring used to split the lineage into ranks.
    unclassified_terms : list of str, default ['unclassified','environmental samples','uncultured','incertae sedis','candidatus','synthetic construct','other']
        Keywords which, if found in the extracted term, force it to be labeled 'unclassified'.

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
    >>> summary = taxon_summary(rad_sam.ndf, 1, column='lineage', separator='>')
    >>> summary.head()
    """
    
    import pandas as pd
    import re

    # Pre-compile regex for faster checks (case-insensitive)
    uncl_pattern = re.compile(
        r'\b(?:' + '|'.join(re.escape(t) for t in unclassified_terms) + r')\b',
        flags=re.IGNORECASE
    )

    def extract_taxonomic_level(lineage):
        # Handle missing or empty
        if pd.isna(lineage) or not isinstance(lineage, str) or not lineage.strip():
            return 'unclassified'

        parts = lineage.split(separator)
        if len(parts) < level:
            return 'unclassified'

        term = parts[level - 1].strip().lower()
        # If the term itself is empty, or matches any placeholder keyword, classify as unclassified
        if not term or uncl_pattern.search(term):
            return 'unclassified'

        return term

    # Extract the specified taxonomic level
    taxonomic_levels = df[column].apply(extract_taxonomic_level)

    # Calculate absolute and normalized counts
    absolute_counts = taxonomic_levels.value_counts()
    normalized_counts = taxonomic_levels.value_counts(normalize=True) * 100

    # Combine results into a single DataFrame
    summary = pd.DataFrame({
        'Count': absolute_counts,
        'Percentage': normalized_counts
    })

    summary.index.name = 'Taxon'
    summary.reset_index(inplace=True)

    return summary

def shannon(self, ignore_gaps=True):
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
