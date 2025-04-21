def summarize_taxonomic_levels(df, column, level):
'''
Howdy there!!

This function aim is to create a taxonomic summary table given a dataframe (normally a ndf)
df asks for the dataframe, column askes for the column with taxonomic info and level asks for the taxonomic level (1, 2, 3...)

I do recommend using epsoares update_lineage function before using this one

Example usage:

rad_sam_tax = summarize_taxonomic_levels(rad_sam.ndf, 'lineage', 1)
rad_sam_tax.head() will output something like:

          Count  Percentage
bacteria  11350   93.162604
archaea     658    5.400969
            175    1.436428
'''
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
