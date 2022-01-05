from rotifer.genome.io import parse
from rotifer.genome.io import GenbankAssemblyFromFilename

@classmethod
def from_genbank(handle, informat, *args, assembly=GenbankAssemblyFromFilename, **kwargs):
    """
    Load genomic data from a list of genbank files.
    Input files may be compressed or downloaded with Bio.Entrez.

    Arguments
    ---------
        handle   : one or more file handles or filenames
        informat : either gff or a Bio.SeqIO supported format
        assembly : rule for setting genome assembly accession ID
                   Assembly identifiers are not part of the INSDC 
                   (http://www.insdc.org/) standard but can be
                   found in some NCBI entries under the DBLINK
                   tag, such as in many RefSeq entries:

                   DBLINK    BioProject: PRJNA224116
                             BioSample: SAMEA3138382
                             Assembly: GCF_000210855.2 <-- here!

                   This parameter allows the user to automate the definition
                   of genome identifiers. You can use either a string or a
                   lambda function that takes two arguments: file path and
                   sequence record object.

                   See GenbankAssemblyFromFilename and SequenceNameToAssembly
                   for examples of functions compatible with this parameter.

    """
    from rotifer.genome.utils import seqrecords_to_dataframe
