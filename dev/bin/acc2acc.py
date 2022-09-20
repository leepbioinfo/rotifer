#!/usr/bin/env python3

import rotifer.core.cli as corecli
import rotifer.core.functions as rcf
import numpy as np
import pandas as pd
__version__ = 0.001
__authors__ = 'Gilberto Kaihami; Rodolfo Alvarenga'

def parse_cli():
    parser = corecli.parser(description = 'Convert an accession to another')

    parser.add(':cli.acc')

    # Add another options here

    parser.add( long_arg = '--target',
                short_arg = '-t',
                dest = 'target',
                arg_type = str,
                helper = 'Select target (e.g. GI)',
                action = "append"
                )

    args = parser.parse_args()

    return args

def connect(table_name):
    from rotifer.alchemy.connect import clickhouse

    uri = 'clickhouse://default:@localhost/rotifer'
    db = clickhouse(uri = uri,
                    table_name = table_name)
    conn = db.conn
    return conn

def types():

    dc= {{
        'UniProt':
            'UniProtKB': 'AC/ID' ACC+ID
    UniProtKB AC    ACC
    UniProtKB ID    ID
    UniParc UPARC
    UniRef50    NF50
    UniRef90    NF90
    UniRef100   NF100
    Gene name   GENENAME

    Category: Sequence databases
    EMBL/GenBank/DDBJ   EMBL_ID
    EMBL/GenBank/DDBJ CDS   EMBL
    Entrez Gene (GeneID)    P_ENTREZGENEID
    GI number   P_GI
    PIR PIR
    RefSeq Nucleotide   REFSEQ_NT_ID
    RefSeq Protein  P_REFSEQ_AC
    UniGene UNIGENE_ID
    Category: 3D structure databases
    PDB PDB_ID
    DisProt DISPROT_ID

    Category: Protein-protein interaction databases
    BioGrid BIOGRID_ID
    DIP DIP_ID
    MINT    MINT_ID
    STRING  STRING_ID
    Category: Chemistry
    ChEMBL  CHEMBL_ID
    DrugBank    DRUGBANK_ID
    GuidetoPHARMACOLOGY GUIDETOPHARMACOLOGY_ID
    SwissLipids SWISSLIPIDS_ID
    Category: Protein family/group databases
    Allergome   ALLERGOME_ID
    ESTHER  ESTHER_ID
    MEROPS  MEROPS_ID
    mycoCLAP    MYCOCLAP_ID
    PeroxiBase  PEROXIBASE_ID
    REBASE  REBASE_ID
    TCDB    TCDB_ID

    Category: Polymorphism and mutation databases
    BioMuta BIOMUTA_ID
    DMDM    DMDM_ID

    Category: 2D gel databases
    World-2DPAGE    WORLD_2DPAGE_ID

    Category: Protocols and materials databases
    DNASU   DNASU_ID

    Category: Genome annotation databases
    Ensembl ENSEMBL_ID
    Ensembl Protein ENSEMBL_PRO_ID
    Ensembl Transcript  ENSEMBL_TRS_ID
    Ensembl Genomes ENSEMBLGENOME_ID
    Ensembl Genomes Protein ENSEMBLGENOME_PRO_ID
    Ensembl Genomes Transcript  ENSEMBLGENOME_TRS_ID
    GeneDB  GENEDB_ID
    GeneID (Entrez Gene)    P_ENTREZGENEID
    KEGG    KEGG_ID
    PATRIC  PATRIC_ID
    UCSC    UCSC_ID
    VectorBase  VECTORBASE_ID
    WBParaSite  WBPARASITE_ID

    Category: Organism-specific databases
    ArachnoServer   ARACHNOSERVER_ID
    Araport ARAPORT_ID
    CCDS    CCDS_ID
    CGD CGD
    ConoServer  CONOSERVER_ID
    dictyBase   DICTYBASE_ID
    EchoBASE    ECHOBASE_ID
    EcoGene ECOGENE_ID
    euHCVdb EUHCVDB_ID
    EuPathDB    EUPATHDB_ID
    FlyBase FLYBASE_ID
    GeneCards   GENECARDS_ID
    GeneReviews GENEREVIEWS_ID
    H-InvDB H_INVDB_ID
    HGNC    HGNC_ID
    HPA HPA_ID
    LegioList   LEGIOLIST_ID
    Leproma LEPROMA_ID
    MaizeGDB    MAIZEGDB_ID
    MGI MGI_ID
    MIM MIM_ID
    neXtProt    NEXTPROT_ID
    Orphanet    ORPHANET_ID
    PharmGKB    PHARMGKB_ID
    PomBase POMBASE_ID
    PseudoCAP   PSEUDOCAP_ID
    RGD RGD_ID
    SGD SGD_ID
    TubercuList TUBERCULIST_ID
    WormBase    WORMBASE_ID
    WormBase Protein    WORMBASE_PRO_ID
    WormBase Transcript WORMBASE_TRS_ID
    Xenbase XENBASE_ID
    ZFIN    ZFIN_ID

    Category: Phylogenomic databases
    eggNOG  EGGNOG_ID
    GeneTree    GENETREE_ID
    HOGENOM HOGENOM_ID
    HOVERGEN    HOVERGEN_ID
    KO  KO_ID
    OMA OMA_ID
    OrthoDB ORTHODB_ID
    TreeFam TREEFAM_ID

    Category: Enzyme and pathway databases
    BioCyc  BIOCYC_ID
    Reactome    REACTOME_ID
    UniPathway  UNIPATHWAY_ID

    Category: Gene expression data
    CollecTF    COLLECTF_ID

    Category: Other
    ChiTaRS CHITARS_ID
    GeneWiki    GENEWIKI_ID
    GenomeRNAi  GENOMERNAI_ID

def find_in_uniprotdb(accs, source = '', target = ''):
    '''
    accs = list of accessions
    source = input types
    target = output types
    '''

    conn_source = connect('uniprotdb_id')
    conn_target = connect('uniprotdb_acc')

    if not source:
        source = []


if __name__ == '__main__':
    args = parse_cli()

    accs = args.accession

    split_size = 2000

    conn, cols = connect()

    targets = args.target

    if not targets:
        targets = ['P_REFSEQ_AC']

    filters = ",".join(["'" + x + "'" for x in targets])

    acc_splited = [accs[x:x+split_size] for x in range(0, len(accs), split_size) ]

    for acc in acc_splited:

        formated = ",".join(["'" + x + "'" for x in acc])
        first = conn.execute(f'''
            SELECT uniprotid, id from uniprotdb where id in ({formated})
''').fetchall()
        print(first)
        formated2 = ",".join(["'" + x + "'" for x,y in first])
        print(formated2)
        second = f'''
SELECT * from uniprotdb prewhere uniprotid in ({formated2})
'''
#         second = f'''
# SELECT S1.uniprotid, S2.id_type, S2.id
# (SELECT uniprotid from uniprotdb where id in ({formated})) S1
# LEFT JOIN (select * from uniprotdb) S2
# using (uniprotid)
# '''
        print(second)
        print(conn.execute(second).fetchall())
        # df = pd.DataFrame(second.fetchall(), columns = ['uniprotid', 'id_type', 'id']
        #                   )
        #
        # print(df)


#         # S1.id query
#         # S2.
#         selected = f'''
# select S1.id, S2.id_type, S2.id from
#         (select * from uniprotdb
#             where uniprotid in (select uniprotid from (select * from uniprotdb where id in ({formated})
#                                                       )
#                                )
#         ) S1
#         left join (
#                    select * from uniprotdb
#                             where uniprotid in (select uniprotid from (select * from uniprotdb where id in ({formated})
#                                                 )
#                                 AND
#                                     id_type in ({filters})
#                   ) S2
#         using (uniprotid)
#
#         '''
#
#         print(selected)
        # results = conn.execute(selected).fetchall()
        # results = np.unique(results, axis=0)
        # results = conn.execute(f"""SELECT DISTINCT id, id_type from uniprotdb where uniprotid in ({second}) and id_type in ({filters}) """)
        # for res in results:
        #     print(f'{res[0]}\t{res[1]}\t{res[2]}')


