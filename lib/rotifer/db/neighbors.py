#!/usr/bin/env python3

import argparse # Not our command line parser of choice but better than none
import os
import psycopg2
from psycopg2.extras import DictCursor
import sys
import rotifer
logger = rotifer.logging.getLogger(__name__)

class rneighbors:
    def __init__(self, above = 3, below = 3,
                 outformat = 'gi2operons', input_file ='',
                 gacc = [], asm = [],
                 database = 'rotifer', port = '5433',
                 host = '10.1.1.1', user = 'rotifer',
                 ):

        self.above = above # get above
        self.below = below # get below
        self.outformat = outformat # Output format
        self.input_file = input_file # Input file with Acc
        self.gacc = gacc # Use gacc
        self.asm = asm # use asm

        self.database = database # DB name
        self.port = port # DB port
        self.host = host # DB host
        self.user = user # DB user

    def connect(self):
        try:
            dsn = "dbname=" + self.database + " user=" + self.user + " host=" + self.host + ' port=' + self.port
            self.db  = psycopg2.connect(dsn, cursor_factory=DictCursor)

        except:
            logger.critical("I am unable to connect to the database")
            exit(1)

    def lsqlAsm(self):
        move = [
            '''
            INSERT INTO uasm(assembly) VALUES {0}
            '''.format(*(self.asm, )),
            '''
            INSERT INTO uquery (
                assembly, genomic_accession, feature,
                feature_class, product_accession, genomic_order
            )
            SELECT
                ft.assembly, ft.genomic_accession, ft.feature,
                ft.feature_class, ft.product_accession, ft.genomic_order
              FROM ftable AS ft
             INNER JOIN upid as p USING (product_accession)
             ORDER BY assembly, genomic_accession, genomic_order
            ''']
        return move

    def lsqlGacc(self):
        move = [
            '''
            insert into ugacc(genomic_accession) values {0}
            '''.format(*(self.gacc,)),
            '''
            insert into uquery (
                assembly, genomic_accession, feature,
                feature_class, product_accession, genomic_order
            )
            select
                ft.assembly, ft.genomic_accession, ft.feature,
                ft.feature_class, ft.product_accession, ft.genomic_order
              from ftable as ft
             inner join ugacc as a using (genomic_accession)
             inner join upid as p using (product_accession)
             order by assembly, genomic_accession, genomic_order
            '''
            ]
        return move

    def baseSql(self):
        self.ids=",".join(map(lambda x: str("('" + x + "')"), self.input_file))

        lsql = [
            '''
            CREATE TEMPORARY TABLE upid (product_accession varchar(64))
            ''',
            '''
            CREATE TEMPORARY TABLE uquery (
                assembly varchar(64), genomic_accession varchar(64),
                feature varchar(64), feature_class varchar(64),
                product_accession varchar(64), genomic_order integer
            );
            ''',
            '''
            INSERT INTO upid(product_accession) VALUES {0}
            '''.format(*(self.ids,))]
        return lsql

    def basicSql(self):
        move = [
            '''
            INSERT INTO uquery (
                assembly, genomic_accession, feature,
                feature_class, product_accession, genomic_order
            )
            SELECT
                ft.assembly, ft.genomic_accession, ft.feature,
                ft.feature_class, ft.product_accession, ft.genomic_order
              FROM ftable AS ft
             INNER JOIN upid as p USING (product_accession)
             ORDER BY assembly, genomic_accession, genomic_order
            ''']
        return move

    def post_ids(self):
        '''
        Send queris (like NCBI epost)
        Build SQL statements
        Decide if filter by ASM, GACC or nothing
        '''
        if self.asm:
            lsql = self.baseSql() + self.lsqlAsm()
        if self.gacc:
            lsql = self.baseSql() + self.lsqlGacc()
        else:
            lsql = self.baseSql() + self.basicSql()
        cursor = self.db.cursor()
        for sql in lsql:
            cursor.execute(sql)
        cursor.execute("SELECT count(*) from upid")
        logger.debug(f'## Number of protein accession numbers in input: {cursor.fetchone()[0]')
        cursor.execute('''SELECT count(distinct product_accession) as proteins from uquery''')
        logger.debug(f'## Number of proteins found in the SQL database: {cursor.fetchone()[0]}')

    def runRneighbors(self, stop = True):
        self.connect() # connect to DB
        self.post_ids() # post to DB

        cursor = self.db.cursor()

        # Some outformat
        if self.outformat == 'missing': self.missing(cursor)
        if self.outformat == 'query': self.query(cursor)
        # Formats for neighborhoods
        else:
            intervals = self.db.cursor()
            print_header = 1
            for nt in self.fetch_contigs():
                self.fetch_intervals(intervals, nt['genomic_accession'], self.above, self.below)
                old = intervals.fetchone()
                for new in intervals:
                    if not isinstance(new[4],(int)):
                        logger.info(
                                "Last upstream neighbor for " + old["product_accession"]
                                + " is not an integer: " + str(type(old[4])) + "\n"
                        )
                    if int(old[4])+1 >= int(new[3]):
                        old[4] = new[4] # merge old and new intervals
                    else: # print old, put new in old
                        if self.outformat == 'gi2operons':
                            self.gi2operons(cursor, nt, old)
                        elif self.outformat == 'table':
                            print_header = self.table(cursor, nt, old, print_header)
                        old = new
                if old != []:
                    if self.outformat == 'gi2operons':
                        self.gi2operons(cursor, nt, old)
                    elif self.outformat == 'table':
                        print_header = self.table(cursor, nt, old, print_header)
                    else:
                        logger.warning("Unsupported output format: " + self.outformat)
            intervals.close()
        cursor.close()
        self.db.close()
        if stop:
            sys.exit(0)
        else:
            pass

    def missing(self, cursor):
        lost = self.fetch_missing_queries(cursor)
        lost = list(map(lambda x: str(x[0]), lost))
        if len(lost):
            logger.warning('\n'.join(lost))

    def query(self, cursor):
        lost = self.fetch_queries(cursor)
        lost = list(map(lambda x: str(x[0], lost)))
        if len(lost):
            logger.warning('\n',join(lost))

    def gi2operons(self, cursor, contig, interval):
        # Prepare main query to fetch a block
        sql = '''
            SELECT
              CASE WHEN q.genomic_order IS NULL THEN '.' ELSE '-->' END AS is_query,
              ft.gstart || '..' || gend as cds,
              strand as dir,
              product_length as len,
              CASE WHEN ft.product_accession IS NULL THEN '.' ELSE ft.product_accession END as pid,
              CASE WHEN ft.feature = 'CDS' AND ft.product_accession IS NULL THEN 'PSE' ELSE ft.feature END as type,
              CASE WHEN ft.symbol IS NULL THEN '.' ELSE ft.symbol END as gene,
              CASE WHEN ft.locus_tag IS NULL THEN '.' ELSE ft.locus_tag END as locus,
              '.' as gi,
              CASE WHEN ft.product_name IS NULL THEN '.' ELSE ft.product_name END as product
            FROM ftable AS ft LEFT OUTER JOIN uquery AS q USING (genomic_accession, genomic_order)
            WHERE
                ft.genomic_accession = '{0}'
                AND ft.feature != 'gene'
                AND genomic_order BETWEEN {1} AND {2}
            ORDER BY genomic_order
        '''
        sql  = sql.format(*(contig['genomic_accession'], interval['gomin'], interval['gomax']))

        # Fetch column lengths
        sql2 = '''
            SELECT
                3 as is_query,
                max(length(cds)) as cds,
                3 as dir,
                CASE WHEN max(length(len::text)) > 3 THEN max(length(len::text)) ELSE 3 END as len,
                max(length(pid)) as pid,
                CASE WHEN max(length(type))  > 4 THEN max(length(type))  ELSE 4 END as type,
                CASE WHEN max(length(gene))  > 4 THEN max(length(gene))  ELSE 4 END as gene,
                CASE WHEN max(length(locus)) > 5 THEN max(length(locus)) ELSE 5 END as locus,
                2 as gi,
                CASE WHEN max(length(product)) > 7 THEN max(length(product)) ELSE 7 END as product
            FROM ({0}) as gi2op
        '''
        sql2 = sql2.format(*(sql,))
        cursor.execute(sql2)
        collen = cursor.fetchone() # collen is a psycopg2.extras.DictRow

        # Print header
        print("ORGANISM",contig['organism_name'],"accession no is",
              contig['genomic_accession'],"Protein is",interval['product_accession'])
        header = [".","cds","dir","len","pid","type","gene","locus","gi"]
        header = [ header[i].ljust(collen[i]) for i in range(0,len(header)) ]
        header.append('product')
        print("  ".join(header))

        # Fetch and format the block
        cursor.execute(sql)
        for row in cursor:
            last = len(row)-1
            row[last] = str(row[last])
            for col in range(0,last):
                if row[col] is None:
                    row[col] = '.'
                else:
                    row[col] = str(row[col])
                row[col] = row[col].ljust(collen[col])
            print("  ".join(row))
        print("---------------------------------------")

    def table(self, cursor, contig, interval, print_header):
        sql = '''
            SELECT
                ft.genomic_accession as nucleotide,
                ft.gstart as start, ft.gend as end,
                CASE WHEN ft.strand = '+' THEN 1 ELSE -1 END as strand,
                CASE WHEN ft.feature = 'CDS' AND ft.product_accession IS NULL THEN 'PSE' ELSE ft.feature END as type,
                '{0}' as query, ft.locus_tag as locus, ft.seq_type, ft.assembly,
                CASE WHEN ft.product_accession IS NULL THEN '' ELSE ft.product_accession END as pid,
                CASE WHEN ft.symbol IS NULL THEN '' ELSE ft.symbol END as gene,
                CASE WHEN ft.product_name      IS NULL THEN '' ELSE ft.product_name      END as product,
                a.organism_name, a.lineage as classification
            FROM ftable AS ft
                LEFT OUTER JOIN uquery AS q USING (genomic_accession, genomic_order)
                INNER JOIN assembly_summary AS a ON (ft.assembly = a.assembly)
            WHERE
                ft.genomic_accession = '{1}'
                AND ft.feature != 'gene'
                AND genomic_order BETWEEN {2} AND {3}
            ORDER BY genomic_order
        '''
        sql = sql.format(*(interval['product_accession'],contig['genomic_accession'],
                           interval['gomin'], interval['gomax']))
        cursor.execute(sql)
        for row in cursor:
            if print_header:
                print("\t".join(row.keys()))
                print_header = 0
            row = [ str(row[i] or '') for i in range(0,len(row)) ]
            print("\t".join(row))
        return print_header

    def fetch_intervals(self, cursor, ntid, above, below):
        args=(ntid, above, below)
        sql = '''
        SELECT genomic_accession, product_accession, genomic_order,
           first_upstream_neighbor(genomic_accession, genomic_order, {1}) as gomin,
           last_downstream_neighbor(genomic_accession, genomic_order, {2}) as gomax
        FROM uquery AS q
        WHERE q.genomic_accession = '{0}'
        ORDER BY genomic_accession, genomic_order, gomin, gomax
        '''
        cursor.execute(sql.format(*args))
        return cursor

    def fetch_queries(self,cursor):
        cursor.execute('''
            SELECT DISTINCT p.product_accession AS missing_queries
            FROM upid AS p LEFT OUTER JOIN uquery AS q USING (product_accession)
            WHERE q.product_accession IS NOT NULL;
        ''')
        return cursor.fetchall()

    def fetch_missing_queries(self,cursor):
        cursor.execute('''
            SELECT DISTINCT p.product_accession AS missing_queries
            FROM upid AS p LEFT OUTER JOIN uquery AS q USING (product_accession)
            WHERE q.product_accession IS NULL;
        ''')
        return cursor.fetchall()

    def fetch_contigs(self):
       cursor = self.db.cursor()
       sql = '''
       SELECT DISTINCT assembly, genomic_accession, organism_name, lineage
       FROM uquery INNER JOIN assembly_summary using (assembly)
       ORDER BY lineage, organism_name, assembly, genomic_accession
       '''
       cursor.execute(sql)
       return cursor

