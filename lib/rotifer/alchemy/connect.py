#!/usr/bin/env python3
from sqlalchemy import create_engine, MetaData, Table
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker


class clickhouse:
    def __init__(self,
                 uri = 'clickhouse://default:@localhost:8123/rotifer',
                 table_name = ''
                 ):

        '''Connect using SQLAlchemy to a database'''

        engine = create_engine(uri)
        Session = sessionmaker(bind=engine)
        Base = declarative_base()
        metadata = MetaData(bind = engine)
        metadata.reflect()
        self.conn = engine.connect()
        self.table = metadata.tables[table_name]
        self.session = Session()

    def table_name(self):
        return (self.table.key)

    def table_columns(self):
        return (self.table.columns.keys())



# Statement
# s = table.select().where(
# table.c.sequence.in_(list here))

# Execute Statement
# conn.execute(s)


