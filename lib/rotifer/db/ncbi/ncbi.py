# Dependencies
import sys
import pandas as pd

# Access configuration
from rotifer.db.ncbi import NcbiConfig

# Import submodules used by _DispatchMethod
import rotifer.db.ncbi.fetch as fetcher
import rotifer.db.ncbi.parse as parser
import rotifer.db.ncbi.read as reader
_MAP = {
    'fetch' : fetcher,
    'parse' : parser,
    'read'  : reader
    }

class ncbi:
    def __init__(self, query=None, email=NcbiConfig['email'], api_key=NcbiConfig['api_key']):
        '''
        Search, download and parse NCBI data.

        Usage:
          a = ncbi(query=['WP_013925940.1','AMX08087.1','WP_009227829.1'])
          b = a.read(method=['blastdbcmd','entrez'], format='fasta')

        Parameters:
          query   : list of NCBI accessions, UIDs, terms, etc.
          email   : e-mail address to use for opening connections
          api_key : NCBI API key

        Configuration: ~/.rotifer/db/ncbi.yml

        The following parameters may be loaded from the YAML
        configuration file: email, api_key
        '''
        self.email        = email
        self.files        = [] # List of files or file objects
        self.__api_key    = api_key
        self.__connection = {} # Open connections
        self.__missing    = [] # List of missing targets
        self.__query      = [] # List of targets
        if query:
            self.submit(query)

    def reset(self):
        '''
        Refresh the object by erasing the list of queries and related data.

        Returns:
          The list of queries, if previously set
        '''
        self.files        = [] # List of files or file objects
        self.__connection = {} # Open connections
        self.__missing    = [] # List of missing targets
        query = self.__query if self.__query else None
        self.__query      = [] # List of targets
        return query

    def api_key(self, value=None):
        '''
        Retrieve first available NCBI API key.

        Search path is:
          - Environment variable NCBI_API_KEY
          - User configuration at ~/.rotifer/etc/db/ncbi.yml
        
        The api_key entry in the ncbi.yaml should look like

        api_key: ncbi_api_key
        '''
        if value:
            self.__api_key = value
            return value
        elif self.__api_key:
            return self.__api_key
        elif 'NCBI_API_KEY' in os.environ and os.environ['NCBI_API_KEY'] != '':
            return os.environ['NCBI_API_KEY']
        elif 'api_key' in NcbiConfig:
            if isinstance(NcbiConfig['api_key'],dict):
                if os.environ['USER'] in NcbiConfig['api_key']:
                    return NcbiConfig['api_key'][os.environ['USER']]
                else:
                    for k,v in NcbiConfig['api_key'].items():
                        return v
            else:
                return NcbiConfig['api_key']
        else:
            print("No NCBI API key found! Set environment variable NCBI_API_KEY or add your key to ~/.rotifer/db/ncbi.yml.")
            exit(1)

    def submit(self, target=None, replace=True):
        '''
        Register, update or get the list of queries

        Parameters:
          target  : list of queries
          replace : reset list if True, extend if False
        '''
        if target == None:
            return self.__query
        elif isinstance(target,set):
            target = list(target)
        elif not isinstance(target,list):
            target = [ target ]
        if replace:
            self.__query   = target
            self.__missing = [] # List of missing targets
        else:
            self.__query.extend(target)

    def missing(self, target=None, replace=True):
        '''
        Get or set list of missing queries.
        The list of missing targets is usually update internally.

        Parameters:
          target  : list of missing queries
          replace : reset list if True, extend if False
        '''
        if target == None:
            return self.__missing
        elif isinstance(target,set):
            target = list(target)
        elif not isinstance(target,list):
            target = [ target ]
        if replace:
            self.__missing = target
        else:
            self.__missing.extend(target)

    def read(self, method=[], concat=True, *args, **kwargs):
        '''
        Generic function to access or retrieve NCBI data and load
        all data to memory data structures. Data may be downloaded
        or read from a local source, depending on the method.

        The type of data structure returned depends on the query and
        the (optional) access and/or parsing method.

        Usage:

            from rotifer.db.ncbi import ncbi
            a = ncbi()
            b = a.read(method='assembly_reports') # b is a pandas DataFrame
            a.submit(["AKT35709.1","NP_250592.1"])
            c = a.read(method='entrez', db='protein') # c is a list

        Returns:
          Depends on the method used (see rotifer.db.ncbi.read)

        Parameters:
          method : (list of) access methods.
          concat : concatenate results from different methods
                   If concat == False, read() always returns a list.
          
          Additional parameters are defined by the access method.
          See rotifer.db.ncbi.read for details on each method.

        Available methods (see rotifer.db.ncbi.read):
         * assembly_reports
             Load asssembly_summary_*.txt files from directory or FTP
         * entrez
             Download data from NCBI using Bio.Entrez.
        '''
        return self.__MethodDispatcher(_name='read', method=method, concat=concat, *args, **kwargs)

    def parse(self, method=[], concat=True, *args, **kwargs):
        '''
        Generic function to apply parsers to local or remote NCBI data.

        This routine is similar to the read routine because it retrieves
        data from its source using fetch() methods and stores it locally
        but it saves memory by returning some sort of iterator that
        will only load the data into memory when requested.

        The type of data structure returned depends on the query and
        the (optional) access and/or parsing method.

        Usage:

            from rotifer.db.ncbi import ncbi
            a = ncbi()
            a.submit(['GCF_900504695.1', 'GCF_004636045.1', 'GCF_902726645.1'])
            c = a.parse(method='genomes') # c is a generator / iterator

        Returns:
          Depends on the method used (see rotifer.db.ncbi.read)

        Parameters:
          method : (list of) access methods.
          concat : concatenate results from different methods
                   If concat == False, read() always returns a list.
          
          Additional parameters are defined by the access method.
          See rotifer.db.ncbi.read for details on each method.

        Available methods (see rotifer.db.ncbi.read):
         * genomes
             Retrieve parsed genomes from NCBI
         * blastdbcmd
             Iterate over a list of Bio.SeqRecords from a local BLAST database.
         * entrez
             Download data from NCBI using Bio.Entrez.
        '''
        return self.__MethodDispatcher(_name='parse', method=method, concat=concat, *args, **kwargs)

    def fetch(self, method=[], concat=False, *args, **kwargs):
        '''
        Generic function to download data from NCBI.

        This method can decompress known formats automatically but
        doesn't make any attempt to parse the data.

        Usage:

            from rotifer.db.ncbi import ncbi
            a = ncbi()
            b = a.fetch(method='assembly_reports') # b is a file path
            a.submit(["AKT35709.1","NP_250592.1"])
            c = a.fetch(method='entrez', db='protein', format='fasta')

        Returns:
          List of fully qualified file paths or open filehandles

        Parameters:
          method : (list of) access methods.
          concat : if possible, generate a single iterator for all data

        Available methods (see rotifer.db.ncbi.fetch):
         * entrez
             Download data from NCBI using Bio.Entrez.
         * ftp
             Fetch data from a local BLAST database.
        '''
        return self.__MethodDispatcher(_name='fetch', method=method, concat=concat, *args, **kwargs)

    def __len__(self):
        '''
        Number of queries.
        '''
        return len(self.submit())

    def __MethodDispatcher(self, _name=None, method=[], *args, **kwargs):
        '''
        Internal method that selects access methods for each calller routine.
        '''

        # Applying each method
        if not method:
            print(f'No methods defined!', file=sys.stderr)
            return None
        elif type(method) != list:
            method = [ method ]
        results = []
        firstType = None
        concat = 1
        for m in method:
            selectedMethod = getattr(_MAP[_name], m, None)
            if callable(selectedMethod):
                if 'verbose' in kwargs and kwargs['verbose'] == True:
                    print(f'{__name__}: Dispatching method {m}...', file=sys.stderr)
                results.append(selectedMethod(self, *args, **kwargs))
                if 'verbose' in kwargs and kwargs['verbose'] == True:
                    print(f'{__name__}: Method {m} executed!', file=sys.stderr)
                curtype = type(results[-1])
                firstType = type(results[0])
                concat = concat and (curtype == firstType) and not getattr(selectedMethod,'_never_concatenate',False)
            else:
                print(f'Unknown error for method {m} in {_name} {_MAP[_name]}', file=sys.stderr)

        # Flatten results from different methods
        if len(results) == 1:
            results = results[0]
        elif firstType == pd.DataFrame:
            results = pd.concat(results, ignore_index=True, axis=0)
        elif firstType == list:
            results = [ x for y in results for x in y ]
        elif firstType == dict:
            results = { k: v for x in results for k,v in x.items() }

        # Return
        return results

# END
if __name__ == '__main__':
    pass
