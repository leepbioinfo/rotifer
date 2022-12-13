# External libraries
from copy import deepcopy

# Rotifer
import rotifer.db.core
import rotifer
logger = rotifer.logging.getLogger(__name__)

class DelegatorCursor(rotifer.db.core.BaseCursor):
    def __init__(self, methods=[], progress=True, tries=3, batch_size=None, threads=10, *args, **kwargs):
        super().__init__(progress=progress, *args, **kwargs)
        self.methods = methods
        self.tries = tries
        self.batch_size = batch_size
        self.threads = threads
        self.load_cursor_modules()
        self.reset_cursors()

    def load_cursor_modules(self):
        import inspect
        import importlib
        mymodule = inspect.getmodule(self)
        self._cursor_modules = dict()

        # Check configuration
        try:
            myconfig = getattr(mymodule,'config')
        except:
            error = f'Module {mymodule.__name__} has no configuration! Blame the developer!'
            logger.error(error)
            raise ValueError(f'No attribute "config" in module {mymodule.__name__}')
        if 'cursor_methods' not in myconfig:
            error = f'Missing dictionary of methods in configuration of module {mymodule.__name__}'
            logger.error(error)
            raise ValueError(error)

        # Load modules
        for module in self.methods:
            if module not in myconfig['cursor_methods']:
                error = f'Method {module} not found in {mymodule.__name__}.config["cursor_methods"]'
                logger.error(error)
                raise ValueError(error)
            module_name = myconfig['cursor_methods'][module]
            try:
                self._cursor_modules[module] = importlib.import_module(module_name)
            except:
                logger.error(f'Unable to load module {module_name}: %s.', exc_info=1)
                raise ImportError(f'Unable to load module {module_name}')

    def reset_cursors(self):
        myname = str(type(self)).split("'")[1].split(".")[-1]
        if hasattr(self,"_shared_attributes"):
            kwargs = { x: getattr(self,x) for x in self._shared_attributes }
        else:
            kwargs = dict()
        self.cursors = dict()
        for modulename in self._cursor_modules:
            module = self._cursor_modules[modulename]
            try:
                cursorClass = getattr(module,myname)
            except:
                logger.error(f'Module {module.__name__} does not define a {myname} class')
                continue
            self.cursors[modulename] = cursorClass(**kwargs)

    def __setattr__(self, name, value):
        super().__setattr__(name, value)
        if hasattr(self,'cursors') and hasattr(self,'_shared_attributes') and name in self._shared_attributes:
            for cursor in self.cursors.values():
                if hasattr(cursor,name):
                    cursor.__setattr__(name,value)

class SequentialDelegatorCursor(DelegatorCursor):
    def __init__(self, methods=[], progress=True, tries=3, batch_size=None, threads=10, *args, **kwargs):
        super().__init__(methods=methods, progress=progress, tries=tries, batch_size=batch_size, threads=threads, *args, **kwargs)

    def __getitem__(self, accessions, *args, **kwargs):
        """
        Dictionary-like access to data.

        Usage
        -----
        General template of the __getitem__ method,
        for any BaseCursor child class:

        >>> import rotifer.db.ncbi as ncbi
        >>> tc = ncbi.TaxonomyCursor(progress=True)
        >>> t = tc[2599]

        Parameters
        ----------
        accessions: list of strings
          Database identifiers.

        Returns
        -------
        Same value as the delegated cursors
        See documentation for the cursors attribute
        """
        targets = self.parse_ids(accessions)

        # Call cursors
        data = []
        todo = deepcopy(targets)
        for i in range(0,len(self.methods)):
            if len(todo) == 0:
                break
            cursorName = self.methods[i]
            if cursorName in self.cursors:
                cursor = self.cursors[cursorName]
            else:
                continue
            result = cursor.__getitem__(todo, *args, **kwargs)
            found = self.getids(result, *args, **kwargs)
            done = targets.intersection(found)
            for j in range(0,i+1):
                c = self.methods[j]
                if c not in self.cursors:
                    continue
                self.cursors[c].remove_missing(done)
            self.remove_missing(done)
            for missing in cursor.missing.iterrows():
                self.update_missing(missing[0], *missing[1].drop("class").to_list())
            if isinstance(result, list):
                data.extend(result)
            else:
                data.append(result)
            todo = todo - done
        if len(targets) == 1 and len(data) == 1:
            data = data[0]
        return data

    def fetchone(self, accessions, *args, **kwargs):
        """
        Get a generator to fetch data for iteratively.

        Parameters
        ----------
        accessions: list of strings
          Database identifiers.

        Returns
        -------
        Generator of Pandas dataframes
        """
        targets = self.parse_ids(accessions)

        # Call cursors
        todo = deepcopy(targets)
        for i in range(0,len(self.methods)):
            if len(todo) == 0:
                break
            cursorName = self.methods[i]
            if cursorName in self.cursors:
                cursor = self.cursors[cursorName]
            else:
                continue
            for result in cursor.fetchone(todo, *args, **kwargs):
                found = self.getids(result, *args, **kwargs)
                done = todo.intersection(found)
                for j in range(0,i+1):
                    c = self.methods[j]
                    if c not in self.cursors:
                        continue
                    self.cursors[c].remove_missing(done)
                self.remove_missing(done)
                for missing in cursor.missing.iterrows():
                    self.update_missing(missing[0], *missing[1].drop("class").to_list())
                todo = todo - done
                yield result

    def fetchall(self, accessions, *args, **kwargs):
        """
        Fetch data for all accessions.

        Parameters
        ----------
        accessions: list of database identifiers
          Database identifiers.

        Returns
        -------
        Pandas dataframe
        """
        targets = self.parse_ids(accessions)
        stack = []
        for data in self.fetchone(targets, *args, **kwargs):
            if isinstance(data,list):
                stack.extend(data)
            else:
                stack.append(data)
        return stack

