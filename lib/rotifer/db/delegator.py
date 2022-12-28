# External libraries
import types
from copy import deepcopy

# Rotifer
import rotifer.db.core
import rotifer
logger = rotifer.logging.getLogger(__name__)

class DelegatorCursor(rotifer.db.core.BaseCursor):
    def __init__(self, readers=[], writers=[], progress=True, tries=3, batch_size=None, threads=10, *args, **kwargs):
        super().__init__(progress=progress, *args, **kwargs)
        self.readers = readers
        self.writers = writers
        self.tries = tries
        self.batch_size = batch_size
        self.threads = threads
        self.reset_cursors()

    @property
    def _cursor_modules(self):
        import inspect
        import importlib
        mymodule = inspect.getmodule(self)
        cursor_modules = dict()

        # Check configuration
        try:
            myconfig = getattr(mymodule,'config')
        except:
            error = f'Module {mymodule.__name__} has no configuration! Blame the developer!'
            logger.error(error)
            raise ValueError(f'No attribute "config" in module {mymodule.__name__}')
        if 'readers' not in myconfig:
            error = f'Missing dictionary of reader modules in configuration of module {mymodule.__name__}'
            logger.error(error)
            raise ValueError(error)
        if 'writers' not in myconfig:
            error = f'Missing dictionary of writer modules in configuration of module {mymodule.__name__}'
            logger.error(error)
            raise ValueError(error)

        # Load modules
        for module in set(self.readers + self.writers):
            if "readers" in myconfig and module in myconfig['readers']:
                module_name = myconfig['readers'][module]
            elif "writers" in myconfig and module in myconfig['writers']:
                module_name = myconfig['writers'][module]
            else:
                error = f'Missing module name "{module}" for writers or readers in {mymodule.__name__}.config'
                logger.error(error)
                raise ValueError(error)
            try:
                cursor_modules[module] = importlib.import_module(module_name)
            except:
                logger.error(f'Unable to load module {module_name}: %s.', exc_info=1)
                raise ImportError(f'Unable to load module {module_name}')

        return cursor_modules

    def reset_cursors(self):
        myname = str(type(self)).split("'")[1].split(".")[-1]
        if hasattr(self,"_shared_attributes"):
            kwargs = { x: getattr(self,x) for x in self._shared_attributes }
        else:
            kwargs = dict()
        self.cursors = dict()
        cursor_modules = self._cursor_modules
        for modulename in cursor_modules:
            module = cursor_modules[modulename]
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
    def __init__(self, readers=[], writers=[], progress=True, tries=3, batch_size=None, threads=10, *args, **kwargs):
        super().__init__(readers=readers, writers=writers, progress=progress, tries=tries, batch_size=batch_size, threads=threads, *args, **kwargs)

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
        for i in range(0,len(self.readers)):
            if len(todo) == 0:
                break
            cursorName = self.readers[i]
            if cursorName in self.cursors:
                cursor = self.cursors[cursorName]
            else:
                continue
            result = cursor.__getitem__(todo, *args, **kwargs)
            found = self.getids(result, *args, **kwargs)
            done = targets.intersection(found)
            for j in range(0,i+1):
                c = self.readers[j]
                if c not in self.cursors:
                    continue
                self.cursors[c].remove_missing(done)
            self.remove_missing(done)
            self.update_missing(data=cursor._missing)
            if isinstance(result,types.NoneType):
                continue
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
        for i in range(0,len(self.readers)):
            if len(todo) == 0:
                break
            cursorName = self.readers[i]
            if cursorName in self.cursors:
                cursor = self.cursors[cursorName]
            else:
                continue
            for result in cursor.fetchone(todo, *args, **kwargs):
                found = self.getids(result, *args, **kwargs)
                done = todo.intersection(found)
                for j in range(0,i+1):
                    c = self.readers[j]
                    if c not in self.cursors:
                        continue
                    self.cursors[c].remove_missing(done)
                for j in self.writers:
                    if j not in self.cursors:
                        continue
                    self.cursors[j].insert(result)
                self.remove_missing(done)
                self.update_missing(data=cursor._missing)
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

