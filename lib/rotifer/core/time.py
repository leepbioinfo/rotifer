#!/usr/bin/env python3

from time import time
from functools import wraps

def rtime(func):
    """
    :param func: Decorated function
    :return: Execution time for the decorated function
    """

    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time()
        result = func(*args, **kwargs)
        end = time()
        print(f'{func.__name__} executed in {end - start:.4f} s')
        return result
    return wrapper


if __name__ == '__main__':
    @rtime
    def sum(a,b):
        return a+b

    sum(1,2)

